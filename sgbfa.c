/*
 * sgbfa.c
 *
 * This routine is the C version of the FORTRAN routine used in LINPACK
 * but converted to C by f2c. The purpose of this conversion is to allow
 * significant changes to be made to the code while maintaining it's
 * API to the other FORTRAN routines in which it might be linked.
 *
 * The goal of this is to make these routines capable of using SMP machines
 * so that a more powerful machine can solve the problem faster. Whether
 * or not this is reached, we do not know, but that's the goal.
 *
 * The first pass at this routine is to make use of threads and split up
 * the problem into several threads so that the operating system can
 * put the different threads on different processors.
 */

/* 	System Headers 			*/

/* 	Third-Party Headers 	*/

/* 	Other Headers 			*/
#include "f2c.h"

/* 	Forward Declarations 	*/
/*
 * These aren't 'forward', but in other files in the distribution and
 * need to be included here so that the compiler knows what to expect.
 */
extern int sscal_();
extern int saxpy_();

/* 	Public Data Types 		*/

/* 	Public Data Constants 	*/
/*
 * This is needed because the functions pass by reference, and some
 * times we need to pass in a constant 1, but need a reference to
 * it. That's why this is here.
 */
const integer one = 1;


/* ----------------------------------------------------------------------------
 *
 *     sgbfa factors a real band matrix by elimination
 *
 *     sgbfa is usually called by sgbco, but it can be called
 *     directly with a saving in time if rcond is not needed.
 *
 *     on entry
 *
 *        abd     real(lda, n)
 *                contains the matrix in band storage. The columns
 *                of the matrix are stored in the columns of abd and
 *                the diagonals of the matrix are stored in rows
 *                ml+1 through 2*ml+mu+1 of abd.
 *                see the comments below for details.
 *
 *        lda     integer
 *                The leading dimension of the array abd.
 *                lda must be .ge. 2*ml + mu + 1
 *
 *        n       integer
 *                the order of the original matrix.
 *
 *        ml      integer
 *                number of diagonals below the main diagonal.
 *                0 .le. ml .lt. n
 *
 *        mu      integer
 *                number of diagonals above the main diagonal.
 *                0 .le. mu .lt. n
 *                more efficient if ml .le. mu
 *     on return
 *
 *        abd     an upper triangular matrix in band storage and
 *                the multipliers which were used to obtain it.
 *                the factorization can be written a = l*u where
 *                l is a product of permutation and unit lower
 *                triangular matricies and u is upper triangular.
 *
 *        ipvt    integer(n)
 *                in integer vactor of pivot indices.
 *
 *        info    integer
 *                = 0  normal value
 *                = k  if u(k,k) .eq. 0.0 - this is not an error
 *                condition for this subroutine, but it does
 *                indicate that sgbsl will divide by zero if
 *                called. Use rcond in sgbco for a reliable
 *                indication of singularity.
 *
 *     band storage
 *
 *           if 'a' is a band matrix, the following program segment
 *           will set up the input.
 *
 *                   ml = (band width below the diagonal)
 *                   mu = (band width above the diagonal)
 *                   m = ml + mu + 1
 *                   do 20 j = 1, n
 *                     i1 = max0(1, j-mu)
 *                     i2 = min0(n, j+ml)
 *                     do 10 i = i1, i2
 *                       k = i - j + m
 *                       abd(k,j) = a(i,j)
 *                10   continue
 *                20 continue
 *
 *           this uses rows ml+1 through 2*ml+mu+1 of abd.
 *           in addition, the first ml rows in abd are used for
 *           elements generated during triangularization.
 *           the total number of rows needed in abd is 2*ml+mu+1.
 *           the ml+mu by ml+mu upper left triangle and the
 *           ml by ml lower right triangle are not referenced.
 *
 *     linpack. this version dated 08/14/78
 *     cleve moler, university of new mexico, argonne national lab
 *
 *     subroutines and functions
 *
 *     blas saxpy, sscal, isamax
 *     fortran max0, min0
 *
 *     internal variables
 *
 *               don't interchange rows for this code
 *     integer isamax
 */

int sgbfa_(doublereal *abd, integer *lda, integer *n, integer *ml, integer *mu, integer *ipvt, integer *info) {
	// do these for speed
	integer sld = *ml;
	integer sud = *mu;
	integer m = sld + sud + 1;
	integer sfda = *lda;

	// Local variables
	integer i, j, k, l;
	doublereal t;
	integer i0, j0, j1, imin, imax;
	integer lm, mm, ju, jz, kp1, nm1;

	/* Parameter adjustments */
	abd -= 1 + sfda * 1;
	--ipvt;

	/* Function Body */
	*info = 0;

	/*
	 * zero initial fill-in columns
	 */
	j0 = sud + 2;
	j1 = min(*n, m) - 1;
	if (j1 >= j0) {
		for (jz = j0; jz <= j1; jz++) {
			imin = m + 1 - jz + jz * sfda;
			imax = sld + jz * sfda;
			for (i = imin; i <= imax; i++) {
				abd[i] = 0.0;
			}
		}
	}
    jz = j1;
    ju = 0;

	/*
	 * gaussian elimination with partial pivoting
	 */
	nm1 = *n - 1;
	if (nm1 >= 1) {
		for (k = 1; k <= nm1; k++) {
			kp1 = k + 1;

			/*
			 * zero next fill-in column
			 */
			++jz;
			if ((sld >= 1) || (jz <= *n)) {
				imin = 1 + jz * sfda;
				imax = imin - 1 + sld;
				for (i = imin; i <= imax; i++) {
					abd[i] = 0.0;
				}
			}

			/*
			 * find l = pivot index
			 */
			// computing MIN
			lm = min(sld, (*n - k));
//			l = isamax(lm+1,abd(m,k),1) + m - 1
//               don't interchange rows
			l = m;
			ipvt[k] = k;

			/*
			 * zero pivot implies this column already triangularized
			 */
			if (abd[l + k * sfda] == 0.0) {
				*info = k;
				continue;
			}

			/*
			 * interchange if necessary
			 */
//			if (l != m) {
//				t = abd[l + k * sfda];
//				abd[l + k * sfda] = abd[m + k * sfda];
//				abd[m + k * sfda] = t;
//			}

			/*
			 * compute multipliers
			 */
			t = -1.0/abd[m + k * sfda];
			sscal_(&lm, &t, &abd[m + 1 + k * sfda], &one);

			/*
			 * row eliminination with column indexing
			 */
			// computing MIN
			// computing MAX
			ju = min(max(ju, (sud + ipvt[k])) , *n);
			mm = m;
			if (ju >= kp1) {
				for (j = kp1; j <= ju; j++) {
					l--;
					mm--;
					t = abd[l + j * sfda];
					if (l != mm) {
						abd[l + j * sfda] = abd[mm + j * sfda];
						abd[mm + j * sfda] = t;
					}
					saxpy_(&lm, &t, &abd[m + 1 + k * sfda], &one, 
							&abd[mm + 1 + j * sfda], &one);
				}
			}
		}
	}
	ipvt[*n] = *n;
	if (abd[m + *n * sfda] == 0.0) {
		*info = *n;
	}
	return 0;
}
