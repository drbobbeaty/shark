/*
 * saxpy.c
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

/* 	Public Data Types 		*/

/* 	Public Data Constants 	*/

/* ----------------------------------------------------------------------------
 *
 *   computer            - cdc/single
 *
 *   latest revision     - january 1, 1978
 *
 *   purpose             - compute a constant times a vector plus
 *                         a vector, all single precision
 *
 *   usage               - call saxpy (n,sa,sx,incx,sy,incy)
 *
 *   arguments    n      - length of vectors x and y. (input)
 *                sa     - real scalar. (input)
 *                sx     - real vector of length max(n*iabs(incx),1).
 *                         (input)
 *                incx   - displacement between elements of sx. (input)
 *                         x(i) is defined to be..
 *                         sx(1+(i-1)*incx) if incx .ge. 0 or
 *                         sx(1+(i-n)*incx) if incx .lt. 0.
 *                sy     - real vector of length max(n*iabs(incy),1).
 *                         (input/output)
 *                         saxpy replaces y(i) with sa*x(i)+y(i) 
 *                         for i=1,...,n.
 *                         x(i) and y(i) refer to specific elements
 *                         of sx and sy, respectively. see incx and 
 *                         incy argument descriptions.
 *                incy   - displacment between elements of sy. (input)
 *                         y(i) is defined to be..
 *                         sy(1+(i-1)*incy) if incy .ge. 0 or
 *                         sy(1+(i-n)*incy) if incy .lt. 0.
 *
 *   precision/hardware  - single/all
 *
 *   reqd. imsl routines - none required
 *
 *   notation            - information on special notation and
 *                         conventions is available in the manual
 *                         introduction or through imsl routine uhelp
 *
 *   copyright           - 1978 by imsl, inc. all rights reserved
 *
 *   warranty            - imsl warrants only that imsl testing has been
 *                         applied to this code. no other warranty,
 *                         expressed or implied, is applicable.
 *
 * ----------------------------------------------------------------------------
 */
int saxpy_(integer *n, doublereal *sa, doublereal *sx, integer *incx, doublereal *sy, integer *incy) {
	int			retval = 0;
	/* dereference the pointers for faster access later on */
	integer		count = *n;
	doublereal	a = *sa;
	integer		dix = *incx;
	integer		diy = *incy;

	/*
	 * First, see if we have anything to do
	 */
	if ((count > 0) && (a != 0.0)) {
		/*
		 * See if it's an easy sub-case of the general class
		 */
		if ((dix == 1) && (diy == 1)) {
			/*
			 * OK... we have the simple, common, case where we have
			 * equal, unit steps and we need to do this as fast as
			 * possible. It turns out that the FORTRAN code was right -
			 * we need to do this in batches. So, finish off the odd
			 * ones, and then batch the rest.
			 */
			integer		m = 0;
			integer		mm = count - (count/4)*4;

			/*
			 * OK... let's do the odd balls first...
			 */
			for (; m < mm; m++) {
				sy[m] += a * sx[m];
			}
			/*
			 * ...now let's finish off the rest in blocks
			 */
			if (count >= 4) {
				for (; m < count; m += 4) {
					sy[m] += a * sx[m];
					sy[m+1] += a * sx[m+1];
					sy[m+2] += a * sx[m+2];
					sy[m+3] += a * sx[m+3];
				}
			}
		} else {
			integer		m = 0;
			integer		ix = 0;
			integer 	iy = 0;

			/*
			 * This is the general routine
			 */
			for (; m < count; m++) {
				sy[iy] += a * sx[ix];
				ix += dix;
				iy += diy;
			}
		}
	}
	return retval;
}
