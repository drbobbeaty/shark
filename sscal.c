/*
 * sscal.c
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

/* 	Forward Declarations 	*/

/* 	Public Data Types 		*/

/* 	Public Data Constants 	*/


/* ----------------------------------------------------------------------------
 *
 *   computer            - cdc/single
 *
 *   latest revision     - january 1, 1978
 *
 *   purpose             - compute a single precision constant
 *                         times a single precision vector
 *
 *   usage               - call sscal (n,sa,sx,incx)
 *
 *   arguments    n      - length of vector x. (input)
 *                sa     - real scalar. (input)
 *                sx     - real vector of length max(n*iabs(incx),1).
 *                         (input)
 *                         sscal replaces x(i) with sa*x(i) for
 *                         i=1,...,n.
 *                         x(i) refers to a specific element of sx.
 *                         see incx argument description.
 *                incx   - displacement between elements of sx. (input)
 *                         x(i) is defined to be sx(1+(i-1)*incx).
 *                         incx must be greater than zero.
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
int sscal_(long int *n, double *sa, double *sx, long int *incx) {
	int		retval = 0;

	/*
	 * First, see if we have anything to do...
	 */
    if (*n > 0) {
    	long int		i = 0;
		long int		di = *incx;
		long int		ns = (*n) * di;
		double		 	a = *sa;

		for (; i < ns; i += di) {
			sx[i] *= a;
   	 	}
	}
    return retval;
}
