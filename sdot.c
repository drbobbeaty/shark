/*
 * sdot.c
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
 *   purpose             - compute single precision dot product
 *
 *   usage               - function sdot (n,sx,incx,sy,incy)
 *
 *   arguments    sdot   - sum from i=1 to n of x(i)*y(i). (output)
 *                         x(i) and y(i) refer to specific elements
 *                         of sx and sy, respectively. see incx and
 *                         incy argument descriptions.
 *                n      - length of vectors x and y. (input)
 *                sx     - real vector of length max(n*iabs(incx),1).
 *                         (input)
 *                incx   - displacement between elements of sx. (input)
 *                         x(i) is defined to be..
 *                         sx(1+(i-1)*incx) if incx .ge. 0 or
 *                         sx(1+(i-n)*incx) if incx .lt. 0.
 *                sy     - real vector of length max(n*iabs(incy),1).
 *                         (input)
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
doublereal sdot_(integer *n, doublereal *sx, integer *incx, doublereal *sy, integer *incy) {
    doublereal 	retval = 0;
	/* dereference some arguments for speed of access later */
	integer		count = *n;
	integer		dix = *incx;
	integer 	diy = *incy;

	/*
	 * First, see if we have anything to do
	 */
	if (count > 0) {
		integer	m = 0;
		integer	ix = 0;
		integer iy = 0;

		for (m = 0; m < count; m++) {
			retval += sx[ix] * sy[iy];
			ix += dix;
			iy += diy;
		}
	}
	return retval;
}
