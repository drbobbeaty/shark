      subroutine saxpy (n,sa,sx,incx,sy,incy)
c
c----------------------------------------------------------------------------
c
c   computer            - cdc/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - compute a constant times a vector plus
c                         a vector, all single precision
c
c   usage               - call saxpy (n,sa,sx,incx,sy,incy)
c
c   arguments    n      - length of vectors x and y. (input)
c                sa     - real scalar. (input)
c                sx     - real vector of length max(n*iabs(incx),1).
c                         (input)
c                incx   - displacement between elements of sx. (input)
c                         x(i) is defined to be..
c                         sx(1+(i-1)*incx) if incx .ge. 0 or
c                         sx(1+(i-n)*incx) if incx .lt. 0.
c                sy     - real vector of length max(n*iabs(incy),1).
c                         (input/output)
c                         saxpy replaces y(i) with sa*x(i)+y(i)
c                         for i=1,...,n.
c                         x(i) and y(i) refer to specific elements
c                         of sx and sy, respectively. see incx and
c                         incy argument descriptions.
c                incy   - displacment between elements of sy. (input)
c                         y(i) is defined to be..
c                         sy(1+(i-1)*incy) if incy .ge. 0 or
c                         sy(1+(i-n)*incy) if incy .lt. 0.
c
c   precision/hardware  - single/all
c
c   reqd. imsl routines - none required
c
c   notation            - information on special notation and
c                         conventions is available in the manual
c                         introduction or through imsl routine uhelp
c
c   copyright           - 1978 by imsl, inc. all rights reserved
c
c   warranty            - imsl warrants only that imsl testing has been
c                         applied to this code. no other warranty,
c                         expressed or implied, is applicable.
c
c----------------------------------------------------------------------------
c
c
c               specifications for arguments
      integer    n, incx, incy
      double precision  sa
      double precision  sx(1), sy(1)
c               specifications for local variables
      integer    i, ix, iy, m, mp1, ns
c               first executable statement
      if (n .le. 0 .or. sa .eq. 0.0e0) return
      if (incx .eq. incy) if (incx-1) 5, 15, 35
    5 continue
c               code for non-equal or non-positive
c               increments
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx+1
      if (incy .lt. 0) iy = (-n+1)*incy+1
      do 10 i=1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c               code for both increments equal to 1
c               clean-up loop so remaining vector
c               length is a multiple of 4.
   15 m = n - (n/4)*4
      if (m .eq. 0) go to 25
      do 20 i=1,m
        sy(i) = sy(i) + sa*sx(i)
   20 continue
      if (n .lt. 4) return
   25 mp1 = m + 1
      do 30 i=mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i+1) = sy(i+1) + sa*sx(i+1)
        sy(i+2) = sy(i+2) + sa*sx(i+2)
        sy(i+3) = sy(i+3) + sa*sx(i+3)
   30 continue
      return
c               code for equal, positive, non-unit
c               increments.
   35 continue
      ns = n*incx
      do 40 i=1,ns,incx
        sy(i) = sa*sx(i) + sy(i)
   40 continue
      return
      end
