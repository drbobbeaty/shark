      double precision function sdot (n,sx,incx,sy,incy)
c
c----------------------------------------------------------------------------
c
c   computer            - cdc/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - compute single precision dot product
c
c   usage               - function sdot (n,sx,incx,sy,incy)
c
c   arguments    sdot   - sum from i=1 to n of x(i)*y(i). (output)
c                         x(i) and y(i) refer to specific elements
c                         of sx and sy, respectively. see incx and
c                         incy argument descriptions.
c                n      - length of vectors x and y. (input)
c                sx     - real vector of length max(n*iabs(incx),1).
c                         (input)
c                incx   - displacement between elements of sx. (input)
c                         x(i) is defined to be..
c                         sx(1+(i-1)*incx) if incx .ge. 0 or
c                         sx(1+(i-n)*incx) if incx .lt. 0.
c                sy     - real vector of length max(n*iabs(incy),1).
c                         (input)
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
      integer     n, incx, incy
      double precision   sx(1), sy(1)
c               specifications for local variables
      integer     i, m, mp1, ns, ix, iy
c               first executable statement
      sdot = 0.0e0
      if (n .le. 0) return
      if (incx .eq. incy) if (incx-1) 5, 15, 35
    5 continue
c               code for unequal increments or
c               non-positive increments
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx+1
      if (incy .lt. 0) iy = (-n+1)*incy+1
      do 10 i=1,n
        sdot = sdot + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c               code for both increments equal to 1
c               clean-up loop so remaining vector
c               length is a multiple of 5.
   15 m = n - (n/5)*5
      if (m .eq. 0) go to 25
      do 20 i=1,m
        sdot = sdot + sx(i)*sy(i)
   20 continue
      if (n .lt. 5) return
   25 mp1 = m + 1
      do 30 i = mp1,n,5
        sdot = sdot + sx(i)*sy(i) + sx(i+1)*sy(i+1) + sx(i+2)*sy(i+2) +
     +                sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   30 continue
      return
c               code for positive equal increments
c               .ne. 1.
   35 continue
      ns = n*incx
      do 40 i=1,ns,incx
        sdot = sdot + sx(i)*sy(i)
   40 continue
      return
      end
