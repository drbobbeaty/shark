      subroutine sscal (n,sa,sx,incx)
c
c----------------------------------------------------------------------------
c
c   computer            - cdc/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - compute a single precision constant
c                         times a single precision vector
c
c   usage               - call sscal (n,sa,sx,incx)
c
c   arguments    n      - length of vector x. (input)
c                sa     - real scalar. (input)
c                sx     - real vector of length max(n*iabs(incx),1).
c                         (input)
c                         sscal replaces x(i) with sa*x(i) for
c                         i=1,...,n.
c                         x(i) refers to a specific element of sx.
c                         see incx argument description.
c                incx   - displacement between elements of sx. (input)
c                         x(i) is defined to be sx(1+(i-1)*incx).
c                         incx must be greater than zero.
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
c               specifications for arguments
      integer     incx, n
      double precision   sa
      double precision   sx(1)
c               specifications for local variables
      integer     i, m, mp1, ns
c               first executable statement
      if (n .le. 0) return
      if (incx .eq. 1) go to 10
c               code for increments not eaual to 1
      ns = n*incx
      do 5 i=1,ns,incx
        sx(i) = sa*sx(i)
    5 continue
      return
c               code for increments equal to 1.
c               clean-up loop so remaining vector
c               length is a multiple of 5.
   10 m = n - (n/5)*5
      if (m .eq. 0) go to 20
      do 15 i=1,m
        sx(i) = sa*sx(i)
   15 continue
      if (n .lt. 5) return
   20 mp1 = m + 1
      do 25 i=mp1,n,5
        sx(i) = sa*sx(i)
        sx(i+1) = sa*sx(i+1)
        sx(i+2) = sa*sx(i+2)
        sx(i+3) = sa*sx(i+3)
        sx(i+4) = sa*sx(i+4)
   25 continue
      return
      end
