      integer function isamax (n,sx,incx)
c
c----------------------------------------------------------------------------
c
c   computer            - cdc/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - find the smallest index of the maximum magnitude
c                         of a single precision vector
c
c   usage               - function isamax (n,sa,incx)
c
c   arguments    isamax - the smallest index i such that abs(x(i))
c                         is the maximum of abs(x(j)) for j=1 to n.
c                         (output)
c                         x(i) refers to a specific element of sx.
c                         see incx argument description.
c                n      - length of vector x. (input)
c                sx     - real vector of length n*incx. (input)
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
c
c               specifications for arguments
      integer          n,incx 
      double precision sx(1)
c               specifications for local variables
      integer          i, ii, ns
      double precision smax, xmag
c               first executable statement
      isamax = 0
      if (n .le. 0) return
      isamax = 1
      if (n .le. 1) return
      if (incx .eq. 1) go to 15
c               code for increments not equal to 1
      smax = abs(sx(1))
      ns = n*incx
      ii = 1
      do 10 i=1,ns,incx
        xmag = abs(sx(i))
        if (xmag .le. smax) go to 5
        isamax = ii
        smax = xmag
    5   ii = ii + 1
   10 continue
      return
c               code for increments equal to 1
   15 smax = abs(sx(1))
      do 20 i=2,n
        xmag = abs(sx(i))
        if (xmag .le. smax) go to 20
        isamax = i
        smax = xmag
   20 continue
      return
      end
