#ifdef Cyber
      subroutine linpk( abd, lda, sd, error, info )
#else
      subroutine linpk( abd, lda, sd, ipvt, error, info )
#endif
c
c     This subroutine solves the Ju=F system by using the
c     banded solution routines in LINPACK, or MAGEV on the Cyber.
c
c     These routines are the double precision LINPACK routines
c     for solving matricies, or bdgel() in MAGEV.
c
c     Get the precision of the machine
c
#     ifdef Cyber
#       define PRECISION real
#     else
#       define PRECISION double precision
#     endif
c
c     first the parameters for the sizes of the arrays
c
      integer Nx, Ny
      parameter( Nx=XMAX, Ny=YMAX )
c
c     Now for some parameters that are constants at the machine's
c     precision
c
      PRECISION tenm10
      parameter( tenm10=1.0e-10 )

#     ifdef NP1
#       ifndef NOVECTOR
c
c         The NP1 is supposed to do direct vector statements.
c
#         define NP1vec
#       endif
#     endif
c
c     now for the arguments to this routine
c
      integer lda, sd, info
#     ifdef Cyber
      PRECISION abd( lda, sd )
#     else
      integer ipvt( sd )
      PRECISION abd( lda, sd )
#     endif
      logical error
c
c     next the common blocks for interface to the rest of the routines
c
      common /mesh/ hx(Nx), hy(Ny), xnodes, ynodes, xfrac, auto, dVa
      integer xnodes, ynodes
      PRECISION hx, hy, xfrac, dVa
      character*8 auto
      common /JF/ J( (2*Nx*Ny), 10 ), delu( 2*Nx*Ny ), F( 2*Nx*Ny )
      PRECISION J, delu, F
c
c     now for the local variables used by the LINPACK routines
c
      integer n, ml, mu, m, i, p, i1, i2, k, l
#     ifdef Cyber
      PRECISION det
#     endif
c
c     compute the constants I need
c
      n = 2 * ynodes * xnodes
      ml = 2*min(xnodes,ynodes) + 1
      mu = 2*min(xnodes,ynodes) + 1
#     ifdef Cyber
      m = mu + 1
#     else
      m = ml + mu + 1
#     endif
c
c     we need to move J to abd for use by sgbfa and sgbsl, or bdgel
c
      do 200 k=1,n
        i1 = max( 1, (k-mu) )
        i2 = min( n, (k+ml) )
        do 100 i=i1,i2
          p = i - k + m
c               calculate the l-value
          l = int( (i-1)/2 ) + 1
c               sort out the value in J()
          if ( k .eq. ( 2*(l-min(xnodes,ynodes)) - 1 ) ) then
            abd( p,k ) = J( i,1 )
          elseif ( k .eq. ( 2*(l-min(xnodes,ynodes)) ) ) then
            abd( p,k ) = J( i,2 )
          elseif ( k .eq. ( 2*l - 3 ) ) then
            abd( p,k ) = J( i,3 )
          elseif ( k .eq. ( 2*l - 2 ) ) then
            abd( p,k ) = J( i,4 )
          elseif ( k .eq. ( 2*l - 1 ) ) then
            abd( p,k ) = J( i,5 )
          elseif ( k .eq. ( 2*l ) ) then
            abd( p,k ) = J( i,6 )
          elseif ( k .eq. ( 2*l + 1 ) ) then
            abd( p,k ) = J( i,7 )
          elseif ( k .eq. ( 2*l + 2 ) ) then
            abd( p,k ) = J( i,8 )
          elseif ( k .eq. ( 2*(l+min(xnodes,ynodes)) - 1 ) ) then
            abd( p,k ) = J( i,9 )
          elseif ( k .eq. ( 2*(l+min(xnodes,ynodes)) ) ) then
            abd( p,k ) = J( i,10 )
          else
c               all else is zero
            abd( p,k ) = 0.0
          endif
  100   continue
  200 continue
c
c     Now solve the RHS vector to delu(i) for solution
c
      do 300 i=1,n
  300   delu(i) = F(i)

#     ifdef Cyber
c
c     Now factor and solve the system of equations with MAGEV
c
      call bdgel( abd, delu, n, mu, ml, 1, det, tenm10, .false. )
      if ( det .eq. 0.0 ) then
        error = .true.
        return
      endif
#     else
#     ifdef LAPACK
c
c     Now factor and solve the system of equations with BLAS
c     and dgbsv().
c
      call dgbsv( n, ml, mu, 1, abd, lda, ipvt, delu, n, info )
      if ( info .ne. 0 ) then
        error = .true.
        return
      endif
#     else
c
c     now factor the system of equations
c
      call sgbfa( abd, lda, n, ml, mu, ipvt, info )
      if ( info .ne. 0 ) then
        error = .true.
        return
      endif
c
c     now solve the system
c
      call sgbsl( abd, lda, n, ml, mu, ipvt, delu, 0 )
#     endif
#     endif
c
c     ...and now we can return to the caller
c
      error = .false.
      return
      end



#     ifdef Cyber
      subroutine Lpk1D( abd, lda, sd, error info )
#     else
      subroutine Lpk1D( abd, lda, sd, ipvt, error, info )
#     endif
c
c     This subroutine solves the 1D Ju=F system by using the
c     banded solution routines in LINPACK, or MAGEV on the Cyber.
c
c     These routines are the double precision LINPACK routines
c     for solving matricies, or bdgel() in MAGEV.
c
c     first the parameters for the sizes of the arrays
c
      integer Nx, Ny
      parameter( Nx=XMAX, Ny=YMAX )

#     ifdef NP1
#       ifndef NOVECTOR
c
c         The NP1 is supposed to do direct vector statements.
c
#         define NP1vec
#       endif
#     endif
c
c     now for the arguments to this routine
c
      integer lda, sd, info
#     ifdef Cyber
      PRECISION abd( lda, sd )
#     else
      integer ipvt( sd )
      PRECISION abd( lda, sd )
#     endif
      logical error
c
c     next the common blocks for interface to the rest of the routines
c
      common /mesh/ hx(Nx), hy(Ny), xnodes, ynodes, xfrac, auto, dVa
      integer xnodes, ynodes
      PRECISION hx, hy, xfrac, dVa
      character*8 auto
      common /JF/ J( (2*Nx*Ny), 10 ), delu( 2*Nx*Ny ), F( 2*Nx*Ny )
      PRECISION J, delu, F
c
c     now for the local variables used by the LINPACK routines
c
      integer n, ml, mu, m, i, p, i1, i2, k, l
#     ifdef Cyber
      PRECISION det
#     endif
c
c     compute the constants I need for the 1D solve
c
      n = 2 * ynodes
      ml = 3
      mu = 3
#     ifdef Cyber
      m = mu + 1
#     else
      m = ml + mu + 1
#     endif
c
c     we need to move J to abd for use by sgbfa and sgbsl, or bdgel
c
      do 200 k=1,n
        i1 = max( 1, (k-mu) )
        i2 = min( n, (k+ml) )
        do 100 i=i1,i2
          p = i - k + m
c               calculate the l-value
          l = int( (i-1)/2 ) + 1
c               sort out the value in J()
          if ( k .eq. ( 2*l - 3 ) ) then
            abd( p,k ) = J( i,1 )
          elseif ( k .eq. ( 2*l - 2 ) ) then
            abd( p,k ) = J( i,2 )
          elseif ( k .eq. ( 2*l - 1 ) ) then
            abd( p,k ) = J( i,3 )
          elseif ( k .eq. ( 2*l ) ) then
            abd( p,k ) = J( i,4 )
          elseif ( k .eq. ( 2*l + 1 ) ) then
            abd( p,k ) = J( i,5 )
          elseif ( k .eq. ( 2*l + 2 ) ) then
            abd( p,k ) = J( i,6 )
          else
c               all else is zero
            abd( p,k ) = 0.0
          endif
  100   continue
  200 continue
c
c     Now solve the RHS vector to delu(i) for solution
c
      do 300 i=1,n
  300   delu(i) = F(i)

#     ifdef Cyber
c
c     Now factor and solve the system of equations with MAGEV
c
      call bdgel( abd, delu, n, mu, ml, 1, det, tenm10, .false. )
      if ( det .eq. 0.0 ) then
        error = .true.
        return
      endif
#     else
#     ifdef LAPACK
c
c     Now factor and solve the system of equations with BLAS
c     and dgbsv().
c
      call dgbsv( n, ml, mu, 1, abd, lda, ipvt, delu, n, info )
      if ( info .ne. 0 ) then
        error = .true.
        return
      endif
#     else
c
c     now factor the system of equations
c
      call sgbfa( abd, lda, n, ml, mu, ipvt, info )
      if ( info .ne. 0 ) then
        error = .true.
        return
      endif
c
c     now solve the system
c
      call sgbsl( abd, lda, n, ml, mu, ipvt, delu, 0 )
#     endif
#     endif
c
c     ...and now we can return to the caller
c
      error = .false.
      return
      end
