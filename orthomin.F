      subroutine ortho( error, iter )
c
c     This routine solves a sparse linear system of equations in the
c     Shark format (see notes of R.E.Beaty for this format) by the 
c     iterative Orthomin(k) method.
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
      integer keepmx, Nx, Ny
      parameter( keepmx=6, Nx=XMAX, Ny=YMAX )
c
c     now for the arguments to this routine
c
      integer iter
      logical error
c
c     next the common blocks for interface to the rest of the routines
c
      common /mesh/ hx(Nx), hy(Ny), xnodes, ynodes, xfrac, auto, dVa
      integer xnodes, ynodes
      PRECISION hx, hy, xfrac, dVa
      character*8 auto
      common /params/ kT, nio, Ldi, kox, chmod, doptyp, ksemi,
     +                dt, tstart, tstop, stpmin, tout, delt,
     +                Vs, Vd, vg1, eavg, soltyp, symtyp,
     +                icchan, chtyp, instyp, vmax, debug,
     +                guess, delV, delE, maxit, keep
      integer maxit, keep
      PRECISION kT, nio, Ldi, kox, ksemi, dt, tstart, tstop,
     +     stpmin, tout, Vs, Vd, vg1, eavg, vmax, delt,
     +     delV, delE
      character*8 chmod, doptyp, icchan, chtyp, instyp, debug, guess,
     +     soltyp, symtyp
      common /JF/ J( (2*Nx*Ny), 10 ), delu( 2*Nx*Ny ), F( 2*Nx*Ny )
      PRECISION J, delu, F
c
c     now for the local variables
c
      integer i, k, l, N, pos, cycpos, look(10), t
      PRECISION p( (2*Nx*Ny), keepmx ), q( (2*Nx*Ny), keepmx ),
     +          temp( 2*Nx*Ny ), r( 2*Nx*Ny ), qi( keepmx ), b( keepmx )
      PRECISION LU( (2*Nx*Ny), 10 )
      PRECISION sume, sumo, ai, maxres
c
c     limit keep to keepmx
c
      if ( keep .gt. keepmx ) keep = keepmx
c
c     Make sure we aren't supposed to be in minres()
c
      if ( keep .eq. 0 ) then
        call minres( error, iter )
        return
      endif
c
c     Precompute the number of equations I will be solving
c
      N = 2 * xnodes * ynodes
c
c     Precompute the look-up table for the compressed array
c
      look(1) = -2*min(xnodes,ynodes) - 1
      look(2) = -2*min(xnodes,ynodes)
      look(3) = -3
      look(4) = -2
      look(5) = -1
      look(6) = 0
      look(7) = 1
      look(8) = 2
      look(9) = 2*min(xnodes,ynodes) - 1
      look(10) = 2*min(xnodes,ynodes)
c
c     we need the LU factorization of the system matrix for the
c     preconditioning of the system
c
      do 500 l=2,N,2
        do 500 k=1,10
          pos = l + look(k)
          if ( (pos .ge. 1) .and. (pos .le. N) ) then
c				reset the sum counts
            sumo = 0.0
            sume = 0.0
c				odd-lined elements
            if ( k .le. 6 ) then
              do 200 t=1,min((k-1),4)
                if ( (l+look(t)) .ge. 1 ) then
                  do 100 i=1,10
                    if ( look(i) .eq. 
     +                   (pos - int((l+look(t)+1)/2)*2) ) then
                      sumo = sumo + LU(l-1,t) * LU( (l+look(t)),i )
                      goto 200
                    endif
  100             continue
                endif
  200         continue
            endif
            LU( l-1, k ) = J( l-1, k ) - sumo
c				divide for the elements of U
            if ( k .ge. 6 ) LU(l-1,k) = LU(l-1,k) / LU(l-1,5)

c				even-lined, lower diagonal elements
            if ( k .le. 6 ) then
              do 400 t=1,(k-1)
                if ( (l+look(t)) .ge. 1 ) then
                  do 300 i=1,10
                    if ( look(i) .eq. 
     +                   (pos - int((l+look(t)+1)/2)*2) ) then
                      sume = sume + LU(l,t) * LU( (l+look(t)),i )
                      goto 400
                    endif
  300             continue
                endif
  400         continue
            else
c				even-lined, upper diagonal elements
              sume = LU( (l-1), k ) * LU( l,5 )
            endif
            LU( l, k ) = J( l, k ) - sume
c				divide for the elements of U
            if ( k .ge. 7 ) LU(l,k) = LU(l,k) / LU(l,6)
          endif
  500 continue
c
c     guess the solution to be near convergence
c
      do 600 l=1,N
  600   delu(l) = 0.0
c
c     compute the initial residual without the preconditioning treatment
c
      do 800 l=2,N,2
        sumo = 0.0
        sume = 0.0
        do 700 k=1,10
          pos = l + look(k)
          if ( (pos .ge. 1) .and. (pos .le. N) ) then
            sumo = sumo + J( (l-1), k )*delu( pos )
            sume = sume + J( l, k )*delu( pos )
          endif
  700   continue
        r( l-1 ) = F( l-1 ) - sumo
        r( l ) = F( l ) - sume
  800 continue
c
c     run the initial residual through the preconditioning
c
c				first the forward substitution
      do 1100 l=2,N,2
c				...odd lines first
        do 900 k=1,4
          pos = l + look(k)
          if ( pos .ge. 1 ) then
            r(l-1) = r(l-1) - LU( l-1, k ) * r(pos)
          endif
  900   continue
        r(l-1) = r(l-1) / LU( l-1, 5 )
c				...even lines next
        do 1000 k=1,5
          pos = l + look(k)
          if ( pos .ge. 1 ) then
            r(l) = r(l) - LU( l, k ) * r(pos)
          endif
 1000   continue
        r(l) = r(l) / LU( l, 6 )
 1100 continue
c				now the reverse substitution
      do 1400 l=N,2,-2
c				...even lines first because they are lower
        do 1200 k=7,10
          pos = l + look(k)
          if ( pos .le. N ) then
            r(l) = r(l) - LU( l, k ) * r(pos)
          endif
 1200   continue
c				...odd lines next
        do 1300 k=6,10
          pos = l + look(k)
          if ( pos .le. N ) then
            r(l-1) = r(l-1) - LU( l-1, k ) * r(pos)
          endif
 1300   continue
 1400 continue
c
c     set p(1) = initial residual
c
      do 1500 i=1,N
 1500   p(i,1) = r(i)
c
c     iterate on p(i) and delu() until converged
c
      do 4000 iter=1,maxit
c				get the cyclic position in keep
        cycpos = mod( (iter-1), keep ) + 1
c
c       calculate J*p(i) and store it in q(i) without the preconditioning
c
        do 1700 l=2,N,2
          sumo = 0.0
          sume = 0.0
          do 1600 k=1,10
            pos = l + look(k)
            if ( (pos .ge. 1) .and. (pos .le. N) ) then
              sumo = sumo + J( (l-1), k )*p( pos, cycpos )
              sume = sume + J( l, k )*p( pos, cycpos )
            endif
 1600     continue
          q( (l-1), cycpos ) = sumo
          q( l, cycpos ) = sume
 1700   continue
c
c       run q(i) through the preconditioning
c
c				first the forward substitution
        do 2000 l=2,N,2
c				...odd lines first
          do 1800 k=1,4
            pos = l + look(k)
            if ( pos .ge. 1 ) then
              q( l-1, cycpos ) = q( l-1, cycpos ) - LU( l-1, k )
     +                           * q( pos, cycpos )
            endif
 1800     continue
          q( l-1, cycpos ) = q( l-1, cycpos ) / LU( l-1, 5 )
c				...even lines first
          do 1900 k=1,5
            pos = l + look(k)
            if ( pos .ge. 1 ) then
              q( l, cycpos ) = q( l, cycpos ) - LU( l, k )
     +                         * q( pos, cycpos )
            endif
 1900     continue
          q( l, cycpos ) = q( l, cycpos ) / LU( l, 6 )
 2000   continue
c				now the reverse substitution
        do 2300 l=N,2,-2
c				...even lines first because they are lower
          do 2100 k=7,10
            pos = l + look(k)
            if ( pos .le. N ) then
              q( l, cycpos ) = q( l, cycpos ) - LU( l, k )
     +                         * q( pos, cycpos )
            endif
 2100     continue
c				...odd lines next
          do 2200 k=6,10
            pos = l + look(k)
            if ( pos .le. N ) then
              q( l-1, cycpos ) = q( l-1, cycpos ) - LU( l-1, k )
     +                         * q( pos, cycpos )
            endif
 2200     continue
 2300   continue
c
c       save the value of q(i)**2 after it is through the
c       preconditioning
c
        qi( cycpos ) = 0.0
        do 2400 l=1,N
 2400     qi( cycpos ) = qi( cycpos ) + q(l,cycpos)*q(l,cycpos)
c
c       calculate ai
c
        sumo = 0.0
        do 2500 i=1,N
 2500     sumo = sumo + r(i)*q( i, cycpos )
        ai = sumo/qi( cycpos )
c
c       calculate new solution vector
c
        do 2600 i=1,N
 2600     delu(i) = delu(i) + ai*p( i, cycpos )
c
c       calculate new residual
c
        maxres = 0.0
        do 2700 i=1,N
          r(i) = r(i) - ai*q( i, cycpos )
          maxres = max( maxres, abs(r(i)) )
 2700   continue
c
c       check for convergence
c
        if ( maxres .le. (vmax/100.0) ) then
          error = .false.
          return
        endif
c
c       calculate temp for the calculation of the b's
c
        do 2900 l=2,N,2
          sumo = 0.0
          sume = 0.0
          do 2800 k=1,10
            pos = l + look(k)
            if ( (pos .ge. 1) .and. (pos .le. N) ) then
              sumo = sumo + J( (l-1), k )*r( pos )
              sume = sume + J( l, k )*r( pos )
            endif
 2800     continue
          temp( l-1 ) = sumo
          temp( l ) = sume
 2900   continue
c
c       pass temp(i) through the preconditioning as well
c
c				first the forward substitution
        do 3200 l=2,N,2
c				...odd lines first
          do 3000 k=1,4
            pos = l + look(k)
            if ( pos .ge. 1 ) then
              temp(l-1) = temp(l-1) - LU( l-1, k ) * temp(pos)
            endif
 3000     continue
          temp(l-1) = temp(l-1) / LU( l-1, 5 )
c				...even lines first
          do 3100 k=1,5
            pos = l + look(k)
            if ( pos .ge. 1 ) then
              temp(l) = temp(l) - LU( l, k ) * temp(pos)
            endif
 3100     continue
          temp(l) = temp(l) / LU( l, 6 )
 3200   continue
c				now the reverse substitution
        do 3500 l=N,2,-2
c				...even lines first because they are lower
          do 3300 k=7,10
            pos = l + look(k)
            if ( pos .le. N ) then
              temp(l) = temp(l) - LU( l, k ) * temp(pos)
            endif
 3300     continue
c				...odd lines next
          do 3400 k=6,10
            pos = l + look(k)
            if ( pos .le. N ) then
              temp(l-1) = temp(l-1) - LU( l-1, k ) * temp(pos)
            endif
 3400     continue
 3500   continue
c
c       calculate the b's
c       make sure the loop is performed at least once
c
        do 3700 i=(max(1,(iter-keep+1))),iter
          sumo = 0.0
          cycpos = mod( (i-1), keep ) + 1
          do 3600 k=1,N
 3600       sumo = sumo + temp(k)*q( k, cycpos )
 3700     b( cycpos ) = -sumo/qi( cycpos )
c
c       calculate p(i+1)
c
        pos = mod( iter, keep ) + 1
        do 3900 i=1,N
          sumo = 0.0
          do 3800 k=(max(1,(iter-keep+1))),iter
            cycpos = mod( (k-1), keep ) + 1
 3800       sumo = sumo + b( cycpos )*p( i, cycpos )
 3900     p( i, pos ) = r(i) + sumo
c
c     done with one loop of the iteration on the solution
c
 4000 continue
c
c     failed to converge in maxit iterations
c
      error = .true.
      return
      end
