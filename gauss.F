      subroutine gauss( error, iter )
c
c     This routine solves a sparse linear system of equations in the
c     Shark format (see notes of R.E.Beaty for this format) by the
c     iterative Gauss-Siedel method
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
      PRECISION omega
      parameter( omega=1.55 )
c
c     now for the arguments to this routine
c
      logical error
      integer iter
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
      integer i, k, l, N, pos, look(10)
      PRECISION sume, sumo, maxchg
c
c     Preset a few local variables so they don't give the compiler trouble
c
      i = 0
c
c     precompute the number of equations I will be solving
c
      N = 2 * xnodes * ynodes
c
c     precompute the look-up table for the compressed array
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
c     prepare [J] for the solution of Gauss-Siedel, and
c     guess the solution to be near F()
c
      do 100 l=2,N,2
c				preset the answer to something close
        delu( l-1 ) = F( l-1 )
        delu( l ) = F( l )
c				zero out the diaginal of [J]
        J( (l-1), 5 ) = 0.0
        J( l, 6 ) = 0.0
c				negate the other elements of [J]
        do 100 i=1,10
          J( (l-1), i ) = -J( (l-1), i )
  100     J( l, i ) = -J( l, i )
c
c     iterate on delu(i) until converged, by the Gauss-Seidel method
c
      do 400 iter=1,(10*maxit)
c
c       reset the maximum change for this downward pass on the system
c
        maxchg = 0.0
c
c       update this matrix for one pass by first getting the correction
c
        do 300 l=2,N,2
          sumo = 0.0
          sume = 0.0
          do 200 k=1,10
            pos = l + look(k)
            if ( (pos .ge. 1) .and. (pos .le. N) ) then
              sumo = sumo + J( (l-1), k ) * delu( pos )
              sume = sume + J( l, k ) * delu( pos )
            endif
  200     continue
c
c         Now get the maximum correction
c
          sumo = F( l-1 ) + sumo
          sume = F( l ) + sume
          maxchg = max( maxchg, abs(sumo - delu(l-1)) )
          maxchg = max( maxchg, abs(sume - delu(l)) )
c
c         now correct these two elements
c
          if ( abs(sumo - delu(l-1)) .gt. (vmax/10.0) ) then
            delu( l-1 ) = delu( l-1 ) + omega*(sumo - delu(l-1))
          else
            delu( l-1 ) = sumo
          endif
          if ( abs(sume - delu(i)) .gt. (vmax/10.0) ) then
            delu( l ) = delu( l ) + omega*(sume - delu(l))
          else
            delu( l ) = sume
          endif
  300   continue
c
c       check for convergence on this downward pass
c
        if ( maxchg .le. (vmax * 1.0e-2) ) then
          error = .false.
          return
        endif
c
c     done with one loop of the iteration on the solution
c
  400 continue
c
c     failed to converge in 10*maxit iterations
c
      error = .true.
      return
      end
