      subroutine DC( wVs, wVd, wVg1, wEavg, error )
c
c     This program solves for the applied bias by starting with
c     a small value of dt and increasing it until we have reached
c     a steady state situation. It calls newton() each time increment
c     for the solution of the equations.
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
      PRECISION m1psec, zero
      parameter( m1psec=-1.0e-12, zero=0.0 )
c
c     now for the arguments to this routine
c
      logical error
      PRECISION wVs, wVd, wVg1, wEavg
c
c     next, the common blocks for the necessary parts of the simulation
c
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
      common /struct/ lsorce, tsorce, dopsrc, ldrain, tdrain, 
     +                dopdrn, tox, tepi, length, tsub, dopsub,
     +                line, space, npulse, defoc, gsol, gdol,
     +                dopch, srccor, drncor, uins, Qss(Nx)
      PRECISION lsorce, tsorce, dopsrc, ldrain, tdrain, dopdrn,
     +     tox, tepi, length, tsub, dopsub, line, space, npulse,
     +     defoc, gsol, gdol, dopch, uins, Qss
      character*16 srccor, drncor
c
c     now for some local variables
c
      integer iter
      PRECISION ttotal, dtsave, stpcon
c
c     save the old value of dt because we will be changing it
c
      dtsave = dt
c
c     set the value of the step constant for the moving of dt
c
      stpcon = 2.0
c
c     start the big loop
c
      ttotal = 0.0
      dt = delt
  100 if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
        write(6,'('' '',''Solving for t='',g11.4,'' sec. with dt='',
     +          g11.4,'' sec.'')') ttotal, dt
      endif
#     ifndef Cyber
        call flush(6)
#     endif
  200 call newton( error, iter )
c
c     check to see if we had coupling problems
c
      if ( iter .eq. -1 ) then
        call refine( .false., .true., .true., error )
        if ( error ) then
          write(6,'('' ''.''**** Error! ****''/
     +              '' '',''     Shark could not resolve the coupling'',
     +     '' problems. This problem is not solvable''/
     +              '' '',''     with the current executable of Shark.'',
     +     '' The final output is being written and''/
     +              '' '',''     the solution halted.'')')
          write(9,'('' '',''#eor'')')
#         ifndef Cyber
            call flush(6)
            call dump( m1psec, .true., zero )
            close( 9 )
#         endif
          call depmap( 1,1 )
          stop
        endif
c
c       so try the bias point again...
c
        call mkdev( .false. )
        call doguess( wVs, wVd, wVg1, wEavg, .false., .false., error )
        if ( error ) return
        write(6,'('' '',''...Re-Trying the bias point...'')')
        goto 200
      endif
c
c     see if we have an error, and if so stop
c
      if ( error ) then
        dt = dtsave
        return
      endif
c
c     see if we are in steady state by looking at iter, and stop if
c     we are
c
      if ( iter .eq. 1 ) then
        write(6,'('' '',''We have reached steady state by virtue of'',
     +          '' newton() only taking one iteration.'')')
#       ifndef Cyber
          call flush(6)
#       endif
        dt = dtsave
        error = .false.
        return
      endif
c
c     move up ttotal and dt
c
      ttotal = ttotal + dt
      dt = dt * stpcon
c
c     see if we need to keep going...
c
c      if ( ttotal .le. 10*(length/1.0e7) ) goto 100
      if ( dt .le. 10*(length/1.0e7) ) goto 100
c
c     all done... so return to the caller
c
      dt = dtsave
      error = .false.
      return
      end
