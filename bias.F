      subroutine bias( error )
c
c     This subroutine starts the simulation off by taking the initial
c     guess and moving it up to the desired operating point in such a
c     way so as not to cause  unstable results from newton() or from
c     the physics of the negative differential mobility of the GaAs.
c
c     This routine does two different algorithims depending on the
c     guess that was performed. If a charge neutral guess was done
c     then we start at equilibrium and move the bias voltages up slowly
c     enough so that we get no problems. At each point we call the
c     appropriate solver to keep the system up to date.
c
c     If we did a more advanced initial guess then the voltages are
c     ramped from that initial point to the operating point. This is
c     not presently supported
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
      PRECISION m1psec, zero, half, hund4
      parameter( m1psec=-1.0e-12, zero=0.0, half=0.5, hund4=400.0 )
c
c     now for the arguments to this routine
c
      logical error
c
c     next the common blocks for interface to the rest of the routines
c
      common /device/ n(Nx,Ny), dop(Nx,Ny), v(Nx,Ny), vgate(Nx),
     +                drncon(Ny,3), srccon(Ny,3), joxide
      integer joxide
      PRECISION n, dop, v, vgate, drncon, srccon
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
      common /struct/ lsorce, tsorce, dopsrc, ldrain, tdrain,
     +                dopdrn, tox, tepi, length, tsub, dopsub,
     +                line, space, npulse, defoc, gsol, gdol,
     +                dopch, srccor, drncor, uins, Qss(Nx)
      PRECISION lsorce, tsorce, dopsrc, ldrain, tdrain, dopdrn,
     +     tox, tepi, length, tsub, dopsub, line, space, npulse,
     +     defoc, gsol, gdol, dopch, uins, Qss
      character*16 srccor, drncor
      common /curves/ Do, dDo, Ad, Bd, mobil, Av, Bv, Cv, Dv, Fv
      PRECISION Do, dDo, Ad, Bd, mobil, Av, Bv, Cv, Dv, Fv
      common /Pshape/ Qm, ftw, fwhm, pulse
      PRECISION Qm, ftw, fwhm
      logical pulse
c
c     we need some local variables
c
      integer iter, i, k
      PRECISION wVs, wVd, wVg1, wEavg, vold(Nx,Ny), dVmax,
     +          dVold, dVb, dE, dVs, dVd, dVg1, dEavg,
     +          Vgmax, Vo, V1, Cox, Sxd, Dxd, Sphis, Dphis,
     +          Sxm, Dxm, Sphim, Dphim, nold(Nx,Ny)
c
c     Preset a few of the local variables so the compiler doesn't complain
c
      dVmax = 0.0
      dVb = 0.0
      dE = 0.0
      dVs = 0.0
      dVd = 0.0
      dVg1 = 0.0
      dEavg = 0.0
      Sphim = 0.0
      Dphim = 0.0
c
c     start with a header so that the output looks reasonable
c
      write(6,'('' '')')
      write(6,'('' '')')
      write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                            REB 1988'',
     +     ''    applying biases ===='')')
      write(6,'('' '')')
#     ifndef Cyber
        call flush(6)
#     endif
c
c     If we are doing a pulse, we do not need to do any
c     biasing, or steady-state solving, just the guess
c
      if ( pulse ) then
        call doguess( Vs, Vd, Vg1, eavg, .false., .true., error )
c
c       Save the initial guess for the curious debugger
c
        if ( debug .eq. 'yes' ) then
          call dump ( m1psec, .true., zero )
        endif
        return
      endif
c
c     See what kind of initial guess we have
c
      if ( guess .eq. 'MOS' ) then
c
c       A MOS guess means that the gate voltage can be set to its
c       final value, but the source and drain potentials need to
c       be ramped.
c
        if ( symtyp .eq. 'full' ) then
c
c         This is a MOS and junction initial guess. Non-trivial,
c         and also not presently working.
c
          if ( auto .ne. 'yes' ) then
            Vgmax = Vg1 + eavg*(length-(lsorce-gsol)-(ldrain-gdol))
          else
            Vgmax = Vg1 + eavg*(length+gdol+gsol)
          endif
          wVs = 10.0
          wVd = 10.0
        else
c
c         Get some important MOS parameters that I will
c         need to get the guessing right
c
          Cox = kox * 8.854e-14 / tox
          Vgmax = Vg1 + eavg*length + 1.6022e-19*Qss(xnodes)/Cox
          Vo = 1.6022e-19 * ksemi * 8.854e-14 * dopsub / (Cox*Cox)
          V1 = (1.6022e-19*tepi/Cox)*( dopsub + dopch )
          Sphis = (Vg1 + 1.6022e-19*Qss(1)/Cox) + Vo + V1 - sqrt( Vo*Vo
     +        + 2.0*Vo*( (Vg1 + 1.6022e-19*Qss(1)/Cox) + V1 )
     +        + V1*1.6022e-19*tepi*dopsub/Cox )
          Dphis = Vgmax + Vo + V1 - sqrt( Vo*Vo + 2.0*Vo*( Vgmax+V1 )
     +        + V1*1.6022e-19*tepi*dopsub/Cox )
          Sxd = sqrt( tepi*tepi*(1.0 + dopch/dopsub)
     +            + 2.0*ksemi*8.854e-14*Sphis/( 1.6022e-19 * dopsub) )
          Dxd = sqrt( tepi*tepi*(1.0 + dopch/dopsub)
     +            + 2.0*ksemi*8.854e-14*Dphis/( 1.6022e-19 * dopsub) )
          if ( chmod .eq. 'buried' ) then
            Sxm = tepi - (dopsub/dopch)*(Sxd-tepi)
            Dxm = tepi - (dopsub/dopch)*(Dxd-tepi)
c				see if we need just the extremum
            if ( Sxm .lt. 0.0 ) Sxm = 0.0
            if ( Dxm .lt. 0.0 ) Dxm = 0.0
            Sphim = (1.6022e-19*dopsub/(2.0*ksemi*8.854e-14))*
     +           ( Sxd*Sxd - tepi*tepi - 2.0*Sxm*(Sxd-tepi) ) -
     +           ( 1.6022e-19*dopch/(2.0*ksemi*8.854e-14))*
     +           ( tepi - Sxm )*( tepi - Sxm )
            Dphim = (1.6022e-19*dopsub/(2.0*ksemi*8.854e-14))*
     +           ( Dxd*Dxd - tepi*tepi - 2.0*Dxm*(Dxd-tepi) ) -
     +           ( 1.6022e-19*dopch/(2.0*ksemi*8.854e-14))*
     +           ( tepi - Dxm )*( tepi - Dxm )
          endif
c
c         Start these guys below the minimum of the channel, for each's
c         type of channel
c
          if ( chmod .eq. 'buried' ) then
            if ( soltyp .eq. 'linpk' ) then
c				see if v(E) is linear
              if ( (Av .ge. 1.0e38) .and. (Cv .ge. 1.0e38) ) then
                wVs = Vs
                wVd = Vd
              else
                wVs = 0.75 * Sphim
                wVd = 0.75 * Dphim
              endif
            else
              wVs = Sphim + 1.414/8.0
              wVd = Dphim + 1.414/8.0
            endif
          else
c
c           Then it must be a MODFET channel, and I know I can go to
c           the final answer right away, for a linpk solve.
c           UPDATE: convergence has a lot of trouble when the steps
c           are *too* big, so let's keep the source/drain voltages
c           reasonably small when we start - even for LINPK...
c
            if ( soltyp .eq. 'linpk' ) then
c              wVs = Vs
c              wVd = Vd
              wVs = sign( min( abs(Vs), 5.0), Vs )
              wVd = sign( min( abs(Vd), 5.0), Vd )
            else
              wVs = Vs
              wVd = Vd
c             wVs = Sphis + 1.414/8.0
c             wVd = Dphis + 1.414/8.0
            endif
          endif
        endif
c
c       Set these to something big, it's trivial
c
        dVmax = 1.0e+32
        dVs = 0.0
        dVd = 0.0
        dVg1 = 0.0
        dEavg = 0.0
c
c       Set the final gate biases because our guess will be good for
c       this.
c
        wVg1 = Vg1
        wEavg = eavg
c
c       We need the initial deltas, but since Vg1 is at its final
c       value, we only need to move the phin's. For the center simulation
c       these are decreases from Phis or Phim to the operating
c       point.
c
        dVb = delV
        dE = delE
c
c       Start with the best guess we can muster, and be vocal
c
        call doguess( wVs, wVd, wVg1, wEavg, .false., .true., error )
c
c       ...and jump right into the biasing loop
c
        goto 101

      elseif ( guess .eq. 'charge' ) then
c
c       Good, a charge neutral guess. so apply the equilibrium biases
c       and ramp to the final values
c
c				start at zero only if full simulation
        if ( symtyp .eq. 'full' ) then
          wVs = 0.0
          wVd = 0.0
        else
          wVs = Vs
          wVd = Vd
        endif
        wVg1 = 0.0
        wEavg = 0.0
c				just set it to something big, it's trivial
        dVmax = 1.0e+32
        dVs = 0.0
        dVd = 0.0
        dVg1 = 0.0
        dEavg = 0.0
c				we need the initial deltas
        dVb = delV
        dE = delE
c
c       We need to use SOR to get equilibrium
c
        call doguess( wVs, wVd, wVg1, wEavg, .true., .false., error )
        write(6,'('' '')')
        write(6,'('' '',''Now attempting to use equil() to '',
     +        ''solve for the equilibrium condition with SOR.'')')
        call equil( error )
c
c       Check for errors from the equilibrium SOR solution
c
        if ( error ) then
          write(6,'('' '',''**** Error! ****''/
     +        '' '',''     The equilibrium bias point could not be'',
     +            '' achieved!'')')
#         ifndef Cyber
            call flush(6)
#         endif
          call depmap( 1,1 )
          error = .true.
          return
        endif
      endif
c
c     Now into the general biasing loop
c
  100 call doguess( wVs, wVd, wVg1, wEavg, .true., .false., error )
c
c     Save the initial guess for the curious
c
  101 if ( debug .eq. 'yes' ) then
        call dump( m1psec, .true., zero )
      endif
c
c     Save the old values of the voltage and n for scaling later
c
      do 200 k=1,ynodes
        do 200 i=1,xnodes
          nold(i,k) = n(i,k)
          vold(i,k) = v(i,k)
  200 continue
c
c     Now call the solution system based on the initial condition
c     choice
c
      if ( icchan .eq. 'zero' ) then
        write(6,'('' '')')
        write(6,'('' '',''Now attempting to use newton() to '',
     +                  '' solve for the biases:'')')
        if ( symtyp .ne. 'center' ) then
          write(6,'('' '',''  Vs = '',f7.3,''   Vd = '',f7.3,
     +           ''   Vg1 = '',f7.3,''   Eavg = '',f7.1)')
     +           wVs, wVd, wVg1, wEavg
        else
          write(6,'('' '',''  phiNl = '',f7.3,''   phiNr = '',
     +           f7.3,''   Vg1 = '',f7.3,''   Eavg = '',f7.1)')
     +           wVs, wVd, wVg1, wEavg
        endif
#       ifndef Cyber
          call flush(6)
#       endif
  300   call newton( error, iter )
c
c       check to see if we have coupling problems
c
        if ( iter .eq. -1 ) then
          call refine( .false., .true., .true., error )
          if ( error ) then
            write(6,'('' '',''**** Error! ****''/
     +                '' '',''     Shark could not resolve the '',
     +                      ''coupling problems. This problem is '',
     +                      ''not solvable''/
     +                '' '',''     with the current executable of '',
     +                      ''Shark. The final output is being '',
     +                      ''written and''/
     +                '' '',''     the solution halted.'')')
            write(9,'('' '',''#eor'')')
#           ifndef Cyber
              call flush(6)
              call dump( m1psec, .true., zero )
              close( 9 )
#           endif
            call depmap( 1,1 )
            return
          endif
c
c         try the bias point again...
c
          call mkdev( .false. )
          call doguess( wVs, wVd, wVg1, wEavg, .false., .false., error )
          if ( error ) return
          write(6,'('' '',''...Re-Trying the bias point...'')')
          goto 300
        endif
      elseif ( (icchan .eq. 'DC') .or. (icchan .eq. 'flat') ) then
        write(6,'('' '')')
        write(6,'('' '',''Now attempting to use DC() '',
     +              ''to solve for the biases:'')')
        if ( symtyp .ne. 'center' ) then
          write(6,'('' '',''  Vs = '',f7.3,''   Vd = '',f7.3,
     +              ''   Vg1 = '',f7.3,''   Eavg = '',f7.1)')
     +              wVs, wVd, wVg1, wEavg
        else
          write(6,'('' '',''  phiNl = '',f7.3,''   phiNr = '',
     +              f7.3,''   Vg1 = '',f7.3,''   Eavg = '',f7.1)')
     +              wVs, wVd, wVg1, wEavg
        endif
#       ifndef Cyber
          call flush(6)
#       endif
        call DC( wVs, wVd, wVg1, wEavg, error )
      else
        write(6,'('' '')')
        write(6,'('' '',''Now attempting to use newton() to '',
     +              ''solve for the biases:'')')
        if ( symtyp .ne. 'center' ) then
          write(6,'('' '',''  Vs = '',f7.3,''   Vd = '',f7.3,
     +              ''   Vg1 = '',f7.3,''   Eavg = '',f7.1)')
     +              wVs, wVd, wVg1, wEavg
        else
          write(6,'('' '',''  phiNl = '',f7.3,''   phiNr = '',
     +              f7.3,''   Vg1 = '',f7.3,''   Eavg = '',f7.1)')
     +              wVs, wVd, wVg1, wEavg
        endif
#       ifndef Cyber
          call flush(6)
#       endif
  400   call newton( error, iter )
c
c       check to see if we had coupling problems
c
        if ( iter .eq. -1 ) then
          call refine( .false., .true., .true., error )
          if ( error ) then
            write(6,'('' '',''**** Error! ****''/
     +                '' '',''     Shark could not resolve the '',
     +                      ''coupling problems. This problem is '',
     +                      ''not solvable''/
     +                '' '',''     with the current executable of '',
     +                      ''Shark. The final output is being '',
     +                      ''written and''/
     +                '' '',''     the solution halted.'')')
            write(9,'('' '',''#eor'')')
#           ifndef Cyber
              call flush(6)
              call dump( m1psec, .true., zero )
              close( 9 )
#           endif
            call depmap( 1,1 )
            return
          endif
c
c         try the bias point again...
c
          call mkdev( .false. )
          call doguess( wVs, wVd, wVg1, wEavg, .false., .false., error )
          if ( error ) return
          write(6,'('' '',''...Re-Trying the bias point...'')')
          goto 400
        endif
      endif
c
c     check for errors from the system solution
c
      if ( error ) then
c
c       see if we are at the lower limit for dVb, and if not then
c       halve dVb and try it all again
c
        if ( (dVmax .gt. 0.001) .and.
     +       (max(dVs, dVd, dVg1, dEavg) .gt. 0.0) ) then
          write(6,'('' '',''This bias point could not be achieved,'',
     +              '' so we are trying a smaller bias change,'')')
#         ifndef Cyber
            call flush(6)
#         endif
c				remove half the bias changes
          dVs = dVs/2.0
          dVd = dVd/2.0
          dVg1 = dVg1/2.0
          dEavg = dEavg/2.0
          wVs = wVs - dVs
          wVd = wVd - dVd
          wVg1 = wVg1 - dVg1
          wEavg = wEavg - dEavg
c				the new dVmax is halved
          dVmax = dVmax/2.0
c				...as are the deltas used to get here
          dVb = dVb/4.0
          dE = dE/4.0
c				restore the old voltages
          do 500 k=1,ynodes
            do 500 i=1,xnodes
              n(i,k) = nold(i,k)
              v(i,k) = vold(i,k)
  500     continue
c				go back and try again without a guess
          goto 100
c
c       See if we are at the stopping point, and could not converge
c       due to possible oscillations.
c
        elseif ( abs(wVs-Vs)/kT + abs(wVd-Vd)/kT + abs(wVg1-Vg1)/kT +
     +           abs(wEavg-Eavg) .lt. vmax ) then
          write(6,'('' '',''Shark had trouble achieving the final'',
     +         '' operating point, but this''/
     +              '' '',''could be due to possible oscillations.'',
     +         '' Therefore, Shark is continuing,''/
     +              '' '',''but care should be taken in checking'',
     +         '' the results.'')')
          error = .false.
          return
        else
c
c         Nope, just an error
c
          write(6,'('' '',''**** Error! ****''/
     +              '' '',''     The final bias point could not be'',
     +        '' achieved! The trouble occured for:'')')
          write(6,'('' '',''     Vs = '',f7.3,''   Vd = '',f7.3,
     +              ''   Vg1 = '',f7.3,''   Eavg = '',f7.1)')
     +              wVs, wVd, wVg1, wEavg
#         ifndef Cyber
            call flush(6)
#         endif
          call depmap( 1,1 )
          error = .true.
          return
        endif
      endif
c
c     check to see if we are at the final values, and return to
c     the caller if we are.
c
      if ( abs(wVs-Vs)/kT + abs(wVd-Vd)/kT + abs(wVg1-Vg1)/kT +
     +     abs(wEavg-eavg) .lt. vmax ) then
        error = .false.
        return
      endif
c
c     we need to ramp a little more... so get the new values
c
      dVold = dVmax
      dVmax = -1.0
      if ( abs( Vs - wVs ) .gt. 0.0 ) then
        dVs = sign( min( abs(dVb), abs(Vs-wVs) ), (Vs - wVs) )
        wVs = wVs + dVs
        dVmax = max( abs(dVs), dVmax )
      else
        dVs = 0.0
      endif

      if ( abs( Vd - wVd ) .gt. 0.0 ) then
        dVd = sign( min( abs(dVb), abs(Vd-wVd) ), (Vd - wVd) )
        wVd = wVd + dVd
        dVmax = max( abs(dVd), dVmax )
      else
        dVd = 0.0
      endif

      if ( abs( Vg1 - wVg1 ) .gt. 0.0 ) then
        dVg1 = sign( min( abs(dVb), abs(Vg1-wVg1) ), (Vg1 - wVg1) )
        wVg1 = wVg1 + dVg1
        dVmax = max( abs(dVg1), dVmax )
      else
        dVg1 = 0.0
      endif

      if ( abs( eavg - wEavg ) .gt. 0.0 ) then
        dEavg = sign( min( abs(dE), abs(eavg-wEavg) ), (eavg - wEavg) )
        wEavg = wEavg + dEavg
        dVmax = max( abs(dEavg*length), dVmax )
      else
        dEavg = 0.0
      endif
c
c     now we need to move V up to be a good guess for the next
c     stab at the bias point
c     NOTE: this was a nice idea, but it's really a very bad
c           implementation. The optimal solution might be to
c           use three points and not just two. The problem is
c           that the values at each node don't scale with the
c           bias change and so it turns out to be faster to just
c           use the last solution as the initial guess of the
c           new bias condition. It's not as good as it probably
c           could be, but it's faster than anything else I've
c           tried.
c
c      if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
c        write(6,'('' '','' ...scaling n and v by '',f7.3)')
c     +          dVmax/dVold
c      endif
c      do 600 k=1,ynodes
c        do 600 i=1,xnodes
c          n(i,k) = n(i,k) + 0.1*( n(i,k)-nold(i,k) )*dVmax/dVold
c          v(i,k) = v(i,k) + 0.1*( v(i,k)-vold(i,k) )*dVmax/dVold
c  600 continue
c
c     move dVb and dE up so that we are accelerating the ramp
c
      dVb = sign( min( 2.0*abs(dVb), half ), dVb )
      dE = sign( min( 2.0*abs(dE), hund4 ), dE )
c
c     return and do it all again
c
      goto 100
c
c     all done...
c
      end
