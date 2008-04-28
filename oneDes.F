      subroutine oneDes( jline, fixn, phin, error, verbose, iter )
c
c     This subroutine solves the system for a consistent V
c     assuming phin, and p held in equilibrium, using Newton
c     iteration. The Jacobian is set up, the RHS is set up, and
c     the change to V is then calculated by the tri-diagonal
c     matrix solver. This is a ONE DIMENSIONAL solver routine.
c
c     This solver assumes that the insulator is a perfect one.
c     The potential barrier between the semiconductor and the
c     insulator is infinitely high. The electron concentration
c     in the insulator is held at zero.
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
c     now for the arguments to this routine
c
      integer jline, iter
      logical fixn, error, verbose
      PRECISION phin
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
c
c     now for the function declarations I need
c
      real dtime
c
c     now for some local variables
c
      integer i, inorm, jnorm, info
      real tarray(2), ittime, ntime
      PRECISION vl, vr, vij, hl, hr, nij, pij, norm
      PRECISION jl(Ny), jm(Ny), ju(Ny), f(Ny)
c
c     reset the number of iterations
c
      iter = 0
c
c     reset the elapsed timer for the start of the first pass
c
      ittime = dtime(tarray)
      ntime = 0.0
c
c     start a pass, and reset the norm of delta v for this pass
c 
  100 norm = -1.0
c
c     For every node in the system, collect the data needed and fill
c     in [J] and -F so that we can solve them to get the corrections
c     and, in the end, the result
c
      do 200 i=1,ynodes
c
c       collect the variables I will need for the simplified equations
c
c       get hl, hr (use symmetry at the boundaries)
c
        if ( i .eq. 1 ) then
          hl = hy(i)
        else
          hl = hy(i-1)
        endif

        if ( i .eq. ynodes ) then
          hr = hy(i-1)
        else
          hr = hy(i)
        endif
c
c       get vr, vl, vij (using boundary condition vectors)
c
        if ( i .eq. 1 ) then
c				only specify the gate voltage where
c				it is not there, otherwise it is the
c				same as where we are
          if ( vgate(jline) .eq. 1.0e20 ) then
            vl = v(jline,i)
          else
            vl = vgate(jline)
          endif
        else
          vl = v(jline, (i-1))
        endif

c				LOS at the back
        if ( i .eq. ynodes ) then
          vr = v(jline, (i-1))
        else
          vr = v(jline, (i+1))
        endif

        vij = v(jline, i)
c
c       Get nij based on the status of fixn
c
        if ( i .le. joxide ) then
          nij = 0.0
        else
          if ( fixn ) then
            nij = n(jline, i)
          else
            if ( symtyp .eq. 'Mock' ) then
              nij = exp( vij + log(dopsrc/nio) - phin )
            else
              nij = exp( vij - phin )
            endif
          endif
        endif
c
c       Get pij based on the offset in voltage
c
        if ( symtyp .eq. 'Mock' ) then
          pij = exp( -vij - log(dopsrc/nio) )
        else
          pij = exp( -vij )
        endif
c
c       Now I can fill in the partial w.r.t. vl in [J]
c
        if ( i .eq. 1 ) then
          jl(i) = 0.0
        elseif ( i .eq. joxide ) then
          jl(i) = 2.0/(hl*hr*(1.0 + (ksemi/kox)*(hl/hr)))
        elseif ( i .eq. ynodes ) then
          jl(i) = 2.0/(hl*(hl+hr)) + 2.0/(hr*(hl+hr))
        else
          jl(i) = 2.0/(hl*(hl+hr))
        endif
c
c       Now I can fill in the partial w.r.t. vij in [J]
c
        if ( i .lt. joxide ) then
          jm(i) = -2.0/(hl*hr)
        elseif ( i .eq. joxide ) then
          jm(i) = -2.0/(hl*hr)
        else
          if ( fixn ) then
            jm(i) = -2.0/(hl*hr) - pij
          else
            jm(i) = -2.0/(hl*hr) - pij - nij
          endif
        endif
c
c       Now I can fill in the partial w.r.t. vr in [J]
c
        if ( i .eq. ynodes ) then
          ju(i) = 0.0
        elseif ( i .eq. joxide ) then
          ju(i) = 2.0/(hr*(hl+hr))*(1.0 + (hl/hr))/((kox/ksemi) +
     +                                              (hl/hr))
        else
          ju(i) = 2.0/(hr*(hl+hr))
        endif
c
c       Now I can fill in the section of -F for this node
c
        if ( i .lt. joxide ) then
c				-Fv for this node in the oxide
          f(i) = -( (hr*vl-(hl+hr)*vij+hl*vr)/(0.5*hr*hl*(hl+hr)) )
        elseif ( i .eq. joxide ) then
c				-Fv for this node at the interface
          f(i) = -( 2.0*Qss(jline)/(nio*hr*Ldi) +
     +              (2.0/(hl*hr*(hl+hr)))*( vl*(hl+hr)/(1.0+
     +              (ksemi/kox)*(hl/hr)) - vij*(hl+hr) + vr*hl*
     +              (1.0+hl/hr)/(kox/ksemi + hl/hr) ) )
        else
c				-Fv for this node in the semi.
          f(i) = -( pij + dop(jline,i) - nij +
     +              (hr*vl-(hl+hr)*vij+hl*vr)/(0.5*hr*hl*(hl+hr)) )
        endif
c
c       Now normalize the rows of [J] w.r.t. the diagonal
c       elements to check for diagonal dominance
c
        if ( jm(i) .ne. 0.0 ) then
          ju(i) = ju(i)/jm(i)
          jl(i) = jl(i)/jm(i)
          f(i) = f(i)/jm(i)
          jm(i) = 1.0
        endif
c
c     All done filling out [J] and -F!
c
  200 continue
#     ifdef LAPACK
c
c     We need to shift jl(i) down one because LAPACK doesn't want
c     to see the unnecessary 0 in the first spot. So do a simple
c     shift.
c
      do 300 i=1,ynodes-1
        jl(i) = jl(i+1)
  300 continue
c
c     Now we can call the solution system for this mess
c
      call dgtsv( ynodes, 1, jl, jm, ju, f, ynodes, info )
      if ( info .ne. 0 ) then
        write(6,'('' '',''***** Error in dgtsv! error #'',i3,
     +                  '' *****''/
     +            '' '',''      This system cannot be solved!'')')
     +          info
        error = .true.
        return
      endif
#     else
c
c     Now we can call the solution system for this mess
c
      call sgtsl( ynodes, jl, jm, ju, f, info )
      if ( info .ne. 0 ) then
        write(6,'('' '',''***** Error in sgtsl! error #'',i3,
     +                  '' *****''/
     +            '' '',''      This system cannot be solved!'')')
     +          info
        error = .true.
        return
      endif
#     endif
c
c     Now I need to go through the system again, updating the variables
c     that I have the corrections for by using Brown and Lindsay's 
c     scheme
c
      do 400 i=1,ynodes
c
c       now change f(i) to go along with Brown and Lindsay's
c       convergence criteria for the voltage
c
        if ( (abs(f(i)) .gt. 1.0) .and. (abs(f(i)) .le. 3.7) ) then
          f(i) = sign( abs(f(i))**0.2, f(i) )
        elseif ( abs(f(i)) .gt. 3.7 ) then
          f(i) = sign( log(abs(f(i))), f(i) )
        endif
c
c       update the voltage matrix, and the electron concentration too
c
        v(jline,i) = v(jline,i) + f(i)
        if ( .not. fixn ) then
          if ( i .le. joxide ) then
            n(jline,i) = 0.0
          else
            if ( symtyp .eq. 'Mock' ) then
              n(jline,i) = exp( v(jline,i) + log(dopsrc/nio) - phin )
            else
              n(jline,i) = exp( v(jline,i) - phin )
            endif
          endif
        endif
c
c       save the maximum f(i) for convergence tests, and the
c       location and type of the maximum update
c
        if ( abs(f(i)) .gt. norm ) then
          norm = abs(f(i))
          inorm = jline
          jnorm = i
        endif
c
c     We are all done updating the system
c
  400 continue
c
c     well, we have done another pass, so count it
c
      iter = iter + 1
c
c     get the time for this iteration
c
      ittime = dtime(tarray)
      ntime = ntime + ittime
c
c     do the debugging stuff if need be
c
      if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
c
c       output, for now, the statistics for this pass
c
        if ( verbose ) then
          write(6,'('' '',''Iter. #'',i3,'' took '',f6.2,'' msec. '',
     +                    ''with norm = '',g9.3,'' at V('',i3,'','',
     +                    i3,'').'')') 
     +            iter, 1.0e3*ittime, norm, inorm, jnorm
#         ifndef Cyber
            call flush(6)
#         endif
        endif
      endif
c
c     check iteration count, and add error message, if need be
c
      if ( iter .gt. (2*maxit) ) then
        if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
          write(6,'('' '')')
          write(6,'('' '',''*** Error! The solution failed to '',
     +                    ''converge in '',i4,'' Newton '',
     +                    ''iterations!'')') 2*maxit
          write(6,'('' '')')
          write(6,'('' '',''This unsuccessful solution took a total '',
     +                    ''of '',i5,'':'',i2.2,''.'')') 
     +           int(ntime/60.0), int(ntime-int(ntime/60.0)*60)
#         ifndef Cyber
            call flush(6)
#         endif
        endif
        error = .true.
        return
      endif
c
c     test for convergence of this pass
c
      if ( norm .gt. 500.0 ) then
        if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
          write(6,'('' '')')
          write(6,'('' '',''*** Error! The solution has an '',
     +                    ''adjustment of '',f6.2,'' kT/q!''/
     +              '' '',''    This is divergent, so the solution '',
     +                    ''is being halted!'')') norm
          write(6,'('' '',''This unsuccessful solution took a total '',
     +                    ''of '',i5,'':'',i2.2,''.'')') 
     +             int(ntime/60.0), int(ntime-int(ntime/60.0)*60)
          write(6,'('' '')')
#         ifndef Cyber
            call flush(6)
#         endif
        endif
        error = .true.
        return
      endif
      if ( norm .gt. vmax ) goto 100
c
c     write out a success message
c
      if ( verbose ) then
        write(6,'('' '',''   The Newton solution converged in '',i4,
     +                  '' iterations which took '',f5.2,'' sec.'')')
     +               iter, ntime
#       ifndef Cyber
          call flush(6)
#       endif
      endif
c
c     all done... so return to the caller
c
      error = .false.
      return
      end
