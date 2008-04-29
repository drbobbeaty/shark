      subroutine equil( error )
c
c     This subroutine solves the system for an equilibrium set of
c     V and n (assuming p is held in equilibrium), using SOR.
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
      PRECISION omega
      parameter( Nx=XMAX, Ny=YMAX, omega=1.7 )
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
c
c     now for some local variables
c
      logical LOSl, LOSr
      integer i, k, iter, inorm, jnorm
      real t0, ittime, ntime
      PRECISION vl, vr, vt, vb, vij, hl, hr, ht, hb, norm,
     +          Vob, Voc, dV
c
c     precalculate the bulk potential
c
      if ( chmod .eq. 'buried' ) then
        Voc = sign( log( (abs(dopch) + sqrt(dopch*dopch +
     +              4.0*nio*nio))/(2.0*nio) ), dopch )
      else
        Voc = sign( log( (abs(dopsub) + sqrt(dopsub*dopsub +
     +              4.0*nio*nio))/(2.0*nio) ), dopsub )
      endif
      Vob = -sign( log( (abs(dopsub) + sqrt(dopsub*dopsub +
     +              4.0*nio*nio))/(2.0*nio) ), dopsub )
c
c     reset the number of iterations
c
      iter = 0
c
c     reset the elapsed timer for the start of the first pass
c
      ntime = 0.0
c
c     start a pass, and reset the norm of delta v for this pass
c
  100 norm = 0.0
c
c     get a mark as to the starting time of this pass - we'll diff
c     the other times in this pass against this for elapsed times.
c
      call cpu_time(t0)
c
c     For every node in the system, collect the data needed and
c     calculate the corrections to v(i,k) and n(i,k) and, in the
c     end, the result
c
      do 200 k=1,ynodes
        do 200 i=1,xnodes
c
c         collect the variables I will need for the simplified equations
c
c         get hl, hr, ht, hb (use symmetry at the boundries)
c
          if ( i .eq. 1 ) then
            hl = hx(i)
          else
            hl = hx(i-1)
          endif

          if ( i .eq. xnodes ) then
            hr = hx(i-1)
          else
            hr = hx(i)
          endif

          if ( k .eq. 1 ) then
            ht = hy(k)
          else
            ht = hy(k-1)
          endif

          if ( k .eq. ynodes ) then
            hb = hy(k-1)
          else
            hb = hy(k)
          endif
c
c         get vr, vl, vt, vb, vij (using boundary condition vectors)
c
          if ( i .eq. 1 ) then
c
c           Check for a LOS or a fixed boundary condition
c
            if ( srccon(k,3) .eq. 1.0 ) then
              vl = srccon(k,1)
              LOSl = .false.
            else
              vl = v((i+1),k)
              LOSl = .true.
            endif
          else
            vl = v((i-1),k)
          endif

          if ( i .eq. xnodes ) then
c
c           Check for a LOS or a fixed boundary condition
c
            if ( drncon(k,3) .eq. 1.0 ) then
              vr = drncon(k,1)
              LOSr = .false.
            else
              vr = v((i-1),k)
              LOSr = .true.
            endif
          else
            vr = v((i+1),k)
          endif

          if ( k .eq. 1 ) then
c				only specify the gate voltage where it is
c				not there, otherwise it is the same as
c				where we are.
            if ( vgate(i) .eq. 1.0e20 ) then
              vt = v(i,k)
            else
              vt = vgate(i)
            endif
          else
            vt = v(i,(k-1))
          endif

          if ( k .eq. ynodes ) then
            vb = v(i,(k-1))
          else
            vb = v(i,(k+1))
          endif

          vij = v(i,k)
c
c         Now I can calculate the size of dV
c
          if ( k .lt. joxide ) then
c				dV for this node in the oxide
            dV = (hl*hr*ht*hb/(hl*hr+ht*hb))*(
     +             (hr*vl+hl*vr)/(hl*hr*(hl+hr)) +
     +             (ht*vb+hb*vt)/(ht*hb*(ht+hb)) ) - vij
          elseif ( k .eq. joxide ) then
c				dV for this node at the interface
            dV = (hl*hr*ht*hb/(hl*hr+ht*hb))*(
     +            2.0*Qss(i)/(nio * hb * Ldi) +
     +            (hr*vl+hl*vr)/(hl*hr*(hl+hr)) +
     +            ( vt*(ht+hb)/(1.0+(ksemi/kox)*(ht/hb)) +
     +              vb*ht*(1.0+ht/hb)/(kox/ksemi+ht/hb) )/
     +             (ht*hb*(ht+hb)) ) - vij
          else
c				dV for this node in the semi
            dV = (hl*hr*ht*hb/(hl*hr+ht*hb))*(
     +             exp(-vij) + dop(i,k) - exp(vij) +
     +             (hr*vl+hl*vr)/(hl*hr*(hl+hr)) +
     +             (ht*vb+hb*vt)/(ht*hb*(ht+hb)) ) - vij
          endif
c
c         update the voltage matrix and the electron concentration,
c         with a bit of overrelaxation
c
          if ( abs(dV) .gt. 1.0 ) then
            v(i,k) = v(i,k) + omega*dV
          else
            v(i,k) = v(i,k) + dV
          endif
          if ( k .le. joxide ) then
            n(i,k) = 0.0
          else
            n(i,k) = exp( v(i,k) )
          endif
c
c         save the maximum of dV for convergence tests, and the
c         location.
c
          if ( abs(dV) .gt. norm ) then
            norm = abs(dV)
            inorm = i
            jnorm = k
          endif
c
c     We are all done updating the system for one pass of SOR
c
  200 continue
c
c     Well, we have another pass, so count it
c
      iter = iter + 1
c
c     get the time for this iteration
c
      call cpu_time(ittime)
      ittime = ittime - t0
      ntime = ntime + ittime
c
c     do the debugging stuff if need be
c
      if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
c
c       output, for now, the statistics for this pass
c
        write(6,'('' '',''Iter. #'',i4,'' took '',f5.2,'' msec.'',
     +       '' with norm = '',g9.3,'' at V('',i3,'','',i3,'').'')')
     +       iter, 1.0e3*ittime, norm, inorm, jnorm
#       ifndef Cyber
          call flush(6)
#       endif
c
c       print out the depletion map finely, for debugging purposes,
c       stopping early if we want to
c
        if ( iter .eq. ((iter/(5*maxit/5))*(5*maxit/5)) ) then
          call depmap( 1,1 )
        endif
      endif
c
c     check iteration count and add an error message, if need be
c
      if ( iter .gt. 5*maxit ) then
        if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
          write(6,'('' '')')
          write(6,'('' '',''*** Error! *** The solution failed to '',
     +       ''converge in '',i4,'' iterations!'')') 5*maxit
          write(6,'('' '')')
c				output the entire map for debugging
          call depmap( 1,0 )
          write(6,'('' '',''This unsuccessful solution took a total '',
     +        ''of '',f7.2,'' msec.'')') 1.0e3*ntime
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
      if ( norm .gt. 600.0 ) then
        if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
          write(6,'('' '')')
          write(6,'('' '',''*** Error! The solution has an '',
     +                    ''adjustment of '',f6.2,'' kT/q!''/
     +              '' '',''This is divergent! So the solution is '',
     +                    ''being halted!'')') norm
          write(6,'('' '',''This unsuccessful solution took a total '',
     +        ''of '',f7.2,'' msec.'')') 1.0e3*ntime
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
c     write out  a success message
c
      write(6,'('' '',''   The equilibrium SOR solution '',
     +      ''converged in'',i4,'' iterations which took '',f5.2,
     +      '' msec.'')') iter, 1.0e3*ntime/60.0
#     ifndef Cyber
        call flush(6)
#     endif
c
c     output the charge distribution for this equilibrium solution
c
      if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
        write(6,'('' '')')
        write(6,'('' '')')
        write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                            REB 1988'',
     +     ''         charge map ===='')')
        write(6,'('' '')')
        write(6,'('' '','' The following is the charge map for '',
     +              ''equilibrium:'')')
        write(6,'('' '',''    g = resistive gate'')')
        write(6,'('' '',''    i = insulator, or oxide'')')
        write(6,'('' '',''    n = n-type'')')
        write(6,'('' '',''    p = p-type'')')
        write(6,'('' '',''    blank = depletion region'')')
        write(6,'('' '')')
        write(6,'('' '')')
c
c       call the depletion mapper with a fine set, and stop early
c
        call depmap( 1,1 )
      endif
c
c     all done... so return to the caller
c
      error = .false.
      return
      end
