      subroutine doguess( gVs, gVd, gVg1, gEavg, BConly,
     +                    verbose, error )
c
c     This subroutine computes the initial guess for the simulation.
c     This entails assuming MOS electrostatics under the gate, and
c     pn junction electrostatics under the drain and source.
c
c     The initial condition on the electron concentration in the
c     channel is also specified from MOS inversion-layer equations.
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
c     This is the percentage error in Qn allowed in the calls to
c     Qn2rho() for the running pulse
c
      PRECISION Qtol
      parameter( Qtol = 0.0001 )
c
c     now for the arguments to this routine
c
      PRECISION gVs, gVd, gVg1, gEavg
      logical BConly, verbose, error
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
      common /Pshape/ Qm, ftw, fwhm, pulse
      PRECISION Qm, ftw, fwhm
      logical pulse
c
c     now for some local variables
c
      integer i, j, i1, j1, iter, jsrc
      PRECISION y, Vsbi, Vdbi, xsp, xdp, Cox, Vo, phis, xd,
     +          Vbulk, ytemp, Vdo, Vso, Vob, Voc, V1, Qn, vtemp,
     +          ysrc, compVs, compVd, Qntemp
      PRECISION rho(Nx), Vr(Ny), Vl(Ny), nr(Ny), nl(Ny)
      logical gate(Nx)
c
c     ...and now for some functions that need to be defined
c
      PRECISION Qn2rho
c
c     Start with a header so the output looks reasonable
c
      if ( verbose ) then
        write(6,'('' '')')
        write(6,'('' '')')
        write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                            REB 1988'',
     +     ''      initial guess ===='')')
        write(6,'('' '')')
      endif
c
c     First, we need the initial guess on v(i,j) and n(i,j)
c     to start the iteration process
c
c     reset the variables in the simulation region only if we are
c     doing more than the boundary conditions
c
      if ( .not. BConly ) then
        do 100 i=1,xnodes
  100     gate(i) = .true.
        do 150 j=1,ynodes
          do 150 i=1,xnodes
            v(i,j) = 0.0
            n(i,j) = 0.0
c				also set a flag for the type of guess
            if ( dop(i,j) .gt. 0.0 ) gate(i) = .false.
  150   continue
      endif
c
c     If we are doing a center, or Mock, simulation then we had better
c     be doing a MOS initial guess.
c
      if ( (symtyp .ne. 'full') .and. (guess .ne. 'MOS') ) then
        if ( verbose ) then
          write(6,'('' '',''***** Warning! *****''/
     +              '' '',''      For a center simulation you '',
     +        ''should really be doing an MOS initial guess.''/
     +              '' '',''      Shark is doing this for you now.'')')
        endif
        guess = 'MOS'
      endif
c
c     Now I need to see which type of initial guess I am really
c     going to be computing
c
      if ( guess .eq. 'zero' ) goto 3400
      if ( guess .eq. 'charge' ) goto 2900
c
c     Now see if we are simulating the entire device, or just the center
c
      if ( symtyp .eq. 'center' ) goto 1200
c
c     First, set the biases at the gate, source and drain for use
c     later.
c
c     start by resetting the temporary variables, and then the arrays
c
      do 160 j=1,ynodes
        srccon(j,3) = 0.0
        drncon(j,3) = 0.0
  160 continue
      j = 1
      ytemp = hy(1) * Ldi
c
c     correct the new drain and source voltages by the bulk
c     contact potentials, as well as starting the gate voltage
c     correctly
c
      Vbulk = -sign( kT*log( abs(dopsub/nio)/2.0 + sqrt( (dopsub*
     +        dopsub/(nio*nio) + 4.0)/4.0 )), dopsub )
      Vso = sign( kT*log( abs(dopsrc/nio)/2.0 + sqrt( (dopsrc*
     +        dopsrc/(nio*nio) + 4.0)/4.0 )), dopsrc )
      Vdo = sign( kT*log( abs(dopdrn/nio)/2.0 + sqrt( (dopdrn*
     +        dopdrn/(nio*nio) + 4.0)/4.0 )), dopdrn )
      if ( dopsub .lt. 0.0 ) then
        Vsbi = abs( kT*log( abs(dopsrc/dopsub) ) )
        Vdbi = abs( kT*log( abs(dopdrn/dopsub) ) )
      elseif ( dopsub .gt. 0.0 ) then
        Vsbi = abs( kT*log( dopsrc*dopsub/(nio*nio) ) )
        Vdbi = abs( kT*log( dopdrn*dopsub/(nio*nio) ) )
      else
        Vsbi = abs( kT*log( dopsrc/nio ) )
        Vdbi = abs( kT*log( dopdrn/nio ) )
      endif
      xsp = sqrt( ((2.0*ksemi*8.854e-14)/1.6022e-19)*
     +            (Vsbi - (0.0-(gVs-Vso)) )/dopsub )
      xdp = sqrt( ((2.0*ksemi*8.854e-14)/1.6022e-19)*
     +            (Vdbi - (0.0-(gVd-Vdo)) )/dopsub )
      Cox = kox * 8.854e-14/tox
      Vo = 1.6022e-19*ksemi*8.854e-14*dopsub/(Cox*Cox)
c
c     Get the potentials for the boundary conditions
c
      if ( symtyp .eq. 'full' ) then
        compVs = (gVs + Vso)/kT
        compVd = (gVd + Vdo)/kT
        vtemp = (gVg1 + Vbulk)/kT
      else
c			Sutherland's and Mock's BC on the source and drain
        compVs = (gVs/kT) - 1
        compVd = (gVd/kT) - 1
        vtemp = (gVg1 + Vbulk - Vso)/kT
      endif
c
c     Do the resistive gate that spans the entire
c     length of the active channel
c
      drncon(j,2) = dopdrn/nio
      do 170 i=1,xnodes
        vgate(i) = vtemp
        vtemp = vtemp + (gEavg * hx(i) * Ldi)/kT
  170 continue
c
c     start at the top and work our way down, first skip past the oxide
c
  180 ytemp = ytemp + hy(j)*Ldi
      j = j + 1
      if ( ytemp .le. tox ) goto 180
c
c     save the start of the semiconductor for the drain
c
      jsrc = j
      ysrc = ytemp
c
c     Put the contact voltage in the matrix. Check to see if we are
c     still on the contact.
c
  190 srccon(j,1) = compVs
      srccon(j,3) = 1.0
c				scan down the y direction
      ytemp = ytemp + hy(j)*Ldi
      j = j + 1
      if ( ytemp .le. (tox+tsorce) ) goto 190
c
c     reset the variables to the top of the semiconductor
c
      j = jsrc
      ytemp = ysrc
c
c     Put the contact voltage in the matrix. Check to see if we are
c     still on the contact.
c
  210 drncon(j,1) = compVd
      drncon(j,3) = 1.0
c				scan down the y direction
      ytemp = ytemp + hy(j)*Ldi
      j = j + 1
      if ( ytemp .le. (tox+tdrain) ) goto 210
c
c     If this is a Mock simulation, then we should finish it out
c     with the center simulation guess
c
      if ( symtyp .eq. 'Mock' ) goto 1200
c
c     If we are only to be doing the boundary conditions, then skip the
c     rest of this, and return to the caller.
c
      if ( BConly ) return
c
c     Now do the MOS-and-Junction initial guess
c
      if ( verbose ) then
        write(6,'('' '',''    Computing the initial guess based on'',
     +      '' the MOS-and-Junction electrostatics.'')')
      endif
c
c     For the initial guess we need a complex function broken into
c     three different regions in the p material and 2 different regions
c     in the N+ source and drains
c
c     output the values of the initial guess variables
c
      if ( verbose ) then
        write(6,'('' '')')
        write(6,2) Vbulk
        write(6,3) Vso
        write(6,4) Vdo
        write(6,5) Vsbi
        write(6,6) Vdbi
        write(6,7) xsp*1.0e4
        write(6,8) xdp*1.0e4
        write(6,9) Cox
        write(6,10) Vo
        write(6,'('' '')')
      endif
    2 format(' ','       Substrate contact potential = ',f6.4,' (V).')
    3 format(' ','          Source contact potential = ',f6.4,' (V).')
    4 format(' ','           Drain contact potential = ',f6.4,' (V).')
    5 format(' ','   Source to Bulk Built-in Voltage = ',f6.4,' (V).')
    6 format(' ','    Drain to Bulk Built-in Voltage = ',f6.4,' (V).')
    7 format(' ','    Source to Bulk Depletion Width = ',f5.2,' (um).')
    8 format(' ','     Drain to Bulk Depletion Width = ',f5.2,' (um).')
    9 format(' ','             Insulator Capacitance = ',g8.2,
     +           ' (F/cm**2).')
   10 format(' ','              Bulk-Charge MOS Term = ',g11.4,' (V).')
c
c     scan left to right across the entire device
c
      if ( verbose ) write(6,11)
   11 format(' ','   compute the MOS initial guess at the solution...')
      do 1100 i=1,xnodes
c
c       start at the top of the interesting region
c
        y = 0.0
        j = 1
c
c       see if we are in a line where there is N-type doping
c
        if ( .not. gate(i) ) then
c
c         yes, we are on a line with at least some n-type material
c
c				now get past the oxide because I
c				will set it there later with the
c				seeding. make n=0.
  200     n(i,j) = 0.0
          y = y + hy(j) * Ldi
          j = j + 1
          if ( y .le. tox ) goto 200
c
c         see which side of the device we are on, either the source
c         or the drain.
c
          if ( i .lt. (xnodes/2) ) then
c				now the n source contact voltage
  300       v(i,j) = ( log(dop(i,j))/log(dopsrc/nio) )*(gVs/kT)
            n(i,j) = dop(i,j)
            y = y + hy(j) * Ldi
            j = j + 1
            if ( (dop(i,j) .ge. 0.0) .and. (j .le. ynodes) ) goto 300
            ytemp = y
c				now the p depletion voltage
            if ( xsp .ne. 0.0 ) then
  400         v(i,j) = ( (-Vbulk/(xsp*xsp))*(y - (xsp+ytemp))**2
     +                   + Vbulk )/kT
              y = y + hy(j) * Ldi
              j = j + 1
              if ( y .le. (ytemp + xsp) ) goto 400
            endif
          else
c				now the n drain contact voltage
  500       v(i,j) = ( log(dop(i,j))/log(dopdrn/nio) )*(gVd/kT)
            n(i,j) = dop(i,j)
            y = y + hy(j) * Ldi
            j = j + 1
            if ( (dop(i,j) .ge. 0.0) .and. (j .le. ynodes) ) goto 500
            ytemp = y
c				now the p depletion voltage
            if ( xdp .ne. 0.0 ) then
  600         v(i,j) = ( (-Vbulk/(xdp*xdp))*(y - (xdp+ytemp))**2
     +                   + Vbulk )/kT
              y = y + hy(j) * Ldi
              j = j + 1
              if ( y .le. (ytemp + xdp) ) goto 600
            endif
          endif
c				now the quasi-equilibrium voltage
          do 700 j1=j,(ynodes-j+1)
            v(i,j1) = Vbulk/kT
  700       n(i,j1) = nio*nio / abs(dop(i,j1))
        else
c
c         we are underneath the gate totally
c
          phis = (vgate(i)*kT-Vbulk+1.6022e-19*Qss(i)/Cox) + Vo
     +           - sqrt( Vo*Vo + 2.0*Vo*(vgate(i)*kT-Vbulk+
     +           1.6022e-19*Qss(i)/Cox) )
          xd = sqrt( 2.0*8.854e-14*ksemi*phis/(1.6022e-19*dopsub) )
c				now compute the voltage in the oxide
  800     v(i,j) = ( 1.6022e-19*dopsub*xd*( xd/(2.0*8.854e-14*ksemi)
     +              + (tox - y)/(8.854e-14*kox) ) + Vbulk )/ kT
          y = y + hy(j)*Ldi
          j = j + 1
          if ( y .le. tox ) goto 800
c				now the depletion voltage
  900     v(i,j) = ( (1.6022e-19*dopsub/(2.0*ksemi*8.854e-14))*
     +          ( y*y - 2.0*(tox+xd)*y + (tox+xd)**2 ) + Vbulk )/kT
          y = y + hy(j)*Ldi
          j = j + 1
          if ( y .le. (tox+xd) ) goto 900
c				now the quasi-equilibrium voltage
          do 1000 j1=j,(ynodes-j+1)
            v(i,j1) = Vbulk/kT
 1000       n(i,j1) = nio*nio / abs(dop(i,j))
        endif
 1100 continue
c
c     now by-pass the reset of the initial guess routines
c
      goto 3400
c
c
c     Now do the MOS initial guess for the center simulation
c
c
 1200 if ( verbose ) then
        write(6,'('' '',''    Computing the initial guess based on'',
     +     '' MOS electrostatics for the''/
     +            '' '',''    center section device simulation.'')')
      endif
c
c     get some constants that I will need
c
      Cox = kox * 8.854e-14 / tox
      Vo = 1.6022e-19 * ksemi * 8.854e-14 * dopsub / (Cox*Cox)
      V1 = (1.6022e-19*tepi/Cox)*(dopsub + dopch)
      if ( chmod .eq. 'buried' ) then
        Voc = sign( log( (abs(dopch) + sqrt(dopch*dopch +
     +              4.0*nio*nio))/(2.0*nio) ), dopch )
      else
        Voc = -sign( log( (abs(dopsub) + sqrt(dopsub*dopsub +
     +               4.0*nio*nio))/(2.0*nio) ), dopsub )
      endif
      Vob = -sign( log( (abs(dopsub) + sqrt(dopsub*dopsub +
     +             4.0*nio*nio))/(2.0*nio) ), dopsub )
c
c     Compensate Vob for Sutherland's and Mock's shift in the reference
c     potential
c
      if ( symtyp .eq. 'Mock' ) then
        Vob = Vob - sign( log( (abs(dopsrc) + sqrt(dopsrc*dopsrc +
     +                    4.0*nio*nio))/(2.0*nio) ), dopsrc )
      endif
c
c     If all we are doing is setting new boundary conditions
c     then save the necessary system variables for later replacement.
c
      if ( BConly ) then
        do 1300 j=1,ynodes
          Vl(j) = v(1,j)
          Vr(j) = v(xnodes,j)
          nl(j) = n(1,j)
          nr(j) = n(xnodes,j)
 1300   continue
      endif
c
c     Make sure that tepi=0 if we are doing a MODFET simulation
c
      if ( chmod .eq. 'MODFET' ) then
        tepi = 0.0
      endif
c
c     Set the gate bias to that which I want the guess based
c
      vtemp = gVg1/kT + Vob
      do 1400 i=1,xnodes
        vgate(i) = vtemp
        vtemp = vtemp + (gEavg * hx(i) * Ldi)/kT
 1400 continue
c
c     Check on the type of variation in the quasi-fermi level
c     we need depending on the conditions of the simulation
c
      if ( gEavg .eq. 0.0 ) then
        if ( chmod .eq. 'MODFET' ) then
c
c         Prepare rho(i) to be the quasi-fermi level for the i-th
c         slice, by looking to see if we can do Pao-Sah, or assume
c         a linear variation
c
          if ( verbose ) then
            write(6,'('' '',''    Using Pao-Sah to calculate the'',
     +         '' initial quasi-fermi levels...'')')
          endif
          call paosah( rho, gVs, gVd, gVg1, error )
        else
c
c         Do a simple polynomial fit to the data
c
          if ( verbose ) then
            write(6,'('' '',''    Using simple polynomial fit '',
     +         ''to calculate the initial quasi-fermi levels...'')')
          endif
          do 1500 i=1,xnodes
            rho(i) = ( (gVd-gVs)*((1.0*i)/(1.0*xnodes))**3 + gVs )/kT
 1500     continue
        endif
      else
c
c       See if we are running with a pulse
c
        if ( pulse ) then
c
c         Get the Qn's and keep them in rho(i) because we will end up
c         converting these into quasi-fermi levels
c
          do 1600 i=1,xnodes
            if ( abs(i - (xnodes/2)) .le. (ftw/(2.0*hx(1)*Ldi)) ) then
              rho(i) = Qm
            else
c				this is a Gaussian shaped pulse
c              rho(i) = Qm*exp( -2.0*( 1.0*abs(i-(xnodes/2)) - (ftw
c     +                 /(2.0*hx(1)*Ldi)) )**2.0 / (fwhm/
c     +                 (hx(1)*Ldi))**2.0 )
c				this is a Parabolic shaped pulse
              rho(i) = Qm*( 1.0 - 2.0*( ( abs(i-(xnodes/2)) - (ftw
     +                 /(2.0*hx(1)*Ldi)) )/(fwhm/(hx(1)*Ldi)) )**2.0 )
c				truncated at the zero limit for charge
              rho(i) = max( rho(i), (1.0/nio) )
            endif
 1600     continue
c
c         Convert the Qn's to quasi-fermi levels
c
          do 1700 i=1,xnodes
            Qntemp = rho(i)
            rho(i) = Qn2rho( Qntemp, i, Qtol, iter, error )
            if ( (verbose) .and. (debug .eq. 'yes') ) then
              write(6,'('' '',''Qn('',i3,'') = '',g10.4,
     +                  '' yields a phin = '',f8.4,
     +                  '' and took '',i2,'' search passes.'')')
     +             i, Qntemp, rho(i), iter
#             ifndef Cyber
                call flush(6)
#             endif
            endif
            if ( error ) return
 1700     continue
c
c       See if the operator wants a flat channel, or the DC solution
c
        elseif ( (icchan .eq. 'flat') .or. (icchan .eq. 'DC') ) then
c
c         Do the flat distribution, so we need some numbers
c
c				get the full 1D solution for phis, etc.
          call oneDes( 1, .false., gVs/kT, error, .false., iter )
          call full1D( 1, error, .false., iter )
c				get Qn based on this
          Qn = -1.6022e-19*Qss(i) - Cox*( (vgate(i)-Vob)*kT - 
     +         (v(1,joxide+1)-Vob)*kT - sqrt( 2.0*Vo*(v(1,joxide+1)
     +         -Vob)*kT ) )
c
c         Calculate all the rho()'s based on this
c
          do 2000 i=1,xnodes
            rho(i) = (gVs/kT) - ( v(1,joxide+1)-Vob ) +
     +         ( (vgate(i)-Vob)*kT + (Qn+1.6022e-19*Qss(i))/Cox + Vo
     +         - sqrt( Vo*Vo + 2.0*Vo*( (vgate(i)-Vob)*kT + (Qn +
     +         1.6022e-19*Qss(i))/Cox ) ) )/kT
 2000     continue
c
c       The default is the nice tail (icchan = 'lorenz')
c
        else
          do 2100 i=1,xnodes
 2100       rho(i) = ( gVs*(xnodes-i)+gVd*(i-1) )/(kT*(xnodes-1))
        endif
      endif
c
c     Now guess the left side voltage for prepping into oneDes()
c
c     Get the constants for this slice in x
c
      i = 1
      Cox = kox * 8.854e-14 / tox
      Vo = 1.6022e-19 * ksemi * 8.854e-14 * dopsub / (Cox*Cox)
      V1 = (1.6022e-19*tepi/Cox)*( dopsub + dopch )
      phis = ( vgate(i)-Vob )*kT + 1.6022e-19*Qss(i)/Cox + Vo + V1 -
     +       sqrt( Vo*Vo + 2.0*Vo*( ( vgate(i)-Vob )*kT + 1.6022e-19*
     +       Qss(i)/Cox+V1 ) + V1*1.6022e-19*tepi*dopsub/Cox )
      xd = sqrt( tepi*tepi*(1.0 + dopch/dopsub)
     +           + 2.0*ksemi*8.854e-14*phis/(1.6022e-19 * dopsub) )
c
c     start at the top of the device and run to the bottom
c
      y = hy(1) * Ldi
      do 2200 j=1,ynodes
c
c       assign the voltages down in y for the different regions
c
c       put the correct potential guess into the 1D electrostatics
c
        if ( chmod .eq. 'MODFET' ) then
          if ( j .lt. joxide ) then
c				this point is in the oxide
            v(i,j) = ( phis + (joxide-j)*hy(1)*Ldi*(ksemi/kox)*
     +               (1.6022e-19*dopsub*xd)/(ksemi*8.854e-14) )/kT 
     +               + Vob
          else if ( (y - tox) .le. xd ) then
c				this is in the depletion region
            v(i,j) = -(1.6022e-19*dopsub/(ksemi*8.854e-14))*( xd*
     +               (y-tox) - (xd*xd + (y-tox)*(y-tox))/2.0 )/kT 
     +               + Vob
          else
c				this is the neutral region
            v(i,j) = Vob
          endif
        else
c
c         we have a buried channel
c
          if ( j .lt. joxide ) then
c				this point is in the oxide
            v(i,j) = ( phis + (joxide-j)*hy(1)*Ldi*(ksemi/kox)*
     +               ( (1.6022e-19*dopsub*xd)/(ksemi*8.854e-14) -
     +               1.6022e-19*tepi*(dopsub+dopch)/(ksemi*8.854e-14)
     +               ) )/kT + Vob
          else if ( (y - tox) .le. tepi ) then
c				this is in the depleted channel
            v(i,j) = ( (1.6022e-19*dopsub/(2.0*ksemi*8.854e-14))*
     +               ( xd*xd - tepi*tepi - 2.0*(y-tox)*(xd-tepi) ) -
     +               (1.6022e-19*dopch/(2.0*ksemi*8.854e-14))*
     +               (tepi - (y-tox))*(tepi - (y-tox)) )/kT + Vob
          else if ( (y - tox) .le. xd ) then
c				this is in the p-type depletion region
            v(i,j) = (1.6022e-19*dopsub/(ksemi*8.854e-14))*( -xd *
     +               (y-tox) + ( xd*xd + (y-tox)*(y-tox) )/2.0 )/kT 
     +               + Vob
          else
c				this is in the neutral region
            v(i,j) = Vob
          endif
        endif
c
c       assign the electron concentration in a similar fashion,
c       noting that the left side controls injection
c
        if ( j .gt. joxide ) then
          if ( symtyp .eq. 'Mock' ) then
            n(i,j) = exp( v(i,j) + log(dopsrc/nio) - rho(i) )
          else
            n(i,j) = exp( v(i,j) - rho(i) )
          endif
        else
c				n is zero in the oxide
          n(i,j) = 0.0
        endif
c
c       increment the position in real units for the next pass
c
        y = y + hy(j)*Ldi
 2200 continue
c
c     Now call oneDes to get the actual electrostatics for this first
c     slice, with no fixed electron concentration and gVs/kT for phin
c
      if ( verbose ) then
        write(6,'('' '')')
        write(6,'('' '',''Solving the electrostatics for the'',
     +        '' left boundary...'')')
      endif
      call oneDes( i, .false., rho(i), error, verbose, iter )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional electrostatic'',
     +       '' solver could not solve the left''/
     +            '' '',''      side boundary condition correctly.'',
     +       '' Shark is therefore Halting.'')')
        return
      endif
      if ( verbose ) then
        write(6,'('' '')')
        write(6,'('' '',''Solving for the exact n and V for the'',
     +       '' left boundary...'')')
      endif
      call full1D( i, error, verbose, iter )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional complete'',
     +       '' solver could not solve the left''/
     +            '' '',''      side boundary condition correctly.'',
     +       '' Shark is therefore Halting.'')')
        return
      endif
c
c     copy the results to srccon() for a boundary condition on V, if
c     we are doing a true center simulation
c
      if ( symtyp .eq. 'center' ) then
        do 2300 j=1,ynodes
          srccon(j,1) = v(i,j)
          srccon(j,2) = n(i,j)
          srccon(j,3) = 1.0
 2300   continue
      endif
c
c     Now guess the right side voltage for prepping into oneDes()
c
c     Get the constants for this slice in x
c
      i = xnodes
      Cox = kox * 8.854e-14 / tox
      Vo = 1.6022e-19 * ksemi * 8.854e-14 * dopsub / (Cox*Cox)
      V1 = (1.6022e-19*tepi/Cox)*( dopsub + dopch )
      phis = ( vgate(i)-Vob )*kT + 1.6022e-19*Qss(i)/Cox + Vo + V1 -
     +       sqrt( Vo*Vo + 2.0*Vo*( ( vgate(i)-Vob )*kT + 1.6022e-19*
     +       Qss(i)/Cox+V1 ) + V1*1.6022e-19*tepi*dopsub/Cox )
      xd = sqrt( tepi*tepi*(1.0 + dopch/dopsub)
     +           + 2.0*ksemi*8.854e-14*phis/(1.6022e-19 * dopsub) )
c
c     start at the top of the device and run to the bottom
c
      y = hy(1) * Ldi
      do 2400 j=1,ynodes
c
c       assign the voltages down in y for the different regions
c
c       put the correct potential guess into the 1D electrostatics
c
        if ( chmod .eq. 'MODFET' ) then
          if ( j .lt. joxide ) then
c				this point is in the oxide
            v(i,j) = ( phis + (joxide-j)*hy(1)*Ldi*(ksemi/kox)*
     +               (1.6022e-19*dopsub*xd)/(ksemi*8.854e-14) )/kT 
     +               + Vob
          else if ( (y - tox) .le. xd ) then
c				this is in the depletion region
            v(i,j) = -(1.6022e-19*dopsub/(ksemi*8.854e-14))*( xd*
     +               (y-tox) - (xd*xd + (y-tox)*(y-tox))/2.0 )/kT 
     +               + Vob
          else
c				this is the neutral region
            v(i,j) = Vob
          endif
        else
c
c         we have a buried channel
c
          if ( j .lt. joxide ) then
c				this point is in the oxide
            v(i,j) = ( phis + (joxide-j)*hy(1)*Ldi*(ksemi/kox)*
     +               ( (1.6022e-19*dopsub*xd)/(ksemi*8.854e-14) -
     +               1.6022e-19*tepi*(dopsub+dopch)/(ksemi*8.854e-14)
     +               ) )/kT + Vob
          else if ( (y - tox) .le. tepi ) then
c				this is in the depleted channel
            v(i,j) = ( (1.6022e-19*dopsub/(2.0*ksemi*8.854e-14))*
     +               ( xd*xd - tepi*tepi - 2.0*(y-tox)*(xd-tepi) ) -
     +               (1.6022e-19*dopch/(2.0*ksemi*8.854e-14))*
     +               (tepi - (y-tox))*(tepi - (y-tox)) )/kT + Vob
          else if ( (y - tox) .le. xd ) then
c				this is in the p-type depletion region
            v(i,j) = (1.6022e-19*dopsub/(ksemi*8.854e-14))*( -xd *
     +               (y-tox) + ( xd*xd + (y-tox)*(y-tox) )/2.0 )/kT 
     +               + Vob
          else
c				this is in the neutral region
            v(i,j) = Vob
          endif
        endif
c
c       assign the electron concentration in a similar fashion,
c       noting that the right side controls injection
c
        if ( j .gt. joxide ) then
          if ( symtyp .eq. 'Mock' ) then
            n(i,j) = exp( v(i,j) + log(dopdrn/nio) - rho(i) )
          else
            n(i,j) = exp( v(i,j) - rho(i) )
          endif
        else
c				n is zero in the oxide
          n(i,j) = 0.0
        endif
c
c       increment the position in real units for the next pass
c
        y = y + hy(j)*Ldi
 2400 continue
c
c     Now call oneDes to get the actual electrostatics for this last
c     slice, with no fixed electron concentration and gVs/kT for phin
c
      if ( verbose ) then
        write(6,'('' '')')
        write(6,'('' '',''Solving the electrostatics for the'',
     +        '' right boundary...'')')
      endif
      call oneDes( i, .false., rho(i), error, verbose, iter )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional electrostatic'',
     +       '' solver could not solve the right''/
     +            '' '',''      side boundary condition correctly.'',
     +       '' Shark is therefore Halting.'')')
        return
      endif
      if ( verbose ) then
        write(6,'('' '')')
        write(6,'('' '',''Solving for the exact n and V for the'',
     +       '' right boundary...'')')
      endif
      call full1D( i, error, verbose, iter )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional complete'',
     +       '' solver could not solve the right''/
     +            '' '',''      side boundary condition correctly.'',
     +       '' Shark is therefore Halting.'')')
        return
      endif
c
c     copy the results to drncon() for a boundary condition on V, if
c     we are doing a true center simulation
c
      if ( symtyp .eq. 'center' ) then
        do 2500 j=1,ynodes
          drncon(j,1) = v(i,j)
          drncon(j,2) = n(i,j)
          drncon(j,3) = 1.0
 2500   continue
      endif
c
c     If we are only to be doing the boundary conditions, then skip the
c     rest of this, and return to the caller
c
      if ( BConly ) then
        do 2600 j=1,ynodes
          v(1,j) = Vl(j)
          v(xnodes,j) = Vr(j)
          n(1,j) = nl(j)
          n(xnodes,j) = nr(j)
 2600   continue
        return
      endif
c
c     Now scan from the left side to the right side of the device
c     passing each line to oneDes() and using the previous line
c     as the initial guess.
c
      if ( verbose ) write(6,'('' '')')
c
c     get the results from the BC line to be used as the guess for the
c     first line
c
      do 2700 j=1,ynodes
        v(2,j) = v(1,j)
        n(2,j) = n(1,j)
 2700 continue
      do 2800 i=2,(xnodes-1)
c
c       call oneDes() to get the answer with rho(i) being phin(x)
c
        if ( verbose ) write(6,'('' '',''Solving the electrostatic '',
     +            ''problem for slice i='',i3,''.'')') i
c				rho(i) is the quasi-fermi level
        call oneDes( i, .false., rho(i), error, verbose, iter )
        if ( error ) then
          write(6,'('' '',''***** Error! *****''/
     +              '' '',''      The one dimensional electrostatic'',
     +           '' solver could not solve the slice i='',i3/
     +              '' '',''      problem correctly.'',
     +           '' Shark is therefore Halting.'')') i
          return
        endif

        if ( verbose ) write(6,'('' '',''Solving the exact '',
     +           ''n and V problem for slice i='',i3,''.'')') i
        call full1D( i, error, verbose, iter )
        if ( error ) then
          write(6,'('' '',''***** Error! *****''/
     +              '' '',''      The one dimensional complete'',
     +           '' solver could not solve the slice i='',i3/
     +              '' '',''      problem correctly.'',
     +           '' Shark is therefore Halting.'')') i
          return
        endif
c
c       copy these results to the next line to be used as a guess
c
        if ( (i+1) .lt. xnodes ) then
          do 2750 j=1,ynodes
            v((i+1),j) = v(i,j)
            n((i+1),j) = n(i,j)
 2750     continue
        endif
 2800 continue
c
c     now by-pass the charge neutral initial guess routine
c
      goto 3400
c
c
c     Now do the charge neutral initial guess because it was
c     requested
c
c
 2900 if ( verbose ) then
        write(6,'('' '',''    Computing the initial guess based on'',
     +     '' simple charge-neutrality.'')')
      endif
c
c     scan through each point and assume charge neutrality
c
      do 3000 j=1,ynodes
        do 3000 i=1,xnodes
c
c         get the charge-neutral voltage in normalized units,
c         and compute the electron concentration based on this
c
          if ( dop(i,j) .ne. 0.0 ) then
            v(i,j) = sign( log( abs(dop(i,j))/2.0 + sqrt(
     +               (dop(i,j)*dop(i,j) + 4.0)/4.0 )), dop(i,j) )
            if ( symtyp .eq. 'Mock' ) then
              n(i,j) = exp( v(i,j) + log(dopsrc/nio) )
            else
              n(i,j) = exp( v(i,j) )
            endif
          else
c				we are in the oxide, so zero the
c				voltage and make n = 0
            v(i,j) = 0.0
            n(i,j) = 0.0
          endif
 3000 continue

c
c     When the gate is zero, we should use the source/drain doping to
c     figure out what the potential of the oxide should be.
c
c     rho(i) is the potential we will seed the oxide with.
c     We will do this by scanning the initial voltage until we get
c     non-zero voltages where the source and drain should be.
c
      j = 1
 3100 if ( (v(1,j) .eq. 0.0) .or. (v(xnodes,j) .eq. 0.0) ) then
        j = j + 1
        goto 3100
      endif
      do 3200 i=1,xnodes
 3200   rho(i) = v(i,j)
c				now seed the oxide, only need the
c				lines up to where we stopped above
      do 3300 i=1,j
        do 3300 i1=1,xnodes
          if ( (i .le. joxide) .and. (v(i1,i) .eq. 0.0) ) then
            v(i1,i) = rho(i1)
          endif
 3300 continue

c
c     We need to limit n to 10x the limit on rho=0, but only in the
c     semiconductor
c
 3400 do 3500 j=1,ynodes
        do 3500 i=1,xnodes
          if ( j .le. joxide ) then
            n(i,j) = 0.0
          else
            n(i,j) = max( n(i,j), (1.0/(10.0*nio)) )
          endif
 3500 continue
c
c     finally, call depmap to show the initial guess to the user
c
      if ( verbose ) then
        call depmap( 1,1 )
      endif
c
c     all done... so return to the caller
c
      return
      end



      subroutine paosah( rho, aVs, aVd, aVg1, error )
c
c     This routine evaluates the Pao-Sah equations and sets the quasi-fermi
c     level versus position for the FET initial guess.
c
c     This was first written by J.A.Cooper Jr., converted to be compatible
c     with Shark by R.E.Beaty
c
c     first the parameters for the sizes of the arrays
c
      integer Nx, Ny
      parameter( Nx=XMAX, Ny=YMAX )
c
c     now for the arguments to this routine
c
      logical error
      PRECISION rho(Nx), aVs, aVd, aVg1
c
c     now for the functions I'll be calling in this subroutine
c
      PRECISION potent, charge
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
      common /pao/ es, Cox, ld, uf
      PRECISION es, Cox, ld, uf
c
c     local variables:
c       nxi = number of 'Xi' increments between 0 and Uds
c
      integer i, j, nxi
      PRECISION xi(0:(4*Nx)), us(0:(4*Nx)), y(0:(4*Nx)),
     +          n(0:(4*Nx)), intqn(0:(4*Nx)), q, w, usb,
     +          uds, ug(0:(4*Nx)), l, xPos
c
c     Set up some constants we will need
c
      nxi = 4*Nx
      q = 1.6022e-19
      es = ksemi * 8.854e-14
      Cox = kox * 8.854e-14 / tox
      do 10 i=0,nxi
   10   ug(i) = (aVg1 + q*Qss(int(xnodes*(1.0*i/(nxi+1.0)))+1)/Cox)/kT
      usb = aVs/kT
      uds = (aVd - aVs)/kT
c
c     Check to see if we have any drain to source voltage
c
      if ( uds .eq. 0.0 ) then
        do 50 i=1,xnodes
   50     rho(i) = usb
        return
      endif
c
c     ...uds isn't zero, so go on normally
c
      l = length * 1.0e-4
      ld = sqrt( es*kT/(2.0*q*nio) )
      uf = log(dopsub/nio)
      w = 1.0e-4
c
c     Begin calculation of variables down the channel
c
      do 100 i=0,nxi
        xi(i) = usb + uds*(1.0*i/(1.0*nxi))
        us(i) = potent( xi(i), ug(i) )
        n(i) = charge( us(i), xi(i) )
        intqn(i) = n(i)*uds/(1.0*nxi)
        if ( i .gt. 0 ) intqn(i) = intqn(i) + intqn(i-1)
  100 continue
      do 200 i=0,nxi
        y(i) = intqn(i)/intqn(nxi)*l
        intqn(i) = -intqn(i)*kT**2*w*mobil*es/(2*ld**2*nio*l)
  200 continue
c
c     On the left end of the device, there's no interpolation possible
c     as it's a fixed point.
c
      rho(1) = xi(0)
c
c     Now I need to interpolate on xi(i) and y(i)/l to get
c     rho(2..(xnodes-2)).
c
      xPos = hx(1) * Ldi
      do 500 i=2,(xnodes-1)
c
c       Since it's possible that the value will lie in front
c       of the first data point, we need to see where it is...
c
        if ( (xPos/length) .le. (y(0)/l) ) then
c
c         we have to assume that it's pegged at the source
c
          rho(i) = rho(1)
        else
c
c         find the appropriate bracketing values, and interpolate
c
          do 300 j=0,(nxi-1)
            if ( ((y(j)/l) .lt. (xPos/length)) .and. 
     +           ((y(j+1)/l) .ge. (xPos/length)) ) then
              rho(i) = xi(j) + ((xi(j+1) - xi(j))/(y(j+1) - y(j)))*
     +                          (xPos*l/length - y(j))
              goto 400
            endif
  300     continue
        endif
  400   xPos = xPos + hx(i) * Ldi
  500 continue
c
c     On the right end of the device, there's no interpolation possible
c     as it's a fixed point.
c
      rho(xnodes) = xi(nxi)
      error = .false.
      return
      end



      PRECISION function potent( xi, ug )
c
c     This subroutine calculates us given xi, ug and uf
c
c
c     Now for the arguments to this routine
c
      PRECISION xi, ug
c
c     Next, the common block needed for this routine
c
      common /pao/ es, Cox, ld, uf
      PRECISION es, Cox, ld, uf
c
c     local variables
c
      PRECISION us, dus, f

      us = 0.0
      dus = 100.0
      us = us - dus
  100 us = us + dus
      f = sqrt( max( exp(uf)*( exp(-us) + us - 1 ) + exp(-uf)*
     +          ( exp(us-xi) - us - exp(-xi) ), 0.0 ) )
      if ( (us+es/(Cox*ld)*f) .le. ug ) go to 100
      if ( dus .lt. 0.01 ) then
        potent = us
        return
      endif
      us = us - dus
      dus = dus/10.0
      go to 100
      end



      PRECISION function charge( us, xi )
c
c     This subroutine calculates qn, given us, uf and xi
c
c     first the parameters for the sizes of the arrays
c
      integer Nx, Ny
      parameter( Nx=XMAX, Ny=YMAX )
c
c     now for the arguments to this routine
c
      PRECISION us, xi
c
c     next the common blocks for interface to the rest of the routines
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
      common /pao/ es, Cox, ld, uf
      PRECISION es, Cox, ld, uf
c
c     local variables:
c       nui = number of Ui increments between Uf and Us
c
      integer i, nui
      PRECISION sum, du, u, f

      nui = Nx
      du = (us-uf)/nui
      sum = 0.0
      do 100 i=1,nui
        u = uf + du*(i - 0.5)
        f = sqrt( exp(uf)*( exp(-u) + u - 1.0 ) + exp(-uf)*
     +            ( exp(u - xi) - u - exp(-xi)) )
  100   sum = sum + exp(u - uf - xi) * du / f
      charge = nio * ld * sum
      return
      end



      PRECISION function Qn2rho( Qn, i, tol, iter, error )
c
c     This function returns the value of the quasi-fermi level, 
c     for the given slice i, that will produce the inversion layer 
c     charge Qn. This will be carried out in the binary search 
c     pattern and will be continued until the error in the 
c     estimate of Qn is within the 'tol' percentage of Qn (i.e. 
c     abs(Qest-Qn) .lt. tol*Qn )
c
c
c     first the parameters for the sizes of the arrays
c
      integer Nx, Ny
      parameter( Nx=XMAX, Ny=YMAX )
c
c     Now for some parameters that are constants at the machine's
c     precision
c
      PRECISION oneK
      parameter( oneK=1000.0 )
c
c     now for the arguments to this routine
c
      integer i, iter
      logical error
      PRECISION Qn, tol
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
c     Now I need some local variables
c
      integer j, iter1D
      PRECISION Qntemp, dely, rho, rhoh, rhol, vimax, Qnmax, tau
c
c     Assume the result is zero - then change it from this
c
      Qn2rho = 0.0
c
c     Start with a guess on the limits on phin, and a try at the right
c     one. Assume the lowest phin is zero (equilibrium) and the largest
c     is the surface potential plus a safety margin.
c
      rhol = 0.0
      call oneDes( i, .false., oneK, error, .false., iter1D )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional electrostatic'',
     +      '' solver could not be solved.''/
     +            '' '',''      Shark is therefore Halting.'')')
        return
      endif
c
c     Get the maximum value of the voltage in this line
c
      vimax = -1.0
      do 100 j=1,ynodes
  100   vimax = max( vimax, v(i,j) )
c
c     Set the high limit at this level
c
      rhoh = vimax
c
c     Reset the number of iterations required
c
      iter = 0
c
c     Make sure Qn falls into these limits, by checking rhol first.
c
  200 call oneDes( i, .false., rhol, error, .false., iter1D )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional electrostatic'',
     +      '' solver could not be solved.''/
     +            '' '',''      Shark is therefore Halting.'')')
        return
      endif
      Qntemp = 0.0
      do 300 j=1,ynodes
c				get the delta y for this point
        if ( j .eq. 1 ) then
          dely = hy(j)*Ldi
        elseif ( j .eq. ynodes ) then
          dely = hy(j-1)*Ldi
        else
          dely = 0.5*( hy(j) + hy(j-1) )*Ldi
        endif
  300   Qntemp = Qntemp + n(i,j)*nio*dely
      Qnmax = Qntemp
c
c     Check it, and move limits if need be
c
      if ( Qn .gt. Qntemp ) then
        rhoh = rhol
        rhol = rhol - 2.0
        iter = iter + 1
        goto 200
      endif
c
c     Now check the upper limit
c
  400 call oneDes( i, .false., rhoh, error, .false., iter1D )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional electrostatic'',
     +      '' solver could not be solved.''/
     +            '' '',''      Shark is therefore Halting.'')')
        return
      endif
      Qntemp = 0.0
      do 500 j=1,ynodes
c				get the delta y for this point
        if ( j .eq. 1 ) then
          dely = hy(j)*Ldi
        elseif ( j .eq. ynodes ) then
          dely = hy(j-1)*Ldi
        else
          dely = 0.5*( hy(j) + hy(j-1) )*Ldi
        endif
  500   Qntemp = Qntemp + n(i,j)*nio*dely
c
c     Get the decay constant for the exponential
c
      tau = -rhoh/log(Qntemp/Qnmax)
c
c     Check it, and move limits if need be
c
      if ( Qn .lt. Qntemp ) then
        rhol = rhoh
        rhoh = rhoh + 2.0
        iter = iter + 1
        goto 400
      endif
c
c     Now take the first guess at phin
c
      rho = -tau*log( Qn/Qnmax )
c
c     Get a n(i,y) based on this guess at phin
c
  600 call oneDes( i, .false., rho, error, .false., iter1D )
      if ( error ) then
        write(6,'('' '',''***** Error! *****''/
     +            '' '',''      The one dimensional electrostatic'',
     +      '' solver could not be solved.''/
     +            '' '',''      Shark is therefore Halting.'')')
        return
      endif
c
c     Integrate n(x,j) to get Qn(x)
c
      Qntemp = 0.0
      do 700 j=1,ynodes
c				get the delta y for this point
        if ( j .eq. 1 ) then
          dely = hy(j)*Ldi
        elseif ( j .eq. ynodes ) then
          dely = hy(j-1)*Ldi
        else
          dely = 0.5*( hy(j) + hy(j-1) )*Ldi
        endif
  700   Qntemp = Qntemp + n(i,j)*nio*dely
c
c     Now see if we need to do all this again...
c
      if ( abs(Qntemp - Qn) .gt. (tol*Qn) ) then
        if ( Qntemp .gt. Qn ) then
          rhol = rho
        else
          rhoh = rho
        endif
        rho = 0.5*(rhoh + rhol)
        iter = iter + 1
        goto 600
      endif
c
c     We must have converged, so return the value
c
      Qn2rho = rho
      return
      end
