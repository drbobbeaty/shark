c
c     Shark - General 2-D analysis package for FETs
c
c     This is a two dimensional extension of the time-of-flight
c     simulator first written by J.A. Cooper, Jr. It was found to 
c     be necessary to model the device more accurately so this
c     two dimensional simulator was written by R.E. Beaty.
c
c     The main variables are as follows:
c       n(i,j) - Normalized carrier (electron) concentration. This
c                is normalized by nio, the intrinsic carrier
c                concentration.
c       dop(i,j) - Normalized dopant concentration. This is assumed
c                  to be fully ionized, and given by (Nd-Na). This
c                  is normalized by nio.
c       v(i,j) - Normalized potential in units of kT/q.
c       hx(i),hy(j) - Normalized grid spacings. This is normalized
c                     by Ldi, the intrinsic Debye length.
c       srccon(j,mode), drncon(j,mode) - left and right side boundary
c                     conditions which depends on mode:
c                     srccon(j,1) = V's for the boundary, if fixed
c                     srccon(j,2) = n's for the boundary, if fixed
c                     srccon(j,3) = mode: =1.0 then fixed n, V
c                                         =0.0 then LOS
c
c     All non-normalized distances are in cm.
c
c     The indices are ordered as follows: j, the second index, is
c     the position down towards the substrate starting at the
c     metal-insulator interface, and i, the first index, is the
c     position down the channel from the source to the drain.
c
c     This has been written by R.E. Beaty as part of a requirement 
c     for a PhD in Electrical Engineering at Purdue University,
c     please check before copying or distributing the source.
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
c     These parameters are for the adaptive time stepping
c
      integer bumpdt
      PRECISION dtmax
      parameter( bumpdt=10, dtmax=0.02e-12 )
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
      common /JF/ J( (2*Nx*Ny), 10 ), delu( 2*Nx*Ny ), F( 2*Nx*Ny )
      PRECISION J, delu, F
      common /Pshape/ Qm, ftw, fwhm, pulse
      PRECISION Qm, ftw, fwhm
      logical pulse
c
c     now for the common blocks for interface to the parser
c
      common /parsd1/ kflag, numr, numc, vreal(70)
      integer kflag, numr, numc
      real vreal
      common /parsd2/ ncard, nreal(70), nchar(30),
     +                vchar(30), ititle(10), ibuf(400)
      character*8 ncard, nreal, nchar, vchar, ititle
      character*1 ibuf
c
c     now I need to declare the functions I need
c
      integer iargc
      real etime, dtime
c
c     now the local veriables to the main routine
c
      integer i, k, inperr, iter, gooddt, days, hours, mins, secs,
     +        deltai, ibuff
      real time, tarray(2)
      PRECISION tnow, tcnt, xref, deltax, Qn(Nx),
     +          dely, delx, Qtot
      logical error, fileok
      character*32 lstfil
#     ifdef Cyber
      logical restart
#     endif
c
c     Reset the timer routines so we know where we start
c
      time = dtime(tarray)
c
c     Reset the number of elapsed hours of run time
c
      hours = 0
c
c     Reset the location of the starting location of the channel for
c     the drifting channel simulations
c
      xref = 0.0

#     ifndef Cyber
c
c     Parse the command line
c
      if ( iargc() .ge. 1 ) then
        call getarg( 1, lstfil )
c
c       Do we have a restart?
c
        if ( lstfil .eq. 'restart' ) then
          lstfil = 'shark.core'
          if ( iargc() .ge. 2 ) then
            call getarg( 2, lstfil )
          endif
          write(6,'('' '',''... Restarting from '',a,''...'')') lstfil
c
c         this is the entry point for a restart from a core dump
c
          call load( lstfil, tnow, xref )
          hours = 0
          iter = 0
          gooddt = 0
          tcnt = 0.0
c
c         now let's see if we have an output file specified, and if
c         so let's get it and open it for appending.
c
          if ( iargc() .eq. 3 ) then
            call getarg( 3, lstfil )
#           ifdef Sun3
              open( 6, file=lstfil, status="UNKNOWN", fileopt="eof" )
#           else
              open( 6, file=lstfil, status="UNKNOWN" )
#           endif
          endif
#         ifndef Cyber
c
c         open the output file, at the end if possible
c
          open( 9, file="shark.cplot", status="UNKNOWN" )
#         endif
c
c         now jump into the thick of it...
c
          goto 400
        else
c
c         No restart, so open the input deck
c
          inquire( file=lstfil, exist=fileok )
          if ( fileok ) then
            open( 8, file=lstfil, status="OLD" )
            rewind 8
          else
            write(6,'('' '',''I am sorry, you have to specify an '',
     +            ''existing file for the source deck for Shark.''/
     +                '' '',''Like:''/
     +                '' '',''  $ shark input_deck [output_file]''/
     +                '' '',''or:''/
     +                '' '',''  $ shark restart [core_file] [output_file]''/
     +                '' '',''Where input_deck, or core_file, is a '',
     +            ''valid, existing, file.'')')
            stop
          endif
c
c         now let's see if we have an output file specified, and if
c         so let's get it and open it
c
          if ( iargc() .eq. 2 ) then
            call getarg( 2, lstfil )
            open( 6, file=lstfil, status="UNKNOWN" )
            rewind(6)
          endif
        endif
      else
c
c       No file names - Error!
c
        write(6,'('' '',''I am sorry, you have to specify an input '',
     +                  ''file for the source deck for Shark.''/
     +            '' '',''Like:''/
     +            '' '',''  $ shark input_deck [output_file]''/
     +            '' '',''or:''/
     +            '' '',''  $ shark restart [core_file] [output_file]''/
     +            '' '',''Where input_deck, or core_file, is a '',
     +                  ''valid, existing, file.'')')
        stop
      endif
#     endif
c
c     now give this program a reasonable header, something to remember
c     it by
c
    2 format(' ',a115)
      write(6,2) '**************************************************'//
     + '************************************************************'//
     + '*****'
      write(6,2) '*                                                 '//
     + '                                                            '//
     + '    *'
      write(6,2) '*         SSSSSSSSS    h             '//
     +   '                   k                        2'//
     +   '222222        000000000         *'
      write(6,2) '*        S             h           aa'//
     +   'aaaa               k                       2 '//
     +   '      2      0         0        *'
      write(6,2) '*        S             h             '//
     +   '    a   r rrrrr    k   kkk                   '//
     +   '      2      0         0        *'
      write(6,2) '*         SSSSSSSSS    h hhhhh     aa'//
     +   'aaa a   rr     r   k kkk          v     v   2'//
     +   '222222       0         0        *'
      write(6,2) '*                  S   hh     h   a  '//
     +   '   aa   r          kkk             v   v   2 '//
     +   '             0         0        *'
      write(6,2) '*                  S   h      h   a  '//
     +   '    a   r          k  kkk           v v    2 '//
     +   '         ..  0         0        *'
      write(6,2) '*         SSSSSSSSS    h      h    aa'//
     +   'aaaaa   r          k    kkk          v     22'//
     +   '2222222  ..   000000000         *'
      write(6,2) '*                                                 '//
     + '                                                            '//
     + '    *'
      write(6,2) '*                                                 '//
     + '                                                            '//
     + '    *'
      write(6,2) '*            T W O   D I M E N S I O N A L   T R A'//
     + ' N S I E N T   G a A s   F E T   S I M U L A T O R          '//
     + '    *'
      write(6,2) '*                                                 '//
     + '                                                            '//
     + '    *'
      write(6,2) '*                  Written by Robert E. Beaty as a'//
     + ' part of a Ph.D. in Electrical Engineering at               '//
     + '    *'
      write(6,2) '*                           Purdue University, Wes'//
     + 't Lafayette, Indiana 47906  June 1988                       '//
     + '    *'
      write(6,2) '*                                                 '//
     + '                                                            '//
     + '    *'
      write(6,2) '**************************************************'//
     + '************************************************************'//
     + '*****'
      write(6,2) ' '
      write(6,2) ' '
      write(6,2) ' '
#     ifndef Cyber
        call flush(6)
#     endif

#     ifndef Cyber
c
c     open the output file
c
      open( 9, file="shark.cplot", status="UNKNOWN" )
#     endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     S E T - U P   D E F A U L T S   F O R   S I M U L A T O R
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     The default channel is GaAs, so set the parameters for GaAs
c     ( these arefrom App. H os Sze second edition )
c
      chtyp = 'GaAs'
      nio = 1.79e6
      ksemi = 12.9
c
c     The default insulator is SiN, so set the parameters for SiN
c     ( also from Sze II )
c
      instyp = 'SiN'
      kox = 7.5
      uins = 1.0e-3
c				this is the surface charge density at
c				the insulator-semiconductor interface
      do 150 i=1,Nx
  150   Qss(i) = 0.0
c				this is the value of kT/q at room temp
      kT = 0.0259
c				these are the number of grid points
c				for the device
      xnodes = Nx
      ynodes = Ny
c				this is the fraction of xnodes that will
c				be placed around each N+ region
      xfrac = 0.4
c				this is yes/no for automatic mesh selection
      auto = 'yes'
c				this is the number of nodes to be placed
c				in the oxide
      joxide = 6
c				this is the difference in voltage, in kT/q,
c				for the y nodes to be placed in the y
c				direction, and for refinement purposes
      dVa = 10.0
c
c     The default channel and source/drain dopings are step junctions. So
c     set the thickness of the epi-layer and the thickness of the insulator
c     (oxide). The length of the source and drain are also set here for
c     defaults.
c
c				this is the channel region model itself
      chmod = 'MODFET'
c				this is the doping profile type for the 
c				channel and the source and drain
      doptyp = 'step'
c				this is the doping in the channel region
      dopch = 1.0e17
c				this is the doping in the source and drain
      dopsrc = 1.0e19
      dopdrn = 1.0e19
c				this is the doping in the substrate
      dopsub = 1.0e15
c				this is the thickness of the channel region
      tepi = 0.15e-4
c				this is the thickness of the insulator
      tox = 1000e-8
c				this is the thickness of the substrate
      tsub = 5.0e-4
c				this is the length of the channel in cm
      length = 10e-4
c				this is the length of the source in cm
      lsorce = 0.0
c				this is the thickness of the source in cm
      tsorce = 0.35e-4
c				this is the type of corner on the source
      srccor = 'square'
c				this is the gate to source overlap
      gsol = 0.35e-4
c				this is the length of the drain in cm
      ldrain = 0.0
c				this is the thickness of the drain in cm
      tdrain = 0.35e-4
c				this is the type of corner on the drain
      drncor = 'square'
c				this is the gate to drain overlap
      gdol = 0.35e-4
c				this is the starting time for the simulation
      tstart = 0.0
c				this is the stopping time of the simulation
      tstop = 50.0e-12
c				this is the minimum allowable time step
      stpmin = 2.0e-17
c				this is the default value of the time step
c				during the main section of the program
      dt = 0.02e-12
c				this is the value of the ramping time step
c				in the biasing phase of the program
      delt = 1.0e32
c				this is the time between outputs
      tout = 10.0e-12
c				this is the source bias voltage
      Vs = 0.0
c				this is the drain bias voltage
      Vd = 20.0
c				this is the average electric field in the 
c				channel
      eavg = 7000.0
c				this is the bias potential on the source-end
c				gate contact
      vg1 = 1.5
c				this is the maximum allowable voltage step
c				in bias()
      delV = 0.5
c				this is the maximum allowable Eavg step in
c				bias()
      delE = 100.0
c				this is the line width of the flash grating
      line = 1.0e-4
c				this is the space width of the flash grating
      space = 1.0e-4
c				this is the magnitude of the flash pulse
      npulse = 1.0e12
c				this is the defocus parameter for smearing
c				of pulses
      defoc = 0.05
c				this is the initial condition of the channel
      icchan = 'equil'
c				this is the stopping criteria for newton()
      vmax = 1.0e-6
c				this is the maximum number of iterations for
c				the solution of newton's method
      maxit = 200
c				this is the type of solution that is used to
c				solve the system of equations on each
c				Newton iteration
      soltyp = 'linpk'
c				this is the number of iterations kept in the
c				Orthomin(k) algorithim
      keep = 6
c				this is the flag that will tell us if we are
c				to be really wordy about the solution
      debug = 'no'
c				this will tell us the initial guess to be used
c				when starting the iterating
      guess = 'MOS'
c				this is the maximum value of the drifting
c				pulse in 1/cm^2
      Qm = 8.0e11
c				this is the width of the flat top on the
c				pulse (um)
      ftw = 0.0
c				this is the full width half maximum of the
c				pulse NOT counting the flat top section
      fwhm = 1.0


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read in the deck and see if there are any changes
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#     ifdef Cyber
      call getdk( inperr, restart )
      if ( restart ) then
        write(6,'('' '',''... Restarting from dump file...'')')
c
c       This is the main entry point for a restart from a core dump
c       on the Cyber
c
        call load( lstfil, tnow, xref )
        hours = 0
        iter = 0
        gooddt = 0
        tcnt = 0.0
        goto 400
      endif
#     else
      call getdk( inperr )
#     endif
      if ( inperr .ne. 0 ) goto 9999
c
c     compute some necessary constants that depend on the simulation
c
      Ldi = sqrt( ( ksemi*8.854e-14 * kT )/(1.6022e-19*nio) )
c
c     normalize the mobility here because it makes sense to
c
      uins = uins * kT
c
c     check to see that the parameters read in are valid ones
c
      if ( auto .ne. 'yes' ) then
        if ( xfrac .ge. 0.5 ) then
          write(6,'('' '',''*** Error! ***''/
     +              '' '',''    You have specified that more than '',
     +                    ''half of the nodes in the x-direction be '',
     +                    ''placed in''/
     +              '' '',''    each of the source and drain '',
     +                    ''regions. This is not a possible '',
     +                    ''configuration. Please place''/
     +              '' '',''    more nodes in the channel region. '',
     +                    ''Shark is halting.'')')
          goto 9999
        endif

        if ( (xnodes*(1.0-2.0*xfrac)*1.0e-4/length) .lt. 10.0 ) then
          write(6,'('' '',''*** Warning! ***''/
     +              '' '',''    You have specified less than 10 '',
     +                    ''points per micron in the channel. This '',
     +                    ''is dangerous''/
     +              '' '',''    as the oscillations may not have '',
     +                    ''enough points in one period to define ''
     +                    ''them well''/
     +              '' '',''    enough. You should either increase '',
     +                    ''xnodes or decrease the length.'')')
        endif
c
c       check the rest of the mesh selection parameters
c
        if ( xnodes .gt. Nx ) then
          write(6,'('' '',''*** Error! ***''/
     +              '' '',''    You want more nodes in the '',
     +                    ''x-direction than the simulator is '',
     +                    ''capable of. Please''/
     +              '' '',''   select a value less than '',i3,
     +                    '' for the value of xnodes.'')') Nx
          goto 9999
        endif

        if ( ynodes .gt. Ny ) then
          write(6,'('' '',''*** Error! ***''/
     +              '' '',''    You want more nodes in the '',
     +                    ''y-direction than the simulator is '',
     +                    ''capable of. Please''/
     +              '' '',''   select a value less than '',i3,
     +                    '' for the value of ynodes.'')') Ny
          goto 9999
        endif
      endif
c
c     the rest need to be check even with automatic mesh selection
c
      if ( (lsorce+ldrain) .gt. (length*0.5) ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the drain '',
     +                  ''and source are greater than half the''/
     +            '' '',''    length of the channel! This is '',
     +                  ''dangerous.'')')
      endif

      if ( tepi .gt. tsub ) then
        write(6,'('' '',''*** Error! ***''/
     +            '' '',''    You have specified that the thickness '',
     +                  ''of the channel is greater than the '',
     +                  ''thickness''/
     +            '' '',''    of the substrate. This is impossible. '',
     +                  ''Shark is Halting!'')')
        goto 9999
      endif
      if ( tepi .gt. (0.1*tsub) ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the thickness '',
     +                  ''of the channel is greater than one-tenth '',
     +                  ''the''/
     +            '' '',''    thickness of the substrate. This is '',
     +                  ''dangerous.'')')
      endif

      if ( (tsorce .gt. tsub) .and. (symtyp .eq. 'full') ) then
        write(6,'('' '',''*** Error! ***''/
     +            '' '',''    You have specified that the thickness '',
     +                  ''of the source is greater than the''/
     +            '' '',''    thickness of the substrate. This is '',
     +                  ''impossible. Shark is Halting!'')')
        goto 9999
      endif
      if ( tsorce .gt. (0.5*tsub) ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the thickness '',
     +                  ''of the source is greater than half''/
     +            '' '',''    the thickness of the substrate. This '',
     +                  ''is dangerous.'')')
      endif

      if ( (tdrain .gt. tsub) .and. (symtyp .eq. 'full') ) then
        write(6,'('' '',''*** Error! ***''/
     +            '' '',''    You have specified that the thickness '',
     +                  ''of the drain is greater than the''/
     +            '' '',''    thickness of the substrate. This is '',
     +                  ''impossible. Shark is Halting!'')')
        goto 9999
      endif
      if ( tdrain .gt. (0.5*tsub) ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the thickness '',
     +                  ''of the drain is greater than half''/
     +            '' '',''    the thickness of the substrate. This '',
     +                  ''is dangerous.'')')
      endif

c     depth = 7.04116*max( tsorce, tdrain )
c     thick = 2.0*depth + ((ynodes-6.0)/(ynodes-2.0))*
c    +                    (tox + depth)/3.0
c     if ( tsub .lt. thick ) then
c       write(6,'('' '',''*** Warning! ***''/
c    +            '' '',''    You have specified that the thickness '',
c    +                  ''of the substrate is less than that''/
c    +            '' '',''    suggested by the author of Shark for '',
c    +                  ''safe margins of certainty in the solution.''/
c    +            '' '',''    A safe number would have been '',f6.2,
c    +                  '' um. Shark is now using this value.'')') 
c    +          thick*1.0e4
c       tsub = thick
c     endif

      if ( (gdol .gt. (0.9*ldrain)) .and. (symtyp .eq. 'full') ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the overlap '',
     +                  ''of the gate with the drain is greater '',
     +                  ''than that''/
     +            '' '',''    suggested by the author of Shark '',
     +                  ''for safe margins of certainty in the '',
     +                  ''solution.''/
     +            '' '',''    A result of this is that the solution '',
     +                  ''may not converge due to an ohmic contact '',
     +                  ''being''/
     +            '' '',''    in the depletion region.'')')
      endif

      if ( (gsol .gt. (0.9*lsorce)) .and. (symtyp .eq. 'full') ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the overlap '',
     +                  ''of the gate with the source is greater '',
     +                  ''than that''/
     +            '' '',''    suggested by the author of Shark '',
     +                  ''for safe margins of certainty in the '',
     +                  ''solution.''/
     +            '' '',''    A result of this is that the solution '',
     +                  ''may not converge due to an ohmic contact '',
     +                  ''being''/
     +            '' '',''    in the depletion region.'')')
      endif

      if ( (lsorce .le. tsorce) .and. (symtyp .eq. 'full') ) then
        write(6,'('' '',''*** Error! ***''/
     +            '' '',''    You have specified the length of the '',
     +                  ''source to be less than the depth of the''/
     +            '' '',''    source, so the ion implanted junction '',
     +                  ''cannot be properly tapered. This would '',
     +                  ''cause''/
     +            '' '',''    very high fields to be created at '',
     +                  ''the corner of the junction. Shark is '',
     +                  ''halting.'')')
        goto 9999
      endif

      if ( (ldrain .le. tdrain) .and. (symtyp .eq. 'full') ) then
        write(6,'('' '',''*** Error! ***''/
     +            '' '',''    You have specified the length of the '',
     +                  ''drain to be less than the depth of the''/
     +            '' '',''    drain, so the ion implanted junction '',
     +                  ''cannot be properly tapered. This would '',
     +                  ''cause''/
     +            '' '',''    very high fields to be created at '',
     +                  ''the corner of the junction. Shark is '',
     +                  ''halting.'')')
        goto 9999
      endif

      if ( dopsrc .le. 0.0 ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the doping '',
     +                  ''in the source would be p-type. This is '',
     +                  ''not possible''/
     +            '' '',''    under the conditions of the '',
     +                  ''solution. Shark is taking the source '',
     +                  ''doping to be 10*nio.''/
     +            '' '',''    This is a bare minimum.'')')
        dopsrc = 10.0*nio
      endif
      if ( dopdrn .le. 0.0 ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified that the doping '',
     +                  ''in the drain would be p-type. This is '',
     +                  ''not possible''/
     +            '' '',''    under the conditions of the '',
     +                  ''solution. Shark is taking the drain '',
     +                  ''doping to be 10*nio.''/
     +            '' '',''    This is a bare minimum.'')')
        dopdrn = 10.0*nio
      endif
c
c     make sure that I have at least 2 nodes in the insulator
c
      if ( joxide .lt. 2 ) then
        write(6,'('' '',''*** Warning! ***''/
     +            '' '',''    You have specified less than 2 nodes '',
     +                  ''in the insulator. This is too low. Shark''/
     +            '' '',''    is setting it to 2.'')')
        joxide = 2
      endif
c
c     Make sure that if we are running with a pulse that we do a
c     center simulation
c
      if ( (pulse) .and. (.not. (symtyp .eq. 'center')) ) then
        write(6,'('' '',''If you are simulating a running pulse '',
     +                  ''then you must do a center''/
     +            '' '',''simultation without the source and drain. '',
     +                  ''Shark is now doing this for you.'')')
        symtyp = 'center'
      endif
c
c     If we are doing a pulse simulation, make sure that the source
c     and drain are positive enough that the center simulation will
c     think that phin->infinity so n->0 at the boundaries.
c
      if ( pulse ) then
        Vs = 1000.0
        Vd = 1000.0
      endif
c
c     output the simulation parameters for all to see
c
      write(6,'('' '')')
      write(6,'('' '')')
      write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                            REB 1988'',
     +     ''     run parameters ===='')')
      write(6,'('' '')')
      if ( auto .ne. 'yes' ) then
        write(6,'('' '',''   The length of the simulated device is '',
     +                  f6.2,'' (um).'')') length*1.0e4
      else
        write(6,'('' '',''   The length of the simulated channel is '',
     +                  f6.2,'' (um).'')') length*1.0e4
      endif
      write(6,'('' '',''   The thickness of the simulated device is '',
     +                f6.2,'' (um).'')') (tox+tsub)*1.0e4
      write(6,'('' '',''   The doping density of the substrate, '',
     +                ''given as (Na-Nd), is '',g8.2,
     +                '' (/cm**3).'')') dopsub
      write(6,'('' '')')
      if ( symtyp .eq. 'full' ) then
        write(6,'('' '',''   The structure of the source is:''/
     +            '' '',''        length = '',f4.2,'' (um).''/
     +            '' '',''        thickness = '',f5.3,'' (um).''/
     +            '' '',''        doping density = '',g8.2,
     +                  '' (/cm**3).'')')
     +        lsorce*1.0e4, tsorce*1.0e4, dopsrc
        write(6,'('' '')')
        write(6,'('' '',''   The structure of the drain is:''/
     +            '' '',''        length = '',f4.2,'' (um).''/
     +            '' '',''        thickness = '',f5.3,'' (um).''/
     +            '' '',''        doping density = '',g8.2,
     +                  '' (/cm**3).'')')
     +        ldrain*1.0e4, tdrain*1.0e4, dopdrn
        write(6,'('' '')')
      endif
      write(6,'('' '',''   The insulator is '',a6,'' with a '',
     +                ''thickness of '',f6.1,'' (Angst.).'')')
     +                instyp, tox*1.0e8
c     write(6,'('' '',''   ...and a surface charge density of '',
c    +                g10.4,'' (/cm**2).'')') Qss
      write(6,'('' '')')
      write(6,'('' '',''   The dielectric constants for the '',
     +                ''materials are:''/
     +          '' '',''        semiconductor = '',f5.2/
     +          '' '',''        insulator = '',f5.2)') ksemi, kox
      write(6,'('' '')')
      write(6,'('' '',''   The channel is based on the '',a6,
     +                '' channel.'')') chmod
      if ( chmod .eq. 'buried' ) then
        write(6,'('' '',''   with a thickness of '',f6.2,'' (um),''/
     +            '' '',''   and a doping density of '',g8.2,
     +                  '' (/cm**3).'')') tepi*1.0e4, dopch
      endif
      write(6,'('' '',''   The channel material type is '',a5,
     +                ''.'')') chtyp
      write(6,'('' '')')
      write(6,'('' '',''   The biasing on the terminals of the '',
     +                ''device is:''/
     +          '' '',''        source voltage = '',f6.3,'' (V).''/
     +          '' '',''        drain voltage = '',f6.3,'' (V).''/
     +          '' '',''        source end of the gate = '',f6.3,
     +                '' (V).'')') Vs, Vd, Vg1
      if ( symtyp .eq. 'full' ) then
        if ( auto .ne. 'yes' ) then
          write(6,69) Vg1+eavg*(length-(lsorce-gsol)-(ldrain-gdol))
        else
          write(6,69) Vg1+eavg*(length+gsol+gdol)
        endif
      else
        write(6,69) Vg1+eavg*length
      endif
   69 format(' ','        drain end of the gate = ',f6.3,' (V).')
      write(6,'('' '',''        field in the gate = '',f7.1,
     +                '' (V/cm).'')')
     +      eavg
      write(6,'('' '')')
      write(6,'('' '',''   The parameters for the velocity versus '',
     +                ''field equation are:''/
     +          '' '',''         low-field mobility = '',f7.1,
     +                '' (cm**2/V-s).''/
     +          '' '',''                          A = '',g11.4/
     +          '' '',''                          B = '',g11.4/
     +          '' '',''                          C = '',g11.4/
     +          '' '',''                          D = '',g11.4/
     +          '' '',''                          F = '',g11.4)')
     +         mobil, Av, Bv, Cv, Dv, Fv
      write(6,'('' '')')
      write(6,'('' '',''   The parameters for the diffusivity '',
     +                ''versus field equation are:''/
     +          '' '',''      low-field diffusivity = '',f7.2,
     +                '' (cm**2/sec).''/
     +          '' '',''           Do-to-peak value = '',f7.2,
     +                '' (cm**2/sec).''/
     +          '' '',''                          A = '',f6.3/
     +          '' '',''                          B = '',f7.4)')
     +         Do, dDo, Ad, Bd
      write(6,'('' '')')
      if ( auto .ne. 'yes' ) then
        write(6,'('' '',''   The mesh was selected by the '',
     +                  ''operator to contain:''/
     +            '' '',''        '',i3,'' nodes in the x direction''/
     +            '' '',''        '',f5.1,'' % of the x nodes near '',
     +                  ''the source and an equal number near the '',
     +                  ''drain.''/
     +            '' '',''        '',i3,'' nodes in the y '',
     +                  ''direction'')') xnodes, xfrac*100.0, ynodes
        write(6,'('' '')')
      endif
      if ( icchan .eq. 'DC' ) then
        write(6,'('' '',''   Starting the channel with the steady '',
     +                  ''state electron concentration in the '',
     +                  ''channel.'')')
      elseif ( icchan .eq. 'zero' ) then
        write(6,'('' '',''   Starting the channel with a near zero '',
     +                  ''electron concentration.'')')
      endif
      write(6,'('' '')')
      write(6,'('' '',''   The Newton iterations will halt if the '',
     +                ''largest change in the voltage is less than '',
     +                g8.2,'' kT/q,''/
     +          '' '',''   or '',i4,'' iterations have been '',
     +                ''completed.'')') vmax, maxit
      if ( soltyp .eq. 'minres' ) then
        write(6,'('' '',''   The systems of equations will be '',
     +                  ''solved by the Minimum Residual '',
     +                  ''algorithm.'')')
      elseif ( soltyp .eq. 'gauss' ) then
        write(6,'('' '',''   The systems of equations will be '',
     +                  ''solved by the Gauss-Seidel algorithm.'')')
      elseif ( soltyp .eq. 'ortho' ) then
        write(6,'('' '',''   The systems of equations will be '',
     +                  ''solved by the Orthomin('',i1,
     +                  '') algorithm.'')') keep
      elseif ( soltyp .eq. 'linpk' ) then
        write(6,'('' '',''   The systems of equations will be '',
     +                  ''solved by the small LINPACK algorithm.'')')
      endif
      write(6,'('' '')')
      write(6,'('' '',''   The starting time of the simulation is '',
     +                f5.1,'' (psec).''/
     +          '' '',''   The stopping time of the simulation is '',
     +                f5.1,'' (psec).''/
     +          '' '',''   Output will be generated every '',f5.1,
     +                '' (psec) of the simulation.'')')
     +        tstart*1.0e12, tstop*1.0e12, tout*1.0e12
      write(6,'('' '')')
      write(6,'('' '')')
#     ifndef Cyber
        call flush(6)
#     endif
c
c     now create the mesh for the simulation
c
      call mkmesh( error )
      if ( error ) then
        write(6,'('' '',''*** Error! ***''/
     +            '' '',''    An error occurred in the mesh '',
     +                  ''selection!''/
     +            '' '',''    There is no final output to be '',
     +                  ''written, but the solution is halted.'')')
        goto 9999
      endif
c
c     If we are doing a Sutherland and Mock simultation we need to taper
c     Qss so that we don't have convergence problems.
c
      if ( symtyp .eq. 'Mock' ) then
        do 101 i=1,(xnodes/2 - 3)
  101     Qss(i) = 0.0
        do 102 i=(xnodes/2 + 3),xnodes
  102     Qss(i) = 0.0
      endif
c
c     Now construct the doping profile for the simulation, and the
c     relative dielectric constant for the device. This really does
c     the bulk of the set-up work.
c
      call mkdev( .true. )
c
c     We need to get the device up to the biasing conditions that
c     the operator wishes. This involves an initial guess and
c     improving that guess with bias stepping techniques.
c
      call bias( error )
      if ( error ) then
        goto 9999
      endif
c
c     Now that the bias is set on the device we can reset the time
c     counters and start to do the time iterations
c
      write(6,'('' '')')
      write(6,'('' '')')
      write(6,'('' '',''==== Shark   v2.0                      '',
     +   ''                                            REB 1988'',
     +   '' starting transient ===='')')
      write(6,'('' '')')
#     ifndef Cyber
        call flush(6)
#     endif
  100 tnow = tstart
      tcnt = tstart
      gooddt = 0

c
c     output the starting results
c
      write(6,'('' '',''Outputting the simulation results for t = '',
     +                f6.2,'' (psec).'')') tnow * 1.0e12
      if ( pulse ) then
c
c       Integrate n(x.y) to get Qn(x), and then Qtot
c
        Qtot = 0.0
        do 300 i=1,xnodes
          Qn(i) = 0.0
          do 200 k=1,ynodes
c				get the delta y for this point
            if ( k .eq. 1 ) then
              dely = hy(k)*Ldi
            elseif ( k .eq. ynodes ) then
              dely = hy(k-1)*Ldi
            else
              dely = 0.5*( hy(k) + hy(k-1) )*Ldi
            endif
  200       Qn(i) = Qn(i) + n(i,k)*nio*dely
c				get the delta x for this point
          if ( i .eq. 1 ) then
            delx = hx(i)*Ldi
          elseif ( i .eq. xnodes ) then
            delx = hx(i-1)*Ldi
          else
            delx = 0.5*( hx(i) + hx(i-1) )*Ldi
          endif
          Qtot = Qtot + Qn(i)*delx
  300   continue
c
c       Now write it out for a sanity check of the simulator
c
        write(6,'('' '',''  ...with total charge, Qtot, equalling '',
     +                  g10.4,'' particles = '',g11.4,'' Coul.'')')
     +        Qtot, -1.6022e-19*Qtot
      endif
#     ifndef Cyber
        call flush(6)
#     endif
      call dump( tnow, .true., xref )
c
c     Check to see if it is time to dump the core again, but not
c     to dump any plotting information. Do it every hour of run 
c     time.
c
  400 if ( int( etime(tarray)/(60*60) ) .gt. hours ) then
        hours = int( etime(tarray)/(60*60) )
#       ifndef Cyber
          call dump( tnow, .false., xref )
#       endif
      endif
c
c     solve for the next time step and check for errors in newton()
c
      call newton( error, iter )
c
c     Check to see if we had a coupling problem
c
      if ( iter .eq. -1 ) then
        call refine( .false., .true., .true., error )
        if ( error ) then
          write(6,'('' '',''*** Error! ***''/
     +              '' '',''    Shark could not resolve the '',
     +                    ''coupling problems. This problem is '',
     +                    ''not solvable''/
     +              '' '',''    with the current executable of '',
     +                    ''Shark. the final output is being '',
     +                    ''written and''/
     +              '' '',''    the solution halted.'')')
          write(9,'('' '',''#eor'')')
#         ifndef Cyber
            call flush(6)
            call dump( tnow, .true., xref )
            close(9)
#         endif
          call depmap( 1, 1 )
          goto 9999
        endif
c
c       No error in mesh refinement, so do all the time stepping
c       over again.
c
        write(6,'('' '',''...Restarting the time stepping...'')')
        goto 100
      endif
c
c     check to see if dt was too big
c
      if ( error ) then
c				try to halve it to get it into
c				quadratic convergence
        dt = dt / 2.0
        if ( dt .lt. stpmin ) then
          write(6,'('' '',''*** Error! ***''/
     +              '' '',''    An error occurred in the solution '',
     +                    ''for t = '',f6.2,'' (psec). The Newton '',
     +                    ''iterations could not''/
     +              '' '',''    converge. The final output is '',
     +                    ''being written and the solution '',
     +                    ''halted.'')') tnow*1.0e12
#         ifndef Cyber
            call flush(6)
            call dump( tnow, .true., xref )
            close(9)
#         endif
          call depmap( 1, 1 )
          goto 9999
        endif
c				not too small, so try again
        if ( debug .ne. 'no' ) then
          write(6,'('' '',''...error at present dt, halving and '',
     +                    ''retrying...'')')
#         ifndef Cyber
            call flush(6)
#         endif
        endif
        gooddt = 0
        goto 400
      else
c
c       good dt, so count it
c
        gooddt = gooddt + 1
      endif
c
c     increment the time counters
c
      tnow = tnow + dt
      tcnt = tcnt + dt
c
c     See if we are doing a running pulse simulation, and if it
c     is time to shift the coordinate system
c
      if ( pulse ) then
c
c       Get the important values of Qn
c
        ibuff = min( (xnodes/8), 5 )
        Qn( xnodes-ibuff ) = 0.0
        do 500 k=1,ynodes
c				get the delta y for this point
          if ( k .eq. 1 ) then
            dely = hy(k)*Ldi
          elseif ( k .eq. ynodes ) then
            dely = hy(k-1)*Ldi
          else
            dely = 0.5*( hy(k) + hy(k-1) )*Ldi
          endif
  500     Qn( xnodes-ibuff ) = Qn( xnodes-ibuff ) +
     +                         n( (xnodes-ibuff), k )*nio*dely

        Qn( xnodes ) = 0.0
        do 600 k=1,ynodes
c				get the delta y for this point
          if ( k .eq. 1 ) then
            dely = hy(k)*Ldi
          elseif ( k .eq. ynodes ) then
            dely = hy(k-1)*Ldi
          else
            dely = 0.5*( hy(k) + hy(k-1) )*Ldi
          endif
  600     Qn( xnodes ) = Qn( xnodes ) + n( xnodes, k )*nio*dely
c
c       Now test to see if the pulse is approaching the right side
c
        if ( Qn( xnodes-ibuff ) .gt. 1.0 ) then
c
c         Find out exactly how much we can shift
c
          do 800 i=1,xnodes
            Qn(i) = 0.0
            do 700 k=1,ynodes
c				get the delta y for this point
              if ( k .eq. 1 ) then
                dely = hy(k)*Ldi
              elseif ( k .eq. ynodes ) then
                dely = hy(k-1)*Ldi
              else
                dely = 0.5*( hy(k) + hy(k-1) )*Ldi
              endif
  700         Qn(i) = Qn(i) + n(i,k)*nio*dely
  800     continue

          deltai = 0
          deltax = 0.0
  900     if ( Qn(deltai+ibuff) .lt. 1.0 ) then
            deltai = deltai + 1
            deltax = deltax + hx(deltai)*Ldi
            goto 900
          endif
c
c         shift the gate voltage and the coordinate axis
c
          do 1000 i=1,xnodes
            if ( vgate(i) .lt. 1.0e20 ) then
              vgate(i) = vgate(i) + eavg*deltax*Ldi/kT
            endif
 1000     continue
          xref = xref + deltax
c
c         Shift the voltage and carrier concentrations
c
          do 1100 i=1,(xnodes-deltai)
            do 1100 k=1,ynodes
              n(i,k) = n(i+deltai,k)
              V(i,k) = V(i+deltai,k)
 1100     continue
c
c         Pad the unknown section with the boundary condition and
c         run it through the 1D electrostatics to update it to the
c         new gate voltage.
c
          do 1300 i=(xnodes-deltai+1), xnodes
c
c           Copy the last good one to where I am now, for an initial
c           guess.
c
            do 1200 k=1,ynodes
              n(i,k) = n(i-1,k)
              V(i,k) = V(i-1,k)
 1200       continue
c
c           Update the present line to be correct with the new gate
c           voltage, but we should be OK as far as the electron
c           concentration goes.
c
            call oneDes( i, .true., 0.0, error, .false., iter )
            call full1D( i, error, .false., iter )
 1300     continue
c
c         Now we need to copy into srccon(j) and drncon(j) from the
c         shifted V and n
c
          do 1400 k=1,ynodes
            srccon(k,1) = v(1,k)
            srccon(k,2) = n(1,k)
            drncon(k,1) = v(xnodes,k)
            drncon(k,2) = n(xnodes,k)
 1400     continue
        endif
      endif
c
c     See if we have enough good dt's to try and bump dt up, but not
c     too far, and reset the number of good iterations. Do this only
c     if we are not using a direct solver.
c
      if ( gooddt .ge. bumpdt ) then
        if ( chtyp .eq. 'Si' ) then
          dt = min( dt*2.0, 10.0*dtmax )
        elseif ( chtyp .eq. 'GaAs' ) then
          dt = min( dt*2.0, dtmax )
        endif
        if ( debug .ne. 'no' ) then
          write(6,'('' '',''...at present dt, '',i2,
     +                    '' successful steps done, so trying '',
     +                    ''to double dt...'')') gooddt
#         ifndef Cyber
            call flush(6)
#         endif
        endif
        gooddt = 0
      endif
c
c     See if it is time to output the status of the simulator
c
      if ( tcnt .ge. tout ) then
        tcnt = 0.0
        write(6,'('' '',''Outputting the simulation results for t = '',
     +                  f6.2,'' (psec).'')') tnow*1.0e12
        if ( pulse ) then
c
c         Integrate n(x,y) to get Qn(x), and then Qtot
c
          Qtot = 0.0
          do 1600 i=1,xnodes
            Qn(i) = 0.0
            do 1500 k=1,ynodes
c				get the delta y for this point
              if ( k .eq. 1 ) then
                dely = hy(k)*Ldi
              elseif ( k .eq. ynodes ) then
                dely = hy(k-1)*Ldi
              else
                dely = 0.5*( hy(k) + hy(k-1) )*Ldi
              endif
 1500         Qn(i) = Qn(i) + n(i,k)*nio*dely
c				get the delta x for this point
            if ( i .eq. 1 ) then
              delx = hx(i)*Ldi
            elseif ( i .eq. xnodes ) then
              delx = hx(i-1)*Ldi
            else
              delx = 0.5*( hx(i) + hx(i-1) )*Ldi
            endif
            Qtot = Qtot + Qn(i)*delx
 1600     continue
c
c         Now write it out for a sanity check of the simulator
c
          write(6,'('' '',''  ...with total charge, Qtot, equalling '',
     +                    g10.4,'' particles = '',g11.4,'' Coul.'')')
     +          Qtot, -1.6022e-19*Qtot
        endif
#       ifndef Cyber
          call flush(6)
          call dump( tnow, .true., xref )
#       endif
      endif
c
c     See if we still have some time to go before we stop
c
      if ( tnow .lt. tstop ) goto 400
c
c     If tcnt <> 0, then output the final state of the simulator
c
      if ( tcnt .ne. 0.0 ) then
        tcnt = 0.0
        write(6,'('' '',''Outputting the simulation results for t = '',
     +                  f6.2,'' (psec).'')') tnow*1.0e12
        if ( pulse ) then
c
c         Integrate n(x,y) to get Qn(x), and then Qtot
c
          Qtot = 0.0
          do 1800 i=1,xnodes
            Qn(i) = 0.0
            do 1700 k=1,ynodes
c				get the delta y for this point
              if ( k .eq. 1 ) then
                dely = hy(k)*Ldi
              elseif ( k .eq. ynodes ) then
                dely = hy(k-1)*Ldi
              else
                dely = 0.5*( hy(k) + hy(k-1) )*Ldi
              endif
 1700         Qn(i) = Qn(i) + n(i,k)*nio*dely
c				get the delta x for this point
            if ( i .eq. 1 ) then
              delx = hx(i)*Ldi
            elseif ( i .eq. xnodes ) then
              delx = hx(i-1)*Ldi
            else
              delx = 0.5*( hx(i) + hx(i-1) )*Ldi
            endif
            Qtot = Qtot + Qn(i)*delx
 1800     continue
c
c         Now write it out for a sanity check of the simulator
c
          write(6,'('' '',''  ...with total charge, Qtot, equalling '',
     +                    g10.4,'' particles = '',g11.4,'' Coul.'')')
     +          Qtot, -1.6022e-19*Qtot
        endif
#       ifndef Cyber
          call flush(6)
          call dump( tnow, .true., xref )
#       endif
      endif
c
c     end of Shark... close everything up and exit
c
      write(9,'('' '',''#eor'')')
c
c     print the elapsed time out for this run
c
      time = etime( tarray )
c
c     compute the breakdown of the total time
c
      days = int( time/(60*60*24) )
      hours = int( (time-days*60*60*24)/(60*60) )
      mins = int( (time - (hours+days*24)*60*60)/60 )
      secs = int( time - ((hours+days*24)*60+mins)*60 )
      if ( days .ge. 0 ) then
        write(6,'('' '','' Total CPU time used: '',i2,'' days, '',
     +                  i2.2,'':'',i2.2,'':'',i2.2)') 
     +               days, hours, mins, secs
      else
        write(6,'('' '','' Total CPU time used: '',
     +                  i2.2,'':'',i2.2,'':'',i2.2)') 
     +               hours, mins, secs
      endif
 9999 continue
#     ifndef Cyber
        call flush(6)
#     endif
      stop
      end

#ifdef Cyber
c
c     These timer routines are Unix equivalents for the Cyber 205.
c     They are here for compatibility sake.
c
      real function dtime( tar )
      real tar(2)
c
c     This function returns the elapsed time since the last call to
c     dtime(), or the beginning of execution on the first call
c
      common /time/ first, last, now
      real first, last, now
c
c     Get the time
c
      now = second()
      if ( first .eq. 0.0 ) first = now
      dtime = now - last
      last = now
      return
      end

      real function etime( tar )
      real tar(2)
c
c     This function returns the elapsed time since the beginning
c     of execution.
c
      common /time/ first, last, now
      real first, last, now
c
c     Get the time
c
      etime = second() - first
      return
      end
#endif

