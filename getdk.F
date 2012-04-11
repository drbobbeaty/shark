#     ifdef Cyber
      subroutine getdk( inperr, restart )
#     else
      subroutine getdk( inperr )
#     endif
c
c     This subroutine will read the input deck and update the parameters
c     from the defaults. This is simply calling the parser and decoding
c     it to the possible valid cards.
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
c     Get the input logical unit number based on the type of computer
c     being compiled and run on
c
      integer inlun
#     ifdef Cyber
        parameter( inlun=5 )
#     else
        parameter( inlun=8 )
#     endif
c
c     now for the arguments to this routine
c
      integer inperr
#     ifdef Cyber
        logical restart
#     endif
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
c     now for the common blocks for interface to the parser
c
      common /parsd1/ kflag,numr,numc,vreal(70)
      integer kflag, numr, numc
      real vreal
      common /parsd2/ ncard,nreal(70),nchar(30),
     +                vchar(30),ititle(10),ibuf(400)
      character*8 ncard,nreal,nchar,vchar,ititle
      character*1 ibuf
c
c     now for some local variables
c
      integer i, j
c
c     Valid input deck cards:
c     col 1.
c     !
c     * this is a title card and will appear on output
c     $ this is a comment card
c     source length=###(um) thick=###(um) doping=###(/cm^3) overlap=###(um)
c     + corner=____(square, sloped, curved)
c     drain length=###(um) thick=###(um) doping=###(/cm^3) overlap=###(um)
c     + corner=____(square, sloped, curved)
c     oxide thick=###(A) k=###(relative) type=____(SiN, SiO2, AlGaAs, AlAs)
c     + mobil=###(cm^2/V-sec) Qss=###(/cm^2)
c     channel thick=###(um) type=____(MODFET, buried) doping=###(/cm^3)
c     + doptyp=____(step, tailed) length=###(um) material=____(GaAs, Si)
c     bulk thick=###(um) k=###(relative) doping=###(/cm^3)
c     mesh xnodes=### xfrac=### ynodes=### auto=yes/no joxide=### dVa=###(kT)
c     time tstart=###(psec.) tstop=###(psec.) stpmin=###(psec.)
c     + tout=###(psec.) delt=###(psec.)
c     bias vsource=###(V) vdrain=###(V) eavg=###(V/cm) vg1=###(V)
c     + delV=###(V) delE=###(V/cm) delT=###(psec.)
c     VvsE mobility=###(cm^2/V-sec) A=###(cm/sec) B=### C=###(cm/sec)
c     + D=### F=###
c     DvsE Do=###(cm^2/sec) dDo=###(cm^2/sec) A=### B=###
c     grate line=###(um) space=###(um) npulse=###(/cm^3) defoc=###
c     pulse Qm=###(/cm^2) flat=###(um) fwhm=###(um)
c     solution type=____(ortho, minres, linpk, gauss)
c     + keep=#(0->6)
c     start with=____(zero, DC) vstop=###(V) maxit=###
c     + debug=yes/no guess=____(MOS, charge)
c
c     On the Cyber ONLY:
c     restart
c
c     start with a header so the the output looks reasonable
c
      write(6,'('' '')')
      write(6,'('' '')')
      write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                            REB 1988'',
     +     ''         input deck ===='')')
      write(6,'('' '')')
c
c     Assume that we have no source or drain cards in the deck, and we
c     will be doing the standard center simultation.
c
      symtyp = 'center'
c
c     Assume that we are not going to be doing a pulse
c
      pulse = .false.

#     ifdef Cyber
c
c     Reset the state of the restart flag for the Cyber ONLY
c
      restart = .false.
#     endif

c
c     now get a line from the input deck
c
 1000 call parser(inlun, inperr)
c     ...see if we are at the end of the input deck
      if ( kflag .eq. -1 ) then
#       ifdef Cyber
          close(inlun)
#       endif
        return
      endif
c     ...see if there was an input syntax error in the card just read
      if ( inperr .ne. 0 ) then
        write(6,'('' '',''*** Syntax error found in input deck!'')')
        write(6,'('' '',''Shark aborted!'')')
#       ifndef Cyber
          call flush(6)
#       endif
        return
      endif
c
c     interpret the read card
c
c     ...is it a title card?
c
      if ( ncard .eq. '*' ) then
c				don't do anything for now
        ncard = ncard
      endif

#     ifdef Cyber
c
c     ...is it a restart card?
c
      if ( ncard .eq. 'restart' ) then
        restart = .true.
      endif
#     endif

c
c     ...is it a source card?
c
      if ( ncard .eq. 'source' ) then
        symtyp = 'Mock'
        do 1010 i=1,numr
          if ( nreal(i) .eq. 'length' ) then
            lsorce = vreal(i) * 1.0e-4
            if ( lsorce .gt. 0.0 ) symtyp = 'full'
          endif
          if ( nreal(i) .eq. 'overlap' ) gsol = vreal(i) * 1.0e-4
          if ( nreal(i) .eq. 'thick' ) tsorce = vreal(i) * 1.0e-4
 1010     if ( nreal(i) .eq. 'doping' ) dopsrc = vreal(i)
        do 1015 i=1,numc
 1015     if ( nchar(i) .eq. 'corner' ) srccor = vchar(i)
      endif
c
c     ...is it a drain card?
c
      if ( ncard .eq. 'drain' ) then
        symtyp = 'Mock'
        do 1020 i=1,numr
          if ( nreal(i) .eq. 'length' ) then
            ldrain = vreal(i) * 1.0e-4
            if ( ldrain .gt. 0.0 ) symtyp = 'full'
          endif
          if ( nreal(i) .eq. 'overlap' ) gdol = vreal(i) * 1.0e-4
          if ( nreal(i) .eq. 'thick' ) tdrain = vreal(i) * 1.0e-4
 1020     if ( nreal(i) .eq. 'doping' ) dopdrn = vreal(i)
        do 1025 i=1,numc
 1025     if ( nchar(i) .eq. 'corner' ) drncor = vchar(i)
      endif
c
c     ...is it an oxide card?
c
      if ( ncard .eq. 'oxide' ) then
        do 1030 i=1,numc
          if ( nchar(i) .eq. 'type' ) then
            instyp = vchar(i)
            if ( instyp .eq. 'SiN' ) then
              kox = 7.5
              uins = 1.0e-3
            elseif ( instyp .eq. 'SiO2' ) then
              kox = 3.9
              uins = 1.0e-4
            elseif ( instyp .eq. 'AlGaAs' ) then
              kox = 12.2
              uins = 100.0
            elseif ( instyp .eq. 'AlAs' ) then
              kox = 10.1
              uins = 10.0
            else
              write(6,'('' '',''  Error! Do not have data on '',
     +                        a6,'' as an insulator. Defaulting '',
     +                        ''to SiN'')') instyp
#             ifndef Cyber
                call flush(6)
#             endif
              kox = 7.5
              uins = 1.0e-3
              instyp = 'SiN'
            endif
          endif
 1030   continue
        do 1035 i=1,numr
          if ( nreal(i) .eq. 'thick' ) tox = vreal(i) * 1.0e-8
          if ( nreal(i) .eq. 'mobil' ) uins = vreal(i)
          if ( nreal(i) .eq. 'k' ) kox = vreal(i)
          if ( nreal(i) .eq. 'Qss' ) then
            do 1032 j=1,Nx
 1032         Qss(j) = vreal(i)
          endif
 1035   continue
      endif
c
c     ...is it a channel card?
c
      if ( ncard .eq. 'channel' ) then
        do 1040 i=1,numr
          if ( nreal(i) .eq. 'doping' ) dopch = vreal(i)
          if ( nreal(i) .eq. 'length' ) length = vreal(i) * 1.0e-4
 1040     if ( nreal(i) .eq. 'thick' ) tepi = vreal(i) * 1.0e-4
        do 1050 i=1,numc
          if ( nchar(i) .eq. 'type' ) then
            chmod = vchar(i)
c				if this is a MODFET channel then tepi << 1
            if ( chmod .eq. 'MODFET' ) then
              tepi = 0.0
            endif
          endif
          if ( nchar(i) .eq. 'doptyp' ) doptyp = vchar(i)
          if ( nchar(i) .eq. 'material' ) chtyp = vchar(i)
 1050   continue
        if ( chtyp .eq. 'Si' ) then
          nio = 1.5e10
          ksemi = 11.7
c				this is the low field diffusivity for Si
          Do = 18.13
c				this is the change, at peak, from Do for Si
          dDo = 0.0
c				these are some more parameters for D vs. E
          Ad = 8.13
          Bd = 0.596
c				this is the low field mobility for Si
          mobil = 700.0
c				these are some parameters for vel vs. E curve
          Av = 9.23e6
          Bv = 1.0e38
          Cv = 0.0
          Dv = 1.0
          Fv = 1.0
        elseif ( chtyp .eq. 'GaAs' ) then
          nio = 1.79e6
          ksemi = 12.9
c				this is the low field diffusivity for GaAs
          Do = 220.0
c				this is the change, at peak, from Do for GaAs
          dDo = 624.0
c				these are some more parameters for D vs. E
          Ad = 8.13
          Bd = 0.596
c				this is the low field mobility for GaAs
          mobil = 8400.0
c				these are some parameters for vel vs. E curve
          Av = 9.546e7
          Bv = 1644.0
          Cv = 6.473e7
          Dv = 130.5
          Fv = 0.3164
        endif
      endif
c
c     ...is it a substrate (bulk) card?
c
      if ( ncard .eq. 'bulk' ) then
        do 1060 i=1,numr
          if ( nreal(i) .eq. 'thick' ) tsub = vreal(i) * 1.0e-4
          if ( nreal(i) .eq. 'k' ) ksemi = vreal(i)
 1060     if ( nreal(i) .eq. 'doping' ) dopsub = vreal(i)
      endif
c
c     ...is it a mesh card?
c
      if ( ncard .eq. 'mesh' ) then
        do 1070 i=1,numr
          if ( nreal(i) .eq. 'xnodes' ) xnodes = int(vreal(i))
          if ( nreal(i) .eq. 'xfrac' ) xfrac = vreal(i)
          if ( nreal(i) .eq. 'ynodes' ) ynodes = int(vreal(i))
          if ( nreal(i) .eq. 'joxide' ) joxide = int(vreal(i))
 1070     if ( nreal(i) .eq. 'dVa' ) dVa = vreal(i)
        do 1075 i=1,numc
 1075     if ( nchar(i) .eq. 'auto' ) auto = vchar(i)
      endif
c
c     ...is it a time card?
c
      if ( ncard .eq. 'time' ) then
        do 1080 i=1,numr
          if ( nreal(i) .eq. 'stpmin' ) stpmin = vreal(i) * 1.0e-12
          if ( nreal(i) .eq. 'delt' ) dt = vreal(i) * 1.0e-12
          if ( nreal(i) .eq. 'tstart' ) tstart = vreal(i) * 1.0e-12
          if ( nreal(i) .eq. 'tstop' ) tstop = vreal(i) * 1.0e-12
 1080     if ( nreal(i) .eq. 'tout' ) tout = vreal(i) * 1.0e-12
      endif
c
c     ...is it a bias point card?
c
      if ( ncard .eq. 'bias' ) then
        do 1090 i=1,numr
          if ( nreal(i) .eq. 'vsource' ) Vs = vreal(i)
          if ( nreal(i) .eq. 'vdrain' ) Vd = vreal(i)
          if ( nreal(i) .eq. 'vg1' ) vg1 = vreal(i)
          if ( nreal(i) .eq. 'mobility' ) mobil = vreal(i)
          if ( nreal(i) .eq. 'eavg' ) then
            eavg = abs( vreal(i) )
            if ( vreal(i) .lt. 0.0 ) then
              write(6,'('' '',''*** Warning! *** eavg has to be'',
     +           '' greater than zero for this device''/
     +                  '' '',''    to function properly. Shark'',
     +           '' has taken the absolute value.'')')
            endif
          endif
          if ( nreal(i) .eq. 'delV' ) delV = vreal(i)
          if ( nreal(i) .eq. 'delT' ) delT = vreal(i) * 1.0e-12
 1090     if ( nreal(i) .eq. 'delE' ) delE = vreal(i)
      endif
c
c     ...is it a VvsE card?
c
      if ( ncard .eq. 'VvsE' ) then
        do 1100 i=1,numr
          if ( nreal(i) .eq. 'mobility' ) mobil = vreal(i)
          if ( nreal(i) .eq. 'A' ) Av = vreal(i)
          if ( nreal(i) .eq. 'B' ) Bv = vreal(i)
          if ( nreal(i) .eq. 'C' ) Cv = vreal(i)
          if ( nreal(i) .eq. 'D' ) Dv = vreal(i)
 1100     if ( nreal(i) .eq. 'F' ) Fv = vreal(i)
      endif
c
c     ...is it a DvsE card?
c
      if ( ncard .eq. 'DvsE' ) then
        do 1200 i=1,numr
          if ( nreal(i) .eq. 'Do' ) Do = vreal(i)
          if ( nreal(i) .eq. 'dDo' ) dDo = vreal(i)
          if ( nreal(i) .eq. 'A' ) Ad = vreal(i)
 1200     if ( nreal(i) .eq. 'B' ) Bd = vreal(i)
      endif
c
c     ...is it a grating card?
c
      if ( ncard .eq. 'grating' ) then
        do 1300 i=1,numr
          if ( nreal(i) .eq. 'line' ) line = vreal(i) * 1.0e-4
          if ( nreal(i) .eq. 'space' ) space = vreal(i) * 1.0e-4
          if ( nreal(i) .eq. 'npulse' ) npulse = vreal(i)
 1300     if ( nreal(i) .eq. 'defoc' ) defoc = vreal(i)
      endif
c
c     ...is it a pulse card?
c
      if ( ncard .eq. 'pulse' ) then
        do 1400 i=1,numr
          if ( nreal(i) .eq. 'Qm' ) Qm = vreal(i)
          if ( nreal(i) .eq. 'flat' ) ftw = vreal(i) * 1.0e-4
 1400     if ( nreal(i) .eq. 'fwhm' ) fwhm = vreal(i) * 1.0e-4
        pulse = .true.
      endif
c
c     ...is it a solution card?
c
      if ( ncard .eq. 'solution' ) then
        do 1500 i=1,numc
 1500     if ( nchar(i) .eq. 'type' ) soltyp = vchar(i)
        do 1600 i=1,numr
          if ( nreal(i) .eq. 'keep' ) then
            keep = int(vreal(i))
            if ( keep .eq. 0 ) soltyp = 'minres'
          endif
 1600   continue
      endif
c
c     ...is it a start card?
c
      if ( ncard .eq. 'start' ) then
        do 1700 i=1,numc
          if ( nchar(i) .eq. 'debug' ) debug = vchar(i)
          if ( nchar(i) .eq. 'guess' ) guess = vchar(i)
 1700     if ( nchar(i) .eq. 'with' ) icchan = vchar(i)
        do 1800 i=1,numr
          if ( nreal(i) .eq. 'vstop' ) vmax = vreal(i)
 1800     if ( nreal(i) .eq. 'maxit' ) maxit = int(vreal(i))
      endif
c
c     ...is it the data card for an initial specified electron density?
c
      if ( ncard .eq. 'data' ) then
        write(6,'('' '',''  --> reading in old data points <--'')')
      endif
c
c     ...and keep going until all the cards have been read
c
      goto 1000
c
c     and put an end here to end this procedure
c
      end
