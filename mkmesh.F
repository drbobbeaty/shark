      subroutine mkmesh( error )
c
c     This subroutine constructs the mesh for the simulation based
c     on the distribution scheme of ynodes distributed uniformly down
c     the length of the channel, and nins nodes in the insulator,
c     nsub nodes in the substrate, and xnodes in the epi-layer.
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
      PRECISION hxchan, hxnpl, Mfact
      parameter( hxchan=0.10e-4, hxnpl=0.10e-4, Mfact=2.00 )
c
c     This parameter is for the number of lengths to be used in the
c     pulse simulation. This is the amount the shifting reference
c     frame can go and still be safe at the back bulk boundary condition.
c
      PRECISION noleng
      parameter( noleng=1.5 )
c
c     These parameters are for the location of nodes in the inversion
c     layer. For the surface devices jinv, is the number of nodes in
c     the first dVa, jbuf is the number of nodes in the second dVa.
c     For the buried channel mesh, jchn is the number of nodes per
c     dVa placed in the channel region.
c
      integer jinv, jbuf, jchn
      parameter( jinv=8, jbuf=4, jchn=4 )
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
      common /Pshape/ Qm, ftw, fwhm, pulse
      PRECISION Qm, ftw, fwhm
      logical pulse
c
c     local variables
c
      integer i, j, jtemp
      PRECISION hyox, hyrat, hxrat, Vgmax, hxmax, hxmin, Cox,
     +          Vo, phis, xd, hymax, hymin, Va, V1, xm, phim,
     +          phit, K1, Qssmax
c
c     Preset a few local variables so the compiler doesn't flake out
c
      K1 = 0.0
c
c
c     select the mesh on the scheme:
c
c                               o T
c                                 hy(j-1)
c                              i,j
c            L o                o          o R
c                  hx(i-1)          hx(i)
c
c                                 hy(j)
c                               o B
c
c
c     start with a header so that the output looks reasonable
c
      write(6,'('' '')')
      write(6,'('' '')')
      write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                            REB 1988'',
     +     ''     mesh selection ===='')')
      write(6,'('' '')')
c
c     Convert the delta V to kT for the rest of the program
c
      dVa = dVa * kT
c
c     If we have a MODFET channel make sure that tepi=0
c
      if ( chmod .eq. 'MODFET' ) then
        tepi = 0.0
      endif
c
c     Let's do a little MOS electrostatics to see how deep the depletion
c     region might go.
c
      xd = 0.0
      if ( symtyp .eq. 'full' ) then
        if ( auto .ne. 'yes' ) then
          Vgmax = Vg1 + eavg*(length-(lsorce-gsol)-(ldrain-gdol))
        else
          Vgmax = Vg1 + eavg*(length + gdol + gsol)
        endif
      else
        if ( pulse ) then
          Vgmax = Vg1 + eavg*(noleng*length)
        else
          Vgmax = Vg1 + eavg*length
        endif
      endif
      if ( (Vgmax .gt. 0.0) .or.
     +     ( (Vgmax .ge. 0.0) .and. (chmod .eq. 'buried') ) ) then
        Cox = kox * 8.854e-14 / tox
        Vo = 1.6022e-19 * ksemi * 8.854e-14 * dopsub / (Cox*Cox)
        V1 = (1.6022e-19*tepi/Cox)*( dopsub + dopch )
c				get the positive maximum of Qss
        Qssmax = -1.0e30
        do 10 i=1,xnodes
   10     Qssmax = max( Qssmax, Qss(i) )
        Vgmax = Vgmax + 1.6022e-19*Qssmax/Cox
        phis = Vgmax + Vo + V1 - sqrt( Vo*Vo + 2.0*Vo*( Vgmax + V1 )
     +                   + V1*1.6022e-19*tepi*dopsub/Cox )
        xd = sqrt( tepi*tepi*(1.0 + dopch/dopsub)
     +             + 2.0*ksemi*8.854e-14*phis/(1.6022e-19 * dopsub) )
        write(6,'('' '',''The simple MOS electrostatics for the '',
     +                  ''channel region are:'')')
        write(6,'('' '',''      Vgmax = '',f6.3,'' V'')')
     +                             Vgmax - 1.6022e-19*Qssmax/Cox
        write(6,'('' '','' Vg+Qss/Cox = '',f6.3,'' V'')') Vgmax
        write(6,'('' '',''        Cox = '',g9.3,'' F/cm**2'')') Cox
        write(6,'('' '',''         Vo = '',f6.4,'' V'')') Vo
        write(6,'('' '',''         V1 = '',f6.4,'' V'')') V1
        write(6,'('' '',''       phis = '',f6.3,'' V'')') phis
        if ( chmod .eq. 'buried' ) then
          xm = tepi - (dopsub/dopch)*(xd-tepi)
c				see if we need just the extremum
          if ( xm .lt. 0.0 ) xm = 0.0
          phim = (1.6022e-19*dopsub/(2.0*ksemi*8.854e-14))*
     +           ( xd*xd - tepi*tepi - 2.0*xm*(xd-tepi) ) -
     +           (1.6022e-19*dopch/(2.0*ksemi*8.854e-14))*
     +           ( tepi - xm )*( tepi - xm )
          phit = (1.6022e-19*dopsub/(2.0*ksemi*8.854e-14))*
     +           ( xd - tepi )*( xd - tepi )
          write(6,'('' '',''       phim = '',f6.3,'' V'')') phim
          write(6,'('' '',''         xm = '',f6.4,'' um'')') xm*1.0e4
          write(6,'('' '',''       phit = '',f6.3,'' V'')') phit
          K1 = (dopsub/dopch)*(dopsub/dopch + 1)*( tepi - xd )**2
        endif
      endif
c				limit it on the low end
      xd = max( xd, tox )
      write(6,'('' '',''         xd = '',f6.3,'' um'')') xd*1.0e4
      write(6,'('' '')')
c
c     Now, if we are selecting the mesh (we should be) then move tsub
c     according to a bit greater than the maximum depletion region
c
      if ( auto .eq. 'yes' ) then
c
c       Is is a center simulation?
c
        if ( (symtyp .eq. 'center') .and. (.not. pulse) ) then
          if ( chmod .eq. 'MODFET' ) then
c				good margin for the MODFET channel
            if ( chtyp .eq. 'Si' ) then
              tsub = 3.0 * xd
            else
              tsub = 2.0 * xd
            endif
          else
c				less margin for the buried channel
            tsub = 1.5 * xd
          endif
c
c       ...is it a pulse center simulation?
c
        elseif ( pulse ) then
          if ( chmod .eq. 'MODFET' ) then
c				good margin for the MODFET channel
            if ( chtyp .eq. 'Si' ) then
              tsub = 2.5 * xd
            else
              tsub = 3.0 * xd
            endif
          else
c				less margin for the buried channel
            tsub = 1.25 * xd
          endif
c
c       ...is it a Sutherland and Mock-type simulation?
c
        elseif ( symtyp .eq. 'Mock' ) then
          if ( chmod .eq. 'MODFET' ) then
c				good margin for the MODFET channel
            if ( chtyp .eq. 'Si' ) then
              tsub = 8.0 * ( xd + max(tsorce, tdrain) )
            else
              tsub = 8.0 * ( xd + max(tsorce, tdrain) )
            endif
          else
c				less margin for the buried channel
            tsub = 7.0 * ( xd + max(tsorce, tdrain) )
          endif
c
c       ...is it a full device simulation?
c
        elseif ( symtyp .eq. 'full' ) then
          tsub = 1.75 * ( xd + 3.52058*max(tsorce, tdrain) )
        endif
c				now write it out
        write(6,'('' '',''Because Shark is selecting the mesh, Shark'',
     +                  '' is setting the thickness of the''/
     +            '' '',''substrate to '',f6.3,'' um'')') tsub*1.0e4
        write(6,'('' '')')
      endif
c
c     If we are supposed to make this mesh ourselves, then pick
c     a rough mesh based on the type of simulation being done.
c     But keep the x mesh separate from the y mesh.
c
      if ( auto .eq. 'yes' ) then
c
c       Do we need a y mesh for a MODFET channel?
c
        if ( chmod .eq. 'MODFET' ) then
          write(6,'('' '',''   automatically selecting the y mesh '',
     +                    ''for a MODFET channel.'')')
          write(6,'('' '')')
          Va = phis - dVa * int( phis/dVa )
          if ( (Va .gt. 0.0) .and. (phis .ge. 2.0*dVa) ) then
c
c           place the first four nodes in dVa, the second two in dVa,
c           and the rest dVa apart.
c
            ynodes = joxide + int( (phis/dVa)-2 ) + jinv + jbuf
     +             + nint( (tsub - xd)/
     +             ( sqrt( 2.0*ksemi*8.854e-14/(1.6022e-19*dopsub) )*
     +             ( sqrt( Va+dVa ) - sqrt( Va ) ) ) )
c				set a minimum on ynodes of 3*joxide
          else
            ynodes = 3*joxide
          endif
c
c       See if we need a y mesh for a buried channel
c
        elseif ( chmod .eq. 'buried' ) then
          write(6,'('' '',''   automatically selecting the y mesh '',
     +                    ''for a MESFET, or buried, channel.'')')
          write(6,'('' '')')
          Va = phit - dVa*( int( phit/dVa )-1 )
          ynodes = joxide + int( jchn*(phim-phis)/dVa ) +
     +             int( jchn*(phim-phit)/dVa ) + int( phit/dVa ) +
     +             nint( (tsub - xd)/( sqrt( 2.0*ksemi*8.854e-14/
     +             (1.6022e-19*dopsub) )*( sqrt( Va ) -
     +             sqrt( Va-dVa ) ) ) )
        endif
c
c       See if my chosen ynodes is too big
c
        if ( ynodes .gt. Ny ) then
          write(6,'('' '',''**** Error! ****''/
     +              '' '',''     You have specified a device that is '',
     +                    ''larger than Shark can simulate''/
     +              '' '',''     with assurance of solution. Please '',
     +                    ''choose a smaller gate bias,''/
     +              '' '',''     or increase the bulk doping. Shark '',
     +                    ''cannot contain '',i4,'' ynodes, only '',
     +                    i4,''.'')') ynodes, Ny
          stop
        endif
c
c       See if we need an x mesh for a full device simulation
c
        if ( symtyp .eq. 'full' ) then
          xnodes = int( (lsorce + 2.52058*tsorce)/hxnpl ) +
     +             int(length/hxchan) +
     +             int( (ldrain + 2.52058*tdrain)/hxnpl ) + 1
c
c       See if we need an x mesh for a center channel simulation
c
        elseif ( symtyp .eq. 'center' ) then
          xnodes = int( length/hxchan ) + 1
c
c       See if we need an x mesh for a Mock simulation, which will
c       use the Mock selection method.
c
        elseif ( symtyp .eq. 'Mock' ) then
c
c         get the appropriate number of nodes we will need
c
          K1 = 0.9
          i = 1
  100     if ( (K1 .ge. (0.5*(length/(Mfact*hxchan) - 1.0))) .or.
     +         (i .ge. (Nx/2 - 1)) ) then
            xnodes = 2*( i + 1 )
          else
            i = i + 1
            K1 = K1 + (0.9)**i
            goto 100
          endif
c
c         Now compute the exact value of hx(xnodes/2)
c
          hx( xnodes/2 ) = length/( 1.0 + 2.0*K1 )/Ldi
        endif
c
c       Now see if we went too far in the node distribution for x
c
        if ( xnodes .gt. Nx ) then
          write(6,'('' '',''**** Error! ****''/
     +              '' '',''     You have specified a device that is '',
     +                    ''larger than Shark can simulate''/
     +              '' '',''     with assurance of solution. Please '',
     +                    ''choose a smaller channel length.''/
     +              '' '',''     Shark '',
     +                    ''cannot contain '',i4,'' xnodes, only '',
     +                    i4,''.'')') xnodes, Nx
          stop
        endif
c
c       Done with a rough node count selection!
c
      endif
c
c     Now do the mesh in the y direction, depending on the type of
c     channel we are supposed to be simulating.
c
c     Is it a MODFET channel?
c
      if ( chmod .eq. 'MODFET' ) then
        hyox = tox / joxide
        write(6,'('' '',''   placing y nodes uniformly in tox...'')')
        do 200 j=1, (joxide-1)
  200     hy(j) = hyox/Ldi
        if ( phis .ge. 2.0*dVa ) then
          write(6,'('' '',''   placing y nodes quadratically in '',
     +                    ''xd...'')')
          do 300 j=0, (jinv-1)
            Va = phis - (dVa/jinv) * j
  300       hy( j+joxide ) = sqrt( 2.0*ksemi*8.854e-14/(1.6022e-19*
     +                       dopsub) )*( sqrt( Va+(dVa/jinv) ) - 
     +                       sqrt( Va ) )/Ldi
          do 310 j=0, (jbuf-1)
            Va = phis - dVa - (dVa/jbuf) *j
  310       hy( j+joxide+jinv ) = sqrt( 2.0*ksemi*8.854e-14/
     +                       (1.6022e-19*dopsub) )*( sqrt( Va+
     +                       (dVa/jbuf) ) - sqrt( Va ) )/Ldi
          do 320 j=0, ( int( (phis/dVa)-2 ) - 1 )
            Va = phis - 2.0*dVa - dVa * j
  320       hy( j+joxide+jinv+jbuf ) = sqrt( 2.0*ksemi*8.854e-14/
     +                      (1.6022e-19*dopsub) )*( sqrt( Va+dVa ) 
     +                      - sqrt( Va ) )/Ldi
        endif
        write(6,'('' '',''   placing y nodes uniformly past xd '',
     +                  ''to tsub...'')')
        do 400 j=( joxide+jinv+jbuf+int((phis/dVa)-2) ), ( ynodes-1 )
  400     hy( j ) = hy( joxide+jinv+jbuf+int((phis/dVa)-2) - 1 )
c
c     Is it a MESFET, or buried, channel?
c
      elseif ( chmod .eq. 'buried' ) then
        hyox = tox / joxide
        write(6,'('' '',''   placing y nodes uniformly in tox...'')')
        do 500 j=1, (joxide-1)
  500     hy(j) = hyox/Ldi
        jtemp = joxide
        write(6,'('' '',''   placing y nodes quadratically to '',
     +                  ''xm...'')')
        do 600 j=0, (int( jchn*(phim-phis)/dVa )-1)
          Va = phis + (dVa/jchn) * j
  600     hy( j+jtemp ) = ( -sqrt( K1 - 2.0*ksemi*8.854e-14*( Va+dVa )
     +                    /(1.6022e-19*dopch) ) + sqrt( K1 - 2.0*ksemi*
     +                    8.854e-14*( Va )/(1.6022e-19*dopch) ) )/Ldi
        jtemp = joxide + int( (phim-phis)/dVa )
        write(6,'('' '',''   placing y nodes quadratically to '',
     +                  ''tepi...'')')
        do 700 j=0, (int( jchn*(phim-phit)/dVa )-1)
          Va = phim - (dVa/jchn) * j
  700     hy( j+jtemp ) = ( sqrt( K1 - 2.0*ksemi*8.854e-14*( Va-dVa )
     +                    /(1.6022e-19*dopch) ) - sqrt( K1 - 2.0*ksemi*
     +                    8.854e-14*( Va )/(1.6022e-19*dopch) ) )/Ldi
        jtemp = joxide + int( (phim-phis)/dVa ) + int((phim-phit)/dVa)
        write(6,'('' '',''   placing y nodes quadratically to '',
     +                  ''xd...'')')
        do 800 j=0, (int( phit/dVa )-1)
          Va = phit - dVa * j
  800     hy( j+jtemp ) = sqrt( 2.0*ksemi*8.854e-14/
     +                    (1.6022e-19*dopsub) )*( sqrt( Va ) - 
     +                    sqrt( Va-dVa ) )/Ldi
        jtemp = joxide + int( (phim-phis)/dVa ) + 
     +                   int( (phim-phit)/dVa ) + int( phit/dVa )
        write(6,'('' '',''   placing y nodes uniformly past xd '',
     +                  ''to tsub...'')')
        do 900 j=jtemp, ( ynodes-1 )
  900     hy( j ) = hy( jtemp-1 )
      endif
c
c     Now what kind of mesh do we need in the x direction?
c
c     Is it a full device simulation?
c
      if ( symtyp .eq. 'full' ) then
        write(6,'('' '')')
        write(6,'('' '',''   placing x nodes in the device...'')')
        do 1000 i=1,(xnodes-1)
          if ( ( i .le. int((lsorce+2.52058*tsorce)/hxnpl) ) .or.
     +         ( i .ge. (xnodes-int((ldrain+2.52058*tdrain)
     +                   /hxnpl)) ) ) then
            hx(i) = hxnpl/Ldi
          else
            hx(i) = hxchan/Ldi
          endif
 1000   continue
c
c     Is it a center channel simulation?
c
      elseif ( symtyp .eq. 'center' ) then
        write(6,'('' '')')
        write(6,'('' '',''   placing x nodes uniformly in the '',
     +                  ''channel...'')')
        do 1100 i=1,(xnodes-1)
c				this is a uniform mesh, after all...
 1100     hx(i) = (length/(xnodes-1))/Ldi
c
c     Is it a Mock simulation?
c
      elseif ( symtyp .eq. 'Mock' ) then
        write(6,'('' '')')
        write(6,'('' '',''   placing x nodes non-uniformly '',
     +                  ''in the channel...'')')
c
c       hx(xnodes/2) has already been set, so go from there
c
        do 1200 i=(xnodes/2), 2, -1
 1200     hx(i-1) = 0.9*hx(i)
        do 1210 i=(xnodes/2), (xnodes-2)
 1210     hx(i+1) = 0.9*hx(i)
      endif
c
c     Refine the mesh to what it really should be, if we did not
c     do the automatic mesh selection.
c
      if ( auto .ne. 'yes' ) then
        if ( debug .eq. 'yes' ) then
          write(6,'('' '')')
          write(6,'('' '',''This is the mesh BEFORE refinement:'')')
          call mkdev( .false. )
          call doguess( Vs, Vd, Vg1, eavg, .false., .false., error )
          if ( error ) return
          call depmap( 1, 1 )
        endif
        call refine( .true., .true., .true., error )
        if ( debug .eq. 'yes' ) then
          write(6,'('' '')')
          write(6,'('' '',''This is the mesh AFTER refinement:'')')
          call mkdev( .false. )
          call doguess( Vs, Vd, Vg1, eavg, .false., .false., error )
          if ( error ) return
          call depmap( 1, 1 )
        endif
        if ( error ) return
      endif
c
c     check the limits on the mesh as it stands now
c
      hxmax = -1.0
      hxmin = 1.0e10
      hxrat = -1.0
      hymax = -1.0
      hymin = 1.0e10
      hyrat = -1.0
c				check the limits for the hx()
      do 1300 i=1,(xnodes-1)
        hxmax = max( hxmax, hx(i) )
        hxmin = min( hxmin, hx(i) )
        if ( i .gt. 1 ) then
          if ( (hx(i)/hx(i-1)) .gt. hxrat ) hxrat = hx(i)/hx(i-1)
          if ( (hx(i-1)/hx(i)) .gt. hxrat ) hxrat = hx(i-1)/hx(i)
        endif
 1300 continue
c				check the limits for the hy()
      do 1400 j=1,(ynodes-1)
        hymax = max( hymax, hy(j) )
        hymin = min( hymin, hy(j) )
        if ( j .gt. 1 ) then
          if ( (hy(j)/hy(j-1)) .gt. hyrat ) hyrat = hy(j)/hy(j-1)
          if ( (hy(j-1)/hy(j)) .gt. hyrat ) hyrat = hy(j-1)/hy(j)
        endif
 1400 continue
c
c     Output some interesting statistics for this particular mesh
c
c     ...like the final results
c
      write(6,'('' '',''   setting xnodes = '',i3)') xnodes
      write(6,'('' '',''   setting ynodes = '',i3)') ynodes
      write(6,'('' '',''   maximum hx ratio = '',f6.3,'' from:''/
     +          '' '',''            hx max. = '',f6.4,'' (um)''/
     +          '' '',''            hx min. = '',f6.4,'' (um)'')')
     +                hxrat, hxmax*Ldi*1.0e4, hxmin*Ldi*1.0e4
      write(6,'('' '',''   maximum hy ratio = '',f6.3,'' from:''/
     +          '' '',''            hy max. = '',f6.4,'' (um)''/
     +          '' '',''            hy min. = '',f6.4,'' (um)'')')
     +                hyrat, hymax*Ldi*1.0e4, hymin*Ldi*1.0e4
c
c     return to caller
c
#     ifndef Cyber
        call flush(6)
#     endif
      error = .false.
      return
      end



      subroutine refine( init, reguess, verbose, error )
c
c     This subroutine refines the mesh based on the present values
c     in v(i,j) and the constant dVa. If any two nodes differ in
c     voltage by more than dVa then insert a node there and remake
c     the device and the initial guess at zero bias.
c
c     first the parameters for the sizes of the arrays
c
      integer Nx, Ny
      parameter( Nx=XMAX, Ny=YMAX )
c
c     Now for some parameters that are constants at the machine's
c     precision. Kc is the coupling safety factor for mesh refinement.
c
      PRECISION hmin, Kc
      parameter( hmin=1.0e-20, Kc=0.95 )
c
c     now for the arguments to this routine
c
      logical init, reguess, verbose, error
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
c
c     now for some local variables
c
      integer i, j, i1, j1, inew, jnew
      logical added
      PRECISION vsatr, vsatb, velr, velb, diffr, diffb,
     +          fr, fb, a, b, c, d, hl, hr, ht, hb, 
     +          vl, vr, vt, vb, vij,
     +          frx, fry, fbx, fby, axr, ayr, axb, ayb
c
c     start with a header so that the output looks reasonable
c
      if ( verbose ) then
        write(6,'('' '')')
        write(6,'('' '')')
        write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                             REB 1988'',
     +     ''    mesh refinement ===='')')
        write(6,'('' '')')
      endif
c
c     scan the lines i=1..xnodes, j=joxide+1..ynodes
c     for correct coupling criteria, and insert if needed.
c
  100 i = 1
      j = joxide + 1
      added = .false.
c
c     Now we need to refine the mesh based on what is made and 
c     the charge neutral guess.
c
  200 if ( init ) then
        call mkdev( .false. )
        call doguess( Vs, Vd, Vg1, eavg, .false., .false., error )
        if ( error ) return
      endif
c
c     get hl, hr, ht, hb (use symmetry at the boundaries)
c
  300 if ( i .eq. 1 ) then
        hl = hx(i)
      else
        hl = hx(i-1)
      endif

      if ( i .eq. xnodes ) then
        hr = hx(i-1)
      else
        hr = hx(i)
      endif

      if ( j .eq. 1 ) then
        ht = hy(j)
      else
        ht = hy(j-1)
      endif

      if ( j .eq. ynodes ) then
        hb = hy(j-1)
      else
        hb = hy(j)
      endif
c
c     get vr, vl, vt, vb, vij (using boundary condition vectors)
c
      if ( i .eq. 1 ) then
c
c       Check for LOS or fixed boundary condition
c
        if ( srccon(j,3) .eq. 1.0 ) then
          vl = srccon(j,1)
        else
          vl = v((i+1),j)
        endif
      else
        vl = v((i-1),j)
      endif

      if ( i .eq. xnodes ) then
c
c       Check for LOS or fixed boundary condition
c
        if ( drncon(j,3) .eq. 1.0 ) then
          vr = drncon(j,1)
        else
          vr = v((i-1),j)
        endif
      else
        vr = v((i+1),j)
      endif

      if ( j .eq. 1 ) then
c				only specify the gate voltage where
c				it is not there, otherwise it is the
c				same as where we are
        if ( vgate(i) .eq. 1.0e20 ) then
          vt = v(i,j)
        else
          vt = vgate(i)
        endif
      else
        vt = v(i,(j-1))
      endif

c
c     LOS at the bottom
c
      if ( j .eq. ynodes ) then
        vb = v(i,(j-1))
      else
        vb = v(i,(j+1))
      endif

      vij = v(i,j)
c
c     Now I need to compute the coefficients for the electric
c     field functions
c
      a = ( hl*vr - (hl+hr)*vij + hr*vl )/( hl*hr*(hl+hr) )
      b = ( ht*vb - (ht+hb)*vij + hb*vt )/( ht*hb*(ht+hb) )
      c = ( hl*hl*vr + (hr*hr - hl*hl)*vij - hr*hr*vl )/
     +         ( hl*hr*(hl+hr) )
      d = ( ht*ht*vb + (hb*hb - ht*ht)*vij - hb*hb*vt )/
     +         ( ht*hb*(ht+hb) )
c
c     Now let's compute the actual normalized electric field
c     components for both of the points of interest
c
      frx = -a*hr - c
      fry = -d
      fbx = -c
      fby = -b*hb - d
c
c     We will also need the unit vector decomposition for
c     later decomposing of the velocity
c
      fr = sqrt( frx*frx + fry*fry )
      if ( fr .ne. 0.0 ) then
        axr = frx / fr
        ayr = fry / fr
c				un-normalize the electric field for computing
        fr = fr * kT / Ldi
      else
        axr = (-hr - 1.0)/sqrt( (1.0 + hr)**2 + 1.0 )
        ayr = -1.0/sqrt( (1.0 + hr)**2 + 1.0 )
      endif

      fb = sqrt( fbx*fbx + fby*fby )
      if ( fb .ne. 0.0 ) then
        axb = fbx / fb
        ayb = fby / fb
c				un-normalize the electric field for computing
        fb = fb * kT / Ldi
      else
        axb = -1.0/sqrt( 1.0 + (1.0 + hb)**2 )
        ayb = (-hb - 1.0)/sqrt( 1.0 + (1.0 + hb)**2 )
      endif
c
c     Now compute the velocities and the diffusion coefficients
c
      vsatr = Av*exp(-fr/Bv) + Cv*(1.0/(1.0 + (fr/Dv)**Fv ))
      velr = mobil*fr*(1.0/sqrt(1.0+(mobil*fr/vsatr)**2)) * Ldi
      if ( fr .ne. 0.0 ) then
        diffr = Do + dDo*exp( -((log(fr)-Ad)/Bd)**2 )
      else
        diffr = Do
      endif
      velr = velr*axr

      vsatb = Av*exp(-fb/Bv) + Cv*(1.0/(1.0 + (fb/Dv)**Fv ))
      velb = mobil*fb*(1.0/sqrt(1.0+(mobil*fb/vsatb)**2)) * Ldi
      if ( fb .ne. 0.0 ) then
        diffb = Do + dDo*exp( -((log(fb)-Ad)/Bd)**2 )
      else
        diffr = Do
      endif
      velb = velb*ayb
c
c     See if we need to add a node in the x direction
c
c      if ( abs(hr*velr) .ge. Kc*(2.0*diffr) ) then
      if ( abs(vij-vr) .gt. (1.5*dVa/kT) ) then
c
c       OK... we need to add a node, let's see if
c       we can given the limitations imposed on the
c       system. For example, if this is the last
c       node in the row, we really need to be adding
c       a node "behind" us.
c
        if ( hr .le. (hmin/Ldi) ) then
c				hx already too small for us
          write(6,'('' '',''*** Error! ***''/
     +              '' '',''    Nodes need to be added, but the '',
     +                    ''lower limit on hx has been''/
     +              '' '',''    reached ('',g10.5,'' um). Shark '',
     +                    ''must halt due to ''/
     +              '' '',''coupling problems.'')') hmin*1.0e4
          error = .true.
          return
        else if ( xnodes .ge. Nx ) then
c				no more nodes to add in the x direction
          write(6,'('' '',''**** Error! ****''/
     +              '' '',''     Shark needs more x nodes than the '',
     +                    ''present program allows. Try a smaller''/
     +              '' '',''     channel length, or recompile '',
     +                    ''Shark with Nx increased.'')')
          error = .true.
          return
        else
c				finally... looks OK to add
          if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
            if ( i .lt. xnodes ) then
              write(6,'('' '',''Adding a node to the right of ('',
     +                        i3,'','',i3,'') because:'')') i, j
            else
              write(6,'('' '',''Adding a node to the left of ('',
     +                        i3,'','',i3,'') because:'')') i, j
            endif
            write(6,'('' '',''  Er='',g11.5,'' V/cm  vr='',g11.5,
     +                      '' cm/s  v/E='',f7.1,'' cm^2/V-s  hr='',
     +                      g8.2,'' um'')') fr, velr/Ldi, 
     +                      velr/(fr*Ldi), hr*Ldi*1.0e4
            write(6,'('' '',''  Dr='',f6.1,'' cm^2/sec  2D/v='',
     +                      g11.5,'' um  h*E='',f8.3,'' kT  dV='',
     +                      f8.3,'' kT'')') diffr, 
     +                      abs(2.0*diffr*Ldi/velr)*1.0e4,
     +                      hr*Ldi*fr/kT, abs(vij-vr)
          endif
c
c         We need to adjust the mesh. First, we need to be
c         sure that we can add the nodes, if not, then we'll
c         have to bail out saying "No can do, Boo".
c
c         OK... now we need to adjust the structure to add
c         a node in this location in the x-direction. This
c         means adjusting everything we can to the very best
c         of our abilities. This should be as small a 
c         'disturbing' force as possible.
c
c         Move everything over one in the x-direction to
c         make room for the new node...
c
          do 400 i1=xnodes,min((i+1), xnodes),-1
            if ( i1 .lt. xnodes ) hx( i1+1 ) = hx( i1 )
            vgate( i1+1 ) = vgate( i1 )
            Qss( i1+1 ) = Qss( i1 )
            do 400 j1=1,ynodes
              v( i1+1, j1 ) = v( i1, j1 )
              n( i1+1, j1 ) = n( i1, j1 )
              dop( i1+1, j1 ) = dop( i1, j1 )
  400     continue
c
c         Get the 'new' x node value so I can correct
c         all the new state values
c
          inew = min(i+1, xnodes)
c
c         Now I need to adjust the values in the 'new'
c         nodes for their 'new' position in the device
c
c				readjust the mesh (halve it)
          hx(inew - 1) = hx(inew - 1)/2.0
          hx(inew) = hx(inew - 1)
c				correct the other values
          vgate( inew ) = ( vgate(inew-1) + vgate(inew) )/2.0
          Qss( inew ) = sqrt( Qss(inew-1) * Qss(inew) )
c				readjust the new potential and carriers
          do 500 j1=1,ynodes
            v( inew, j1 ) = ( v(inew-1,j1) + v(inew,j1) )/2.0
            n( inew, j1 ) = sqrt( n(inew-1,j1) * n(inew,j1) )
            if ( dop(inew-1,j1) * dop(inew,j1) .gt. 0.0 ) then
              dop( inew, j1 ) = sign( sqrt( dop(inew-1,j1) * 
     +                               dop(inew,j1) ), dop(inew,j1) )
            else
              dop( inew, j1 ) = exp( ( log(max(dop(inew-1,j1), 
     +                     dop(inew,j1))) - log(-min(dop(inew-1,j1), 
     +                     dop(inew,j1))) )/2.0 )
            endif
  500     continue
c				update the size, and tag it as added
          xnodes = xnodes + 1
          added = .true.
c
c         At this point, it is still possible that the
c         change is too big to keep, and so rather than
c         continue with this error, let's go back to the
c         top and re-process this node and see if there
c         are more changes. The only exception to this
c         is when the node inserted was at the very far
c         right side. In that case, we really need to
c         stop the mesh refinement and let the system
c         readjust. If not, we'll be adding nodes in a
c         loop without changing anything.
c
          if ( inew .eq. (xnodes-1) ) then
            if ( (debug .ne. 'no') .and. (added) ) then
c				write a special message
              write(6,'('' '',''A node was added at the right '',
     +                        ''side, so now xnodes='',i3,
     +                        '' and ynodes='',i3,''.'')') 
     +              xnodes, ynodes
            endif
            goto 800
          else
c				not at the end, so changes are OK to device
            goto 200
          endif
        endif
      endif
c
c     See if we need to add a node in the y direction
c
c      if ( abs(hb*velb) .ge. Kc*(2.0*diffb) ) then
      if ( abs(vij-vb) .gt. (1.5*dVa/kT) ) then
c
c       OK... we need to add a node, let's see if
c       we can given the limitations imposed on the
c       system.
c
        if ( hb .le. (hmin/Ldi) ) then
c				hy already too small for us
          write(6,'('' '',''*** Error! ***''/
     +              '' '',''    Nodes need to be added, but the '',
     +                    ''lower limit on hy has been''/
     +              '' '',''    reached. Shark must halt due to '',
     +                    ''coupling problems.'')')
          error = .true.
          return
        else if ( ynodes .ge. Ny ) then
c				no more nodes to add in the y direction
          write(6,'('' '',''**** Error! ****''/
     +              '' '',''     Shark needs more y nodes than the '',
     +                    ''present program allows. Try a smaller''/
     +              '' '',''     channel length, or recompile '',
     +                    ''Shark with Ny increased.'')')
          error = .true.
          return
        else
c				finally... looks OK to add
          if ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) then
            write(6,'('' '',''Adding a node below ('',
     +                      i3,'','',i3,'') because:'')') i, j
            write(6,'('' '',''  Eb='',g11.5,'' V/cm  vb='',g11.5,
     +                      '' cm/s  v/E='',f7.1,'' cm^2/V-s  hb='',
     +                      g8.2,'' um'')') fb, velb/Ldi, 
     +                      velb/(fb*Ldi), hb*Ldi*1.0e4
            write(6,'('' '',''  Db='',f6.1,'' cm^2/sec  2D/v='',
     +                      g11.5,'' um  h*E='',f8.3,'' kT  dV='',
     +                      f8.3,'' kT'')') diffb, 
     +                      abs(2.0*diffb*Ldi/velb)*1.0e4,
     +                      hb*Ldi*fb/kT, abs(vij-vb)
          endif
c
c         We need to adjust the mesh. First, we need to be
c         sure that we can add the nodes, if not, then we'll
c         have to bail out saying "No can do, Boo".
c
c         OK... now we need to adjust the structure to add
c         a node in this location in the y-direction. This
c         means adjusting everything we can to the very best
c         of our abilities. This should be as small a 
c         'disturbing' force as possible.
c
c         Move everything down one in the y-direction to
c         make room for the new node...
c
          do 700 j1=ynodes,min((j+1), ynodes),-1
            if ( j1 .lt. ynodes ) hy( j1+1 ) = hy( j1 )
            do 600 i1=1,3
              srccon( j1+1, i1 ) = srccon( j1, i1 )
              drncon( j1+1, i1 ) = drncon( j1, i1 )
  600       continue
            do 700 i1=1,xnodes
              v( i1, j1+1 ) = v( i1, j1 )
              n( i1, j1+1 ) = n( i1, j1 )
              dop( i1, j1+1 ) = dop( i1, j1 )
  700     continue
c
c         Get the 'new' y node value so I can correct
c         it's values
c
          jnew = min(j+1, ynodes)
c
c         Now I need to adjust the values in the 'new'
c         nodes for their 'new' position in the device
c
c				readjust the mesh (halve it)
          hy(jnew-1) = hy(jnew-1)/2.0
          hy(jnew) = hy(jnew-1)
c				correct the other values
c         do 800 i1=1,3
c           srccon( jnew, i1 ) = srccon( jnew, i1 )
c           drncon( jnew, i1 ) = drncon( jnew, i1 )
c 800     continue
c				readjust the new potential and carriers
          do 900 i1=1,xnodes
            v( i1, jnew ) = ( v(i1,jnew-1) + v(i1,jnew) )/2.0
            n( i1, jnew ) = sqrt( n(i1,jnew-1) * n(i1,jnew) )
            if ( dop(i1,jnew-1) * dop(i1,jnew) .gt. 0.0 ) then
              dop( i1, jnew ) = sign( sqrt( dop(i1,jnew-1) * 
     +                               dop(i1,jnew) ), dop(i1,jnew) )
            else
              dop( i1, jnew ) = exp( ( log(max(dop(i1,jnew-1), 
     +                     dop(i1,jnew))) - log(-min(dop(i1,jnew-1), 
     +                     dop(i1,jnew))) )/2.0 )
            endif
  900     continue
c				update the size, and tag it as added
          if ( j .lt. joxide ) then
            joxide = joxide + 1
          endif
          ynodes = ynodes + 1
          added = .true.
c
c         At this point, it is still possible that the
c         change is too big to keep, and so rather than
c         continue with this error, let's go back to the
c         top and re-process this node and see if there
c         are more changes.
c
          goto 200
        endif
      endif
c
c     Check to see if we're at the end of this row
c
      if ( i .ge. xnodes ) then
c
c       See if we are doing this from scratch
c
        if ( init ) then
c
c         See if any were added in this pass
c
          if ( added ) then
c
c           If we are doing this from scratch, go back and re-guess
c           the device and try again
c
            if ( debug .ne. 'no' ) then
              write(6,'('' '',''Nodes were added, so now xnodes='',
     +                        i3,'' and ynodes='',i3,''.'')') 
     +                  xnodes, ynodes
            endif
            goto 100
          else
c
c           If initially, we had no additions, then we must be done due
c           to the direction of the sweeping, i.e. towards lower E field.
c
            goto 800
          endif
        else
c
c         If not new, keep checking the entire device for all the 
c         changed rows.
c
          if ( j .lt. ynodes ) then
c				go to the start of the next row
            i = 1
            j = j + 1
            goto 200
          endif
          if ( (debug .ne. 'no') .and. (added) ) then
c				write a message and reset the changed flag
            write(6,'('' '',''Nodes were added, so now xnodes='',
     +                      i3,'' and ynodes='',i3,''.'')') 
     +            xnodes, ynodes
            added = .false.
          endif
c
c         We had a given problem, so re-guess the system and return
c         to the calling routine - if they want it...
c
          if ( reguess ) then
            call mkdev( .false. )
            call doguess( Vs, Vd, Vg1, eavg, .false., .false., error )
          endif
          goto 800
        endif
      else
c				go to the next node in the row and continue
        i = i + 1
        goto 300
      endif
c
c     now return to the caller
c
  800 error = .false.
      return
      end
