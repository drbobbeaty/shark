      subroutine mkdev( verbose )
c
c     This subroutine really does the bulk of the set-up for the
c     simulation. It sets the doping density as a function of
c     the inputted constants and the resulting mesh, as well as
c     setting the dielectric constant.
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
      logical verbose
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
      integer i, j, j1
      PRECISION xtemp, ytemp, ysrc
c
c     start with a header so that the output looks reasonable
c
      if ( verbose ) then
        write(6,'('' '')')
        write(6,'('' '')')
        write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                            REB 1988'',
     +     ''  setting up device ===='')')
        write(6,'('' '')')
      endif
c
c     start by resetting the temporary variables, and the arrays
c
      do 100 j=1,ynodes
        do 100 i=1,xnodes
          dop(i,j) = 0.0
  100 continue
      j = 1
      ytemp = hy(1)*Ldi
c
c     start at the top and work our way down, first the oxide
c
      if ( verbose ) write(6,'('' '',''  setting up the oxide...'')')
  200 do 300 i=1,xnodes
  300   dop(i,j) = 0.0
      ytemp = ytemp + hy(j) * Ldi
      j = j + 1
      if ( j .le. joxide ) goto 200
c
c     ...then the channel
c
      if ( verbose ) write(6,'('' '',''  setting up the channel...'')')
  400 if ( chmod .eq. 'buried' ) then
        if ( (doptyp .eq. 'step') .or.
     +       (doptyp .eq. 'tailed') .and. (ytemp .le. (tox+tepi)) ) then
          do 500 i=1,xnodes
  500       dop(i,j) = dopch/nio
        elseif ( doptyp .eq. 'tailed' ) then
c				this is the tailing part

        endif
      else
c				this is the MODFET channel
        do 600 i=1,xnodes
  600     dop(i,j) = -dopsub/nio
      endif
      ytemp = ytemp + hy(j) * Ldi
      j = j + 1
      if ( ( (chmod .ne. 'buried') .and. (ytemp .le. (tox+tepi)) ) .or.
     +     ( (chmod .eq. 'buried') .and. (doptyp .ne. 'tailed') .and.
     +       (ytemp .le. (tox+tepi)) ) .or.
     +     ( (chmod .eq. 'buried') .and. (doptyp .eq. 'tailed') .and.
     +       (ytemp .le. (tox+tepi*3.52058)) ) ) goto 400
c
c     ...then all the way out to the back of the substrate
c
      if ( verbose ) write(6,'('' '',''  setting the substrate...'')')
      do 700 j1=j,ynodes
        do 700 i=1,xnodes
  700     dop(i,j1) = -dopsub/nio
c
c     See if we are simulating just the center section of the FET,
c     and if so then skip putting in the drain and source
c
      if ( symtyp .ne. 'full' ) goto 1100
c
c     Start by offsetting the values of doping by dopsub so that
c     we will get smooth transisitions between the N+ and P regions
c
      dopsrc = dopsrc + dopsub
      dopdrn = dopdrn + dopsub
c
c     Now go back and do the source scanning in the x direction
c
      if ( verbose ) then
        if ( srccor .eq. 'square' ) then
          write(6,'('' '',''  setting the source with a rectangular '',
     +                    ''shape...'')')
        elseif ( srccor .eq. 'sloped' ) then
          write(6,'('' '',''  setting the source with a sloped '',
     +                    ''shape...'')')
        elseif ( srccor .eq. 'curved' ) then
          write(6,'('' '',''  setting the source with a curved '',
     +                    ''shape...'')')
        endif
      endif
c
c     reset the variables
c
      j = joxide + 1
      ytemp = ( hy(1)*Ldi )*( joxide+1 )
      ysrc = ytemp
      i = 1
      xtemp = 0
c
c     ...do the source
c
c     ...for the square (idealized) source
c
  900 if ( srccor .eq. 'square' ) then
        if ( (xtemp .le. lsorce) .and. (ytemp .le. (tox+tsorce)) ) then
          dop(i,j) = dopsrc/nio
        elseif ( (xtemp .le. lsorce) .and. 
     +           (ytemp .gt. (tox+tsorce)) ) then
          dop(i,j) = (dopsrc/nio)*exp(-(ytemp-(tox+tsorce))**2/
     +                   (0.352963539*tsorce*tsorce) )
        elseif ( (xtemp .gt. lsorce) .and. 
     +           (ytemp .le. (tox+tsorce)) ) then
          dop(i,j) = (dopsrc/nio)*exp(-(xtemp-lsorce)**2/
     +                   (0.352963539*tsorce*tsorce) )
        else
          dop(i,j) = (dopsrc/nio)*exp(-( (ytemp-(tox+tsorce))
     +                   + (xtemp-lsorce) )**2/
     +                   (0.352963539*tsorce*tsorce) )
        endif
c
c     ...for the sloped (ion implanted) source
c
      elseif ( srccor .eq. 'sloped' ) then
        if ( (xtemp .le. (lsorce-(ytemp-ysrc))) .and. 
     +       (ytemp .le. (tox+tsorce)) ) then
          dop(i,j) = dopsrc/nio
        elseif ( (xtemp .le. (lsorce-(ytemp-ysrc))) .and. 
     +           (ytemp .gt. (tox+tsorce)) ) then
          dop(i,j) = (dopsrc/nio)*exp(-(ytemp-(tox+tsorce))**2/
     +                   (0.352963539*tsorce*tsorce) )
        elseif ( (xtemp .gt. (lsorce-(ytemp-ysrc))) .and. 
     +           (ytemp .le. (tox+tsorce)) ) then
          dop(i,j) = (dopsrc/nio)*exp(-(xtemp-(lsorce-
     +                   (ytemp-ysrc)))**2/
     +                   (0.352963539*tsorce*tsorce) )
        else
          dop(i,j) = (dopsrc/nio)*exp(-( (ytemp-(tox+tsorce))
     +                   + (xtemp-(lsorce-(ytemp-ysrc))) )**2/
     +                   (0.352963539*tsorce*tsorce) )
        endif
c
c     ...for the curved (diffused) source
c
      elseif ( srccor .eq. 'curved' ) then
        if ( (xtemp .le. (lsorce-tsorce+sqrt( abs(tsorce**2 -
     +                    (ytemp - ysrc) )**2)) ) .and. 
     +       (ytemp .le. (tox+tsorce)) ) then
          dop(i,j) = dopsrc/nio
        elseif ( (xtemp .le. (lsorce-tsorce)) .and. 
     +           (ytemp .gt. (tox+tsorce)) ) then
          dop(i,j) = (dopsrc/nio)*exp(-(ytemp-(tox+tsorce))**2/
     +                   (0.352963539*tsorce*tsorce) )
        elseif ( (xtemp .gt. (lsorce-tsorce+sqrt( abs(tsorce**2 -
     +                        (ytemp-ysrc) )**2)) ) .and. 
     +           (ytemp .le. (tox+tsorce)) ) then
          dop(i,j) = (dopsrc/nio)*exp(-(xtemp-(lsorce-tsorce+
     +                   sqrt(tsorce**2-(ytemp-ysrc)**2)))**2/
     +                   (0.352963539*tsorce*tsorce) )
        else
          dop(i,j) = (dopsrc/nio)*exp(-(xtemp-(lsorce-tsorce)
     +          + (ytemp-(tox+tsorce)) )**2/
     +            (0.352963539*tsorce*tsorce) )
        endif
      endif
c
c     if it is less than 10**5 then make it zero
c
      if ( dop(i,j) .lt. (1.0e5/nio) ) dop(i,j) = 0.0
c
c     Now we have to remove the offset put into dop(i,j)
c
      dop(i,j) = dop(i,j) - (dopsub/nio)
c				scan out the x direction
      xtemp = xtemp + hx(i)*Ldi
      i = i + 1
      if ( ( (xtemp .le. lsorce+2.52058*tsorce) .and.
     +                            (srccor .eq. 'square') ) .or.
     +     ( (xtemp .le. (lsorce-(ytemp-ysrc)+2.52058*tsorce)) .and.
     +                            (srccor .eq. 'sloped') ) .or.
     +     ( (xtemp .le. (lsorce-tsorce+sqrt((3.52058*tsorce)**2-
     +                                   (ytemp-ysrc)**2)) ) .and.
     +                            (srccor .eq. 'curved') ) )
     +    goto 900
      i = 1
      xtemp = 0
c				scan down the y direction
      ytemp = ytemp + hy(j)*Ldi
      j = j + 1
      if ( ytemp .le. (tox+3.52058*tsorce) ) goto 900

c
c     Now go back and do the drain scanning in the -x direction
c
      if ( verbose ) then
        if ( drncor .eq. 'square' ) then
          write(6,'('' '',''  setting the drain with a rectangular '',
     +                    ''shape...'')')
        elseif ( drncor .eq. 'sloped' ) then
          write(6,'('' '',''  setting the drain with a sloped '',
     +                    ''shape...'')')
        elseif ( drncor .eq. 'curved' ) then
          write(6,'('' '',''  setting the drain with a curved '',
     +                    ''shape...'')')
        endif
      endif
c
c     reset the variables
c
      j = joxide + 1
      ytemp = ( hy(1)*Ldi )*( joxide+1 )
      i = xnodes
      xtemp = 0
c
c     ...do the drain
c
c     ...for the square (idealized) drain
c
 1000 if ( drncor .eq. 'square' ) then
        if ( (xtemp .le. ldrain) .and. 
     +       (ytemp .le. (tox+tsorce)) ) then
          dop(i,j) = dopdrn/nio
        elseif ( (xtemp .le. ldrain) .and. 
     +           (ytemp .gt. (tox+tdrain)) ) then
          dop(i,j) = (dopdrn/nio)*exp(-(ytemp-(tox+tdrain))**2/
     +                   (0.352963539*tdrain*tdrain) )
        elseif ( (xtemp .gt. ldrain) .and. 
     +           (ytemp .le. (tox+tdrain)) ) then
          dop(i,j) = (dopdrn/nio)*exp(-(xtemp-ldrain)**2/
     +                   (0.352963539*tdrain*tdrain) )
        else
          dop(i,j) = (dopdrn/nio)*exp(-( (ytemp-(tox+tdrain))
     +                   + (xtemp-ldrain) )**2/
     +                   (0.352963539*tdrain*tdrain) )
        endif
c
c     ...for the sloped (ion implanted) drain
c
      elseif ( drncor .eq. 'sloped' ) then
        if ( (xtemp .le. (ldrain-(ytemp-ysrc))) .and. 
     +       (ytemp .le. (tox+tdrain)) ) then
          dop(i,j) = dopdrn/nio
        elseif ( (xtemp .le. (ldrain-(ytemp-ysrc))) .and. 
     +           (ytemp .gt. (tox+tdrain)) ) then
          dop(i,j) = (dopdrn/nio)*exp(-(ytemp-(tox+tdrain))**2/
     +                   (0.352963539*tdrain*tdrain) )
        elseif ( (xtemp .gt. (ldrain-(ytemp-ysrc))) .and. 
     +           (ytemp .le. (tox+tdrain)) ) then
          dop(i,j) = (dopdrn/nio)*exp(-(xtemp-(ldrain-
     +                   (ytemp-ysrc)))**2/
     +                   (0.352963539*tdrain*tdrain) )
        else
          dop(i,j) = (dopdrn/nio)*exp(-( (ytemp-(tox+tdrain))
     +                   + (xtemp-(ldrain-(ytemp-ysrc))) )**2/
     +                   (0.352963539*tdrain*tdrain) )
        endif
c
c     ...for the curved (diffused) drain
c
      elseif ( drncor .eq. 'curved' ) then
        if ( (xtemp .le. (ldrain-tdrain+sqrt( abs(tdrain**2 -
     +                    (ytemp - ysrc) )**2)) ) .and. 
     +       (ytemp .le. (tox+tdrain)) ) then
          dop(i,j) = dopdrn/nio
        elseif ( (xtemp .le. (ldrain-tdrain)) .and. 
     +           (ytemp .gt. (tox+tdrain)) ) then
          dop(i,j) = (dopdrn/nio)*exp(-(ytemp-(tox+tdrain))**2/
     +                   (0.352963539*tdrain*tdrain) )
        elseif ( (xtemp .gt. (ldrain-tdrain+sqrt( abs(tdrain**2 -
     +                        (ytemp-ysrc) )**2)) ) .and. 
     +           (ytemp .le. (tox+tdrain)) ) then
          dop(i,j) = (dopdrn/nio)*exp(-(xtemp-(ldrain-tdrain+
     +                   sqrt(tdrain**2-(ytemp-ysrc)**2)))**2/
     +                   (0.352963539*tdrain*tdrain) )
        else
          dop(i,j) = (dopdrn/nio)*exp(-(xtemp-(ldrain-tdrain)
     +          + (ytemp-(tox+tdrain)) )**2/
     +            (0.352963539*tdrain*tdrain) )
        endif
      endif
c
c     if it is less than 10**5 then make it zero
c
      if ( dop(i,j) .lt. (1.0e5/nio) ) dop(i,j) = 0.0
c
c     Now we have to remove the offset put into dop(i,j)
c
      dop(i,j) = dop(i,j) - (dopsub/nio)
c				scan in the x direction
      xtemp = xtemp + hx(i)*Ldi
      i = i - 1
      if ( ( (xtemp .le. ldrain+2.52058*tdrain) .and.
     +                            (drncor .eq. 'square') ) .or.
     +     ( (xtemp .le. (ldrain-(ytemp-ysrc)+2.52058*tdrain)) .and.
     +                            (drncor .eq. 'sloped') ) .or.
     +     ( (xtemp .le. (ldrain-tdrain+sqrt((3.52058*tdrain)**2-
     +                                   (ytemp-ysrc)**2)) ) .and.
     +                            (drncor .eq. 'curved') ) )
     +    goto 1000
      i = xnodes
      xtemp = 0
c				scan down the y direction
      ytemp = ytemp + hy(j)*Ldi
      j = j + 1
      if ( ytemp .le. (tox+3.52058*tdrain) ) goto 1000
c
c     Now that we are done, we have to remove the offset put in
c     initially.
c
      dopsrc = dopsrc - dopsub
      dopdrn = dopdrn - dopsub
c
c     return to the caller
c
#     ifndef Cyber
 1100   if ( verbose ) then
          call flush(6)
        endif
        return
#     else
 1100   return
#     endif
      end
