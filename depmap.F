      subroutine depmap( mode, early )
c
c     This routine prints out a very nice little depletion map
c     for the device under simulation. If mode=1 then we will
c     print out a very fine discretization of the doping levels,
c     if mode=0 the just N and P will be printed out.
c
c     If early=1 then we allow the routine to stop printing early
c     if we are well into the substrate.
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
      integer mode, early
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
      integer skipx, skipy, i, i1, j, prtsub, vi, vj, l
      PRECISION src, sub, Vlo, Vhi, nmin, nmax, dVmax
      character*1 dir
      character*132 conout, subout
c
c     get the minimum and maximum variables for the system
c
      Vhi = -1.0e10
      nmax = -1.0e10
      Vlo = 1.0e10
      nmin = 1.0e10
      dVmax = -1.0
      do 50 j=joxide,ynodes
        do 50 i=1,xnodes
          Vlo = min( Vlo, v(i,j) )
          Vhi = max( Vhi, v(i,j) )
          if ( j .gt. joxide ) then
            nmin = min( nmin, n(i,j) )
            nmax = max( nmax, n(i,j) )
            if ( (i .gt. 1) .and. 
     +           (abs(v(i,j)-v(i-1,j)) .gt. dVmax) ) then
              dVmax = abs(v(i,j)-v(i-1,j))
              vi = i
              vj = j
              dir = '<'
            endif
            if ( (i .lt. xnodes) .and. 
     +           (abs(v(i,j)-v(i+1,j)) .gt. dVmax) ) then
              dVmax = abs(v(i,j)-v(i+1,j))
              vi = i
              vj = j
              dir = '>'
            endif
            if ( (j .gt. (joxide+1)) .and. 
     +           (abs(v(i,j)-v(i,j-1)) .gt. dVmax) ) then
              dVmax = abs(v(i,j)-v(i,j-1))
              vi = i
              vj = j
              dir = '^'
            endif
            if ( (j .lt. ynodes) .and.
     +           (abs(v(i,j)-v(i,j+1)) .gt. dVmax) ) then
              dVmax = abs(v(i,j)-v(i,j+1))
              vi = i
              vj = j
              dir = 'v'
            endif
          endif
   50 continue
c
c     calculate the skip factor for the x direction
c
      write(6,'('' '')')
      skipx = int( (xnodes/115.0) + 0.99999 )
      skipy = int( (ynodes/50.0) + 0.99999 )
      write(6,'('' '','' +------> '',i1,'' x nodes/character'',
     + ''         maximum voltage = '',f7.2,'' kT/q = '',f7.3,
     + '' V'')') skipx, Vhi, Vhi*kT
      write(6,'('' '','' !                    at '',i3,'','',i3,
     + '' going '',a1,'' the local dV = '',f7.2,'' kT/q = '',f7.3,
     + '' V'')') vi, vj, dir, dVmax, dVmax*kT
      write(6,'('' '','' !                           '',
     + ''             elec. conc. = '',g8.2,'' nio = '',g8.2,
     + '' /cm**3'')') nmax, nmax*nio
      write(6,'('' '','' ! '',i1,'' y nodes/character       '',
     + ''         minimum voltage = '',f7.2,'' kT/q = '',f7.3,
     + '' V'')') skipy, Vlo, Vlo*kT
      write(6,'('' '','' !                           '',
     + ''             elec. conc. = '',g8.2,'' nio = '',g8.2,
     + '' /cm**3'')') nmin, nmin*nio
      write(6,'('' '','' V '')')
      write(6,'('' '')')
c
c     do one line of output for each node line in the y direction
c
c				reset the index into the string
      i1 = 1
c				reset the output strings
      conout = ' '
      subout = ' '
c				add the leading space for the contact
      conout(i1:i1) = ' '
      subout(i1:i1) = ' '
      i1 = i1 + 1
c
c     now scan all the way across the device
c
      do 100 i=1,xnodes,skipx
c				decode the constant into a character
        if ( vgate(i) .lt. 1.0e20 ) then
          conout(i1:i1) = 'g'
        else
          conout(i1:i1) = ' '
        endif
c				set subout to be all P's as a test
#       ifdef Cyber
          if ( i1 .eq. i1 ) then
            subout(i1:i1) = 'p'
          endif
#       else
          subout(i1:i1) = 'p'
#       endif
        i1 = i1 + 1
  100 continue
c				add the trailing space for the contact
#     ifdef Cyber
        if ( i1 .eq. i1 ) then
          conout(i1:i1) = ' '
          subout(i1:i1) = ' '
        endif
#     else
        conout(i1:i1) = ' '
        subout(i1:i1) = ' '
#     endif
c				now output it
      if ( i1 .lt. 80 ) then
        write(6,*) conout(:i1)
      else
        write(6,'('' '',117a1)') (conout(l:l),l=1,i1)
      endif
c
c     reset the number of times we have printed out the substrate
c
      prtsub = 0
c
c     now do the bulk of the simulator
c
      do 300 j=1,ynodes,skipy
c				reset the index into the string
        i1 = 1
c				reset the output string
        conout = ' '
c
c       make the first character an s if the source contact is there
c
        if ( (srccon(j,3) .eq. 1.0) .and.
     +       (symtyp .ne. 'center') ) then
          conout(i1:i1) = 's'
        endif
c				up the index count
        i1 = i1 + 1
c
c       now scan across the interesting region
c
        do 200 i=1,xnodes,skipx
c
c         decode the voltage into a character, depending on the region
c         we are in
c
c         check the n-type regions first
c
          if ( dop(i,j) .gt. 0.0 ) then
c            src = exp( v(i,j) - log(dop(i,j)) )*10.0
c            src = exp( v(i,j) - log(min(dopsrc,dopdrn)/nio) )*10.0
            src = ( n(i,j)/dop(i,j) )*10.0
            if ( mode .eq. 1 ) then
              if (src .gt. 1.0 ) then
                conout(i1:i1) = '1'
              endif
              if (src .gt. 2.0 ) then
                conout(i1:i1) = '2'
              endif
              if (src .gt. 3.0 ) then
                conout(i1:i1) = '3'
              endif
              if (src .gt. 4.0 ) then
                conout(i1:i1) = '4'
              endif
              if (src .gt. 5.0 ) then
                conout(i1:i1) = '5'
              endif
              if (src .gt. 6.0 ) then
                conout(i1:i1) = '6'
              endif
              if (src .gt. 7.0 ) then
                conout(i1:i1) = '7'
              endif
              if (src .gt. 8.0 ) then
                conout(i1:i1) = '8'
              endif
              if (src .gt. 9.0 ) then
                conout(i1:i1) = '9'
              endif
              if (src .gt. 9.5 ) then
                conout(i1:i1) = 'n'
              endif
              if (src .gt. 11.0 ) then
                conout(i1:i1) = 'N'
              endif
            else
              if ( src .gt. 9.0 ) then
                conout(i1:i1) = 'n'
              endif
            endif
          else
c
c           check the p-type next, being careful about near zero
c           dopings.
c
            if ( dop(i,j) .lt. 0.0 ) then
c              sub = exp( -v(i,j) - log(-dop(i,j)) )*10.0
              if ( symtyp .eq. 'Mock' ) then
                sub = exp( -v(i,j) - log(dopsub/nio) - 
     +                     log(dopsrc/nio) )*10.0
              else
                sub = exp( -v(i,j) - log(dopsub/nio) )*10.0
              endif
              if ( mode .eq. 1 ) then
                if ( dop(i,j-1) .gt. 0.0 ) then
                  conout(i1:i1) = '^'
                endif
c                if ( sub .gt. 1.0 ) then
c                  conout(i1:i1) = '1'
c                endif
c                if ( sub .gt. 2.0 ) then
c                  conout(i1:i1) = '2'
c                endif
c                if ( sub .gt. 3.0 ) then
c                  conout(i1:i1) = '3'
c                endif
c                if ( sub .gt. 4.0 ) then
c                  conout(i1:i1) = '4'
c                endif
c                if ( sub .gt. 5.0 ) then
c                  conout(i1:i1) = '5'
c                endif
c                if ( sub .gt. 6.0 ) then
c                  conout(i1:i1) = '6'
c                endif
c                if ( sub .gt. 7.0 ) then
c                  conout(i1:i1) = '7'
c                endif
                if ( sub .gt. 8.0 ) then
                  conout(i1:i1) = '8'
                endif
                if ( sub .gt. 9.0 ) then
                  conout(i1:i1) = '9'
                endif
                if ( sub .gt. 9.5 ) then
                  conout(i1:i1) = 'p'
                endif
                if ( sub .ge. 10.0 ) then
                  conout(i1:i1) = 'P'
                endif
c				check for electrons
c                if ( n(i,j) .gt. (1.0e-5*dopsub/nio) ) then
c                  conout(i1:i1) = '.'
c                endif
                if ( n(i,j) .gt. (0.1*dopsub/nio) ) then
                  conout(i1:i1) = '-'
                endif
                if ( n(i,j) .gt. (0.9*dopsub/nio) ) then
                  conout(i1:i1) = 'n'
                endif
                if ( n(i,j) .gt. (dopsub/nio) ) then
                  conout(i1:i1) = 'N'
                endif
              else
                if ( dop(i,j-1) .gt. 0.0 ) then
                  conout(i1:i1) = '^'
                endif
c                if ( sub .gt. 9.0 ) then
c                  conout(i1:i1) = 'p'
c                endif
                if ( sub .gt. 9.5 ) then
                  conout(i1:i1) = 'P'
                endif
c				check for inversion
c                if ( n(i,j) .gt. (1.0e-5*dopsub/nio) ) then
c                  conout(i1:i1) = '.'
c                endif
                if ( n(i,j) .gt. (0.1*dopsub/nio) ) then
                  conout(i1:i1) = '-'
                endif
                if ( n(i,j) .gt. (0.9*dopsub/nio) ) then
                  conout(i1:i1) = 'n'
                endif
                if ( n(i,j) .gt. (dopsub/nio) ) then
                  conout(i1:i1) = 'N'
                endif
              endif
            else
c
c             check for the oxide last
c
              if ( j .le. joxide ) then
                conout(i1:i1) = 'i'
              endif
            endif
          endif
c
c         go to the next position in the output string
c
          i1 = i1 + 1
  200   continue
c
c       make the last character a d if the drain contact is there
c
        if ( (drncon(j,3) .eq. 1.0) .and. (symtyp .ne. 'center') ) then
          conout(i1:i1) = 'd'
        endif
c				now output it
        if ( i1 .lt. 80 ) then
          write(6,*) conout(:i1)
        else
          write(6,'('' '',117a1)') (conout(l:l),l=1,i1)
        endif
c
c       see if this is just the substrate, and if so tally another up.
c       If we get to 5, then we can see about leaving early
c
        if ( conout .eq. subout ) then
          prtsub = prtsub + 1
          if ( (prtsub .ge. 5) .and. (early .eq. 1) ) goto 400
        else
          prtsub = 0
        endif
  300 continue
c
c     ...add the final few lines
c
  400 write(6,'('' '')')
      write(6,'('' '')')
      write(6,'('' '')')
#     ifndef Cyber
        call flush(6)
#     endif
c
c     and return to the caller
c
      return
      end
