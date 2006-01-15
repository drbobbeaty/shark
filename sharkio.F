      subroutine dump( tnow, plot, xref )
c
c     This routine dumps out the status of the simulator into two
c     different files. "shark.core" has all the common block variables
c     in a form suitable for reading back in. "shark.out" has the
c     plottable stuff. If 'plot' is .true. then dump the plot stuff.
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
      logical plot
      PRECISION tnow, xref
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
c     now some local variables for the dumping
c
      integer i, j, k
      PRECISION xtemp, ntemp, dely, Vob
c
c     first, the simplest, dump the contents of each common block
c     that is of any use to me.
c
c     open the "core dump" file
c
#     ifndef Cyber
      open( 10, file="shark.core", status="UNKNOWN" )
#     endif
      rewind( 10 )
c
c     first, dump /device/, which is:
c
c     common /device/ n(Nx,Ny), dop(Nx,Ny), v(Nx,Ny), vgate(Nx),
c    +                drncon(Ny,3), srccon(Ny,3), joxide
c     integer joxide
c     PRECISION n, dop, v, vgate, drncon, srccon
c
      do 100 i=1,Nx
        do 100 j=1,Ny
  100     write(10,'('' '',3(e19.10E3))') n(i,j), dop(i,j), v(i,j)
      do 200 i=1,Nx
  200   write(10,'('' '',e19.10E3)') vgate(i)
      do 300 j=1,Ny
        write(10,'('' '',3(e19.10E3))') (srccon(j,k),k=1,3)
  300   write(10,'('' '',3(e19.10E3))') (drncon(j,k),k=1,3)
      write(10,'('' '',i3)') joxide
c
c     next dump /mesh/, which is:
c
c     common /mesh/ hx(Nx), hy(Ny), xnodes, ynodes, xfrac, auto, dVa
c     integer xnodes, ynodes
c     PRECISION hx, hy, xfrac, dVa
c     character*8 auto
c
      do 400 j=1,Nx
  400   write(10,'('' '',e19.10E3)') hx(j)
      do 500 j=1,Ny
  500   write(10,'('' '',e19.10E3)') hy(j)
      write(10,'('' '',i3,'' '',i3,'' '',f6.4,'' '',e19.10E3)')
     +          xnodes, ynodes, xfrac, dVa
      write(10,'(1x,a)') auto
c
c     next dump /params/, which is:
c
c     common /params/ kT, nio, Ldi, kox, chmod, doptyp, ksemi,
c    +                dt, tstart, tstop, stpmin, tout, delt,
c    +                Vs, Vd, vg1, eavg, soltyp, symtyp,
c    +                icchan, chtyp, instyp, vmax, debug,
c    +                guess, delV, delE, maxit, keep
c     integer maxit, keep
c     PRECISION kT, nio, Ldi, kox, ksemi, dt, tstart, tstop,
c    +     stpmin, tout, Vs, Vd, vg1, eavg, vmax, delt,
c    +     delV, delE
c     character*8 chmod, doptyp, icchan, chtyp, instyp, debug, guess,
c    +     soltyp, symtyp
c
      write(10,'('' '',i6,'' '',i1)') maxit, keep
      write(10,'('' '',4(e19.10E3))') kT, nio, Ldi, kox
      write(10,'('' '',4(e19.10E3))') ksemi, dt, tstart, tstop
      write(10,'('' '',4(e19.10E3))') stpmin, tout, Vs, Vd
      write(10,'('' '',3(e19.10E3))') Vg1, eavg, vmax
      write(10,'('' '',3(e19.10E3))') delt, delV, delE
      write(10,'(1x,a)') chmod
      write(10,'(1x,a)') doptyp
      write(10,'(1x,a)') icchan
      write(10,'(1x,a)') chtyp
      write(10,'(1x,a)') instyp
      write(10,'(1x,a)') debug
      write(10,'(1x,a)') guess
      write(10,'(1x,a)') soltyp
      write(10,'(1x,a)') symtyp
c
c     next dump /struct/, which is:
c
c     common /struct/ lsorce, tsorce, dopsrc, ldrain, tdrain, 
c    +                dopdrn, tox, tepi, length, tsub, dopsub,
c    +                line, space, npulse, defoc, gsol, gdol,
c    +                dopch, srccor, drncor, uins, Qss(Nx)
c     PRECISION lsorce, tsorce, dopsrc, ldrain, tdrain, dopdrn,
c    +     tox, tepi, length, tsub, dopsub, line, space, npulse,
c    +     defoc, gsol, gdol, dopch, uins, Qss
c     character*16 srccor, drncor
c
      write(10,'('' '',4(e19.10E3))') lsorce, tsorce, dopsrc, ldrain
      write(10,'('' '',4(e19.10E3))') tdrain, dopdrn, tox, tepi
      write(10,'('' '',4(e19.10E3))') length, tsub, dopsub, line
      write(10,'('' '',4(e19.10E3))') space, npulse, defoc, gsol
      write(10,'('' '',3(e19.10E3))') gdol, dopch, uins
      do 550 i=1,Nx
  550   write(10,'('' '',e19.10E3)') Qss(i)
      write(10,'(1x,a)') srccor
      write(10,'(1x,a)') drncor
c
c     next dump /curves/, which is:
c
c     common /curves/ Do, dDo, Ad, Bd, mobil, Av, Bv, Cv, Dv, Fv
c     PRECISION Do, dDo, Ad, Bd, mobil, Av, Bv, Cv, Dv, Fv
c
      
      write(10,'('' '',4(e19.10E3))') Do, dDo, Ad, Bd
      write(10,'('' '',4(e19.10E3))') mobil, Av, Bv, Cv
      write(10,'('' '',2(e19.10E3))') Dv, Fv
c
c     finally dump /Pshape/, which is:
c
c     common /Pshape/ Qm, ftw, fwhm, pulse
c     PRECISION Qm, ftw, fwhm
c     logical pulse
c
      write(10,'('' '',3(e19.10E3))') Qm, ftw, fwhm
      write(10,'('' '',l10)') pulse
c
c     and finally, the present time in the simulation, and location
c
      write(10,'('' '',2(e19.10E3))') tnow, xref
c
c     close this file because we are done with it
c
#     ifndef Cyber
      close( 10 )
#     endif
c
c     Now see if we need to write the plottable stuff in cplot format
c
      if ( .not. plot ) return
c
c     We need the offset voltage
c
      Vob = -sign( log( (abs(dopsub) + sqrt(dopsub*dopsub +
     +              4.0*nio*nio))/(2.0*nio) ), dopsub )
c
c     Compensate Vob for Sutherland's and Mock's shift in the reference
c     potential.
c
      if ( symtyp .eq. 'Mock' ) then
        Vob = Vob - sign( log( (abs(dopsrc) + sqrt(dopsrc*dopsrc +
     +                   4.0*nio*nio))/(2.0*nio) ), dopsrc )
      endif
c
c     give me a header
c
      write(9,'('' '','' t = '',f6.2,'' psec.  x, Qn, phis  #3'')')
     +         tnow*1.0e12
c				reset the position
      xtemp = xref
      do 700 i=1,xnodes
c				integrate n(x,j) to get Qn(x)
        ntemp = 0.0
        do 600 j=1,ynodes
c				get the delta y for this point
          if ( j .eq. 1 ) then
            dely = hy(j)*Ldi
          elseif ( j .eq. ynodes ) then
            dely = hy(j-1)*Ldi
          else
            dely = 0.5*( hy(j) + hy(j-1) )*Ldi
          endif
  600     ntemp = ntemp + n(i,j)*nio*dely
c				write out the data point
        write(9,'('' '',f6.3,'' '',e14.6,'' '',f8.4)') 
     +       xtemp*1.0e4, ntemp, (v(i,joxide)-Vob)*kT
c				move the location over
        xtemp = xtemp + hx(i) * Ldi
  700 continue
#     ifndef Cyber
        call flush(9)
#     endif
c
c     all done... return to the caller
c
      return
      end


      subroutine load( filenm, tnow, xref )
c
c     This routine read in, from the source file 'filenm' the
c     "shark.core" style common blocks for restarting when
c     something has gone wrong.
c
c     first the parameters for the sizes of the arrays
c
      integer Nx, Ny
      parameter( Nx=XMAX, Ny=YMAX )
c
c     now for the arguments to this routine
c
      PRECISION tnow, xref
      character*32 filenm
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
c     now some local variables for the loading
c
      integer i, j, k
c
c     load the contents of each common block that is of any use to me.
c
c     open the "core dump" file
c
#     ifndef Cyber
      open( 10, file=filenm, status="UNKNOWN" )
#     endif
      rewind( 10 )
c
c     first, load /device/, which is:
c
c     common /device/ n(Nx,Ny), dop(Nx,Ny), v(Nx,Ny), vgate(Nx),
c    +                drncon(Ny,3), srccon(Ny,3), joxide
c     integer joxide
c     PRECISION n, dop, v, vgate, drncon, srccon
c
      do 100 i=1,Nx
        do 100 j=1,Ny
  100     read(10,*) n(i,j), dop(i,j), v(i,j)
      do 200 i=1,Nx
  200   read(10,*) vgate(i)
      do 300 j=1,Ny
        read(10,*) (srccon(j,k),k=1,3)
  300   read(10,*) (drncon(j,k),k=1,3)
      read(10,*) joxide
c
c     next load /mesh/, which is:
c
c     common /mesh/ hx(Nx), hy(Ny), xnodes, ynodes, xfrac, auto, dVa
c     integer xnodes, ynodes
c     PRECISION hx, hy, xfrac, dVa
c     character*8 auto
c
      do 400 j=1,Nx
  400   read(10,*) hx(j)
      do 500 j=1,Ny
  500   read(10,*) hy(j)
      read(10,*)
     +          xnodes, ynodes, xfrac, dVa
      read(10,'(1x,a)') auto
c
c     next load /params/, which is:
c
c     common /params/ kT, nio, Ldi, kox, chmod, doptyp, ksemi,
c    +                dt, tstart, tstop, stpmin, tout, delt,
c    +                Vs, Vd, vg1, eavg, soltyp, symtyp,
c    +                icchan, chtyp, instyp, vmax, debug,
c    +                guess, delV, delE, maxit, keep
c     integer maxit, keep
c     PRECISION kT, nio, Ldi, kox, ksemi, dt, tstart, tstop,
c    +     stpmin, tout, Vs, Vd, vg1, eavg, vmax, delt,
c    +     delV, delE
c     character*8 chmod, doptyp, icchan, chtyp, instyp, debug, guess,
c    +     soltyp, symtyp
c
      read(10,*) maxit, keep
      read(10,*) kT, nio, Ldi, kox
      read(10,*) ksemi, dt, tstart, tstop
      read(10,*) stpmin, tout, Vs, Vd
      read(10,*) Vg1, eavg, vmax
      read(10,*) delt, delV, delE
      read(10,'(1x,a)') chmod
      read(10,'(1x,a)') doptyp
      read(10,'(1x,a)') icchan
      read(10,'(1x,a)') chtyp
      read(10,'(1x,a)') instyp
      read(10,'(1x,a)') debug
      read(10,'(1x,a)') guess
      read(10,'(1x,a)') soltyp
      read(10,'(1x,a)') symtyp
c
c     next load /struct/, which is:
c
c     common /struct/ lsorce, tsorce, dopsrc, ldrain, tdrain, 
c    +                dopdrn, tox, tepi, length, tsub, dopsub,
c    +                line, space, npulse, defoc, gsol, gdol,
c    +                dopch, srccor, drncor, uins, Qss(Nx)
c     PRECISION lsorce, tsorce, dopsrc, ldrain, tdrain, dopdrn,
c    +     tox, tepi, length, tsub, dopsub, line, space, npulse,
c    +     defoc, gsol, gdol, dopch, uins, Qss
c     character*16 srccor, drncor
c
      read(10,*) lsorce, tsorce, dopsrc, ldrain
      read(10,*) tdrain, dopdrn, tox, tepi
      read(10,*) length, tsub, dopsub, line
      read(10,*) space, npulse, defoc, gsol
      read(10,*) gdol, dopch, uins
      do 550 i=1,Nx
  550   read(10,*) Qss(i)
      read(10,'(1x,a)') srccor
      read(10,'(1x,a)') drncor
c
c     next load /curves/, which is:
c
c     common /curves/ Do, dDo, Ad, Bd, mobil, Av, Bv, Cv, Dv, Fv
c     PRECISION Do, dDo, Ad, Bd, mobil, Av, Bv, Cv, Dv, Fv
c
      
      read(10,*) Do, dDo, Ad, Bd
      read(10,*) mobil, Av, Bv, Cv
      read(10,*) Dv, Fv
c
c     finally load /Pshape/, which is:
c
c     common /Pshape/ Qm, ftw, fwhm, pulse
c     PRECISION Qm, ftw, fwhm
c     logical pulse
c
      read(10,*) Qm, ftw, fwhm
      read(10,*) pulse
c
c     and finally, the present time in the simulation, and location
c
      read(10,*) tnow, xref
c
c     close this file because we are done with it
c
#     ifndef Cyber
      close( 10 )
#     endif
c
c     all done... return to the caller
c
      return
      end
