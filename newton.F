      subroutine newton( error, iter )
c
c     This subroutine solves the system for a consistent set of
c     V and n (assuming p is held in equilibrium), using Newton
c     iteration. The Jacobian is set up, the RHS is set up, and
c     the change to V and phin are then calculated by the chosen
c     matrix solver.
c
c     This solver assumes that the insulator is a perfect one.
c     The potential barrier between the semiconductor and the
c     insulator is infinitely high. The electron concentration
c     in the insulator is held at zero. To keep [J] from going
c     singular I fabricate a simple equation and put it into [J]
c     where the updating for phin of every node in the insulator
c     is to be.
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
      PRECISION one
      parameter( one=1.0 )
c
c     now for the arguments to this routine
c
      logical error
      integer iter
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
c
c     now for the function declarations I need
c
      real dtime
c
c     Now I need to see if the computer I am on can do a
c     direct matrix inversion
c
#     ifdef Cyber
#       define DIRECT
#     endif
#     ifdef NP1
#       define DIRECT
#     endif
#     ifdef Sun3
#       define DIRECT
#     endif
c
c     see if we are to override the default memory model
c
#     ifdef SMALL
#       undef DIRECT
#     endif
c
c     now for some local variables
c
#     ifdef DIRECT
      integer info
#     endif
      logical LOSl, LOSr, LOSb
      integer i, k, siter, inorm, jnorm, l,
     +        LPMl, LPMr, RPMl, RPMr, TPMl, TPMr, BPMl, BPMr,
     +        SPMl, SPMr, Fvl, Fnl
      real tarray(2), sitime, time0, ittime, ntime, sutime, lptime
      PRECISION vl, vr, vt, vb, vij,
     +          hl, hr, ht, hb, vsatl, vsatr, vsatt, vsatb,
     +          nl, nr, nt, nb, nij, norm, Fnorm, a, b, c, d,
     +          fl, flx, fly, fr, frx, fry, ft, ftx, fty, fb, fbx, fby,
     +          Rexp, Lexp, Texp, Bexp, Rexpm1, Lexpm1, Texpm1, Bexpm1
      PRECISION GR, dGRdV, dGRdn, vell, velr, velt, velb,
     + diffl, diffr, difft, diffb, dFndV, dFndVl, dFndVr, dFndVt, 
     + dFndVb,
     + dDldV, dDrdV, dDtdV, dDbdV, dvldV, dvrdV, dvtdV, dvbdV,
     + dDldVl, dDrdVl, dDtdVl, dDbdVl, dvldVl, dvrdVl, dvtdVl, dvbdVl,
     + dDldVr, dDrdVr, dDtdVr, dDbdVr, dvldVr, dvrdVr, dvtdVr, dvbdVr,
     + dDldVt, dDrdVt, dDtdVt, dDbdVt, dvldVt, dvrdVt, dvtdVt, dvbdVt,
     + dDldVb, dDrdVb, dDtdVb, dDbdVb, dvldVb, dvrdVb, dvtdVb, dvbdVb,
     + Ql, Qr, Qt, Qb, dFndn, dFndnl, dFndnr, dFndnt, dFndnb,
     + Jl, Jr, Jt, Jb,
     + dJldV, dJrdV, dJtdV, dJbdV, dJldn, dJrdn, dJtdn, dJbdn,
     + dJldVl, dJrdVl, dJtdVl, dJbdVl, dJldnl, dJrdnl, dJtdnl, dJbdnl,
     + dJldVr, dJrdVr, dJtdVr, dJbdVr, dJldnr, dJrdnr, dJtdnr, dJbdnr,
     + dJldVt, dJrdVt, dJtdVt, dJbdVt, dJldnt, dJrdnt, dJtdnt, dJbdnt,
     + dJldVb, dJrdVb, dJtdVb, dJbdVb, dJldnb, dJrdnb, dJtdnb, dJbdnb,
     + axl, axr, axt, axb, ayl, ayr, ayt, ayb, ul, ur, ut, ub,
     + d2Vdx2, d2Vdy2, rho, expm1
      PRECISION nold(Nx,Ny), nijold, Vob, Voc, diagd, linf
      PRECISION Vso, deln
      logical volt
      character*16 oursol

#     ifdef DIRECT
c
c     Now for the matrix for the direct solvers
c
#     ifdef Cyber
      common abd
      PRECISION abd( (4*Ny+3), (2*Nx*Ny) ), det
#     else
      common abd, ipvt
      integer ipvt( 2*Nx*Ny )
      PRECISION abd( (6*Ny+4), (2*Nx*Ny) )
#     endif
      integer fd, sd
c
c     Calculate the dimensions on the matricies to be passed
c
#     ifdef Cyber
      fd = 4*ynodes + 3
#     else
      fd = 6*ynodes + 4
#     endif
c               keep the size of abd() the same, but don't
c               blow ipvt out of the water
      sd = min( ( (6*Ny+4)*(2*Nx*Ny)/fd ), 2*Nx*Ny )
#     endif
c
c     precalculate the bulk potential
c
      if ( chmod .eq. 'buried' ) then
        Voc = sign( log( (abs(dopch) + sqrt(dopch*dopch +
     +                4.0*nio*nio))/(2.0*nio) ), dopch )
      else
        Voc = -sign( log( (abs(dopsub) + sqrt(dopsub*dopsub +
     +                4.0*nio*nio))/(2.0*nio) ), dopsub )
      endif
      Vob = -sign( log( (abs(dopsub) + sqrt(dopsub*dopsub +
     +                4.0*nio*nio))/(2.0*nio) ), dopsub )
      Vso = sign( log( abs(dopsrc/nio)/2.0 + sqrt( (dopsrc*
     +        dopsrc/(nio*nio) + 4.0)/4.0 )), dopsrc )
c
c     reset the number of iterations and some other loop-used
c     variables
c
      iter = 0
c
c     copy the solver of choice to the solver we will start with
c
      oursol = soltyp
c
c     copy the matrix n(i,j) into nold(i,j) for the time evolution
c
      do 50 k=1,ynodes
        do 50 i=1,xnodes
   50     nold(i,k) = n(i,k)
c
c     reset the elapsed timer for the start of the first pass
c
      ntime = 0.0
      sitime = dtime(tarray)
c
c     Set the locations of the matricies in [J], so as to minimize the
c     bandwidth of [J].
c
      if ( ynodes .lt. xnodes ) then
        LPMl = 1
        LPMr = 2
        TPMl = 3
        TPMr = 4
        SPMl = 5
        SPMr = 6
        BPMl = 7
        BPMr = 8
        RPMl = 9
        RPMr = 10
      else
        TPMl = 1
        TPMr = 2
        LPMl = 3
        LPMr = 4
        SPMl = 5
        SPMr = 6
        RPMl = 7
        RPMr = 8
        BPMl = 9
        BPMr = 10
      endif
c
c     Get the number exp(-1) for speed's sake. We will be using it
c     a lot.
c
      expm1 = exp( -one )
c
c     start a pass, and reset the norm of delta v and F[]
c     for this pass
c
  100 norm = -1.0
      Fnorm = -1.0
c
c     get a mark as to the starting time of this pass - we'll diff
c     the other times in this pass against this for elapsed times.
c
      time0 = dtime(tarray)
c
c     For every node in the system, collect the data needed and fill
c     in [J] and -F so that we can solve them to get the corrections
c     and, in the end, the result.
c
c     Fill it all in, and reset the l infinity norm of [J].
c
      linf = 0.0
      do 200 k=1,ynodes
        do 200 i=1,xnodes
c
c         assume no LOS exists for this node to do
          LOSl = .false.
          LOSr = .false.
          LOSb = .false.
c
c         calculate where I can put the elements into [J] and -F
c
          if ( ynodes .lt. xnodes ) then
            l = (i-1)*ynodes + k
          else
            l = (k-1)*xnodes + i
          endif
          Fnl = 2*l
          Fvl = Fnl - 1
c
c         collect the variables I will need for the simplified equations
c
c         get hl, hr, ht, hb (use symmetry at the boundaries)
c
          hl = hx(i-1)
          if ( i .eq. 1 ) then
            hl = hx(i)
          endif

          hr = hx(i)
          if ( i .eq. xnodes ) then
            hr = hx(i-1)
          endif

          ht = hy(k-1)
          if ( k .eq. 1 ) then
            ht = hy(k)
          endif

          hb = hy(k)
          if ( k .eq. ynodes ) then
            hb = hy(k-1)
          endif
c
c         get vl, vr, vt, vb, vij (using boundary conditions)
c
          vl = v((i-1),k)
          if ( i .eq. 1 ) then
c
c           See if we are simulating a real contact, or a LOS
c
            if ( srccon(k,3) .eq. 1.0 ) then
              vl = srccon(k,1)
            else
              vl = v((i+1),k)
              LOSl = .true.
            endif
          endif

          vr = v((i+1),k)
          if ( i .eq. xnodes ) then
c
c           See if we are simulating a real contact, or a LOS
c
            if ( drncon(k,3) .eq. 1.0 ) then
              vr = drncon(k,1)
            else
              vr = v((i-1),k)
              LOSr = .true.
            endif
          endif

          vt = v(i,(k-1))
          if ( k .eq. 1 ) then
c
c           Only specify the gate voltage where it is not there,
c           otherwise it is the same as where we are
c
            if ( vgate(i) .eq. 1.0e20 ) then
              vt = v(i,k)
            else
              vt = vgate(i)
            endif
          endif

c
c         LOS at the bottom - except for Mock
c
          vb = v(i,(k+1))
          if ( k .eq. ynodes ) then
            if ( symtyp .eq. 'Mock' ) then
              vb = Vob - Vso
            else
              vb = v(i,(k-1))
              LOSb = .true.
            endif
          endif

          vij = v(i,k)
c
c         Get nr, nl, nt, nb, nij, and nijold
c
          nl = n((i-1),k)
          if ( i .eq. 1 ) then
c
c           See if we are simulating a real contact, or a LOS
c
            if ( srccon(k,3) .eq. 1.0 ) then
c               n+ region ohmic contact
              if ( symtyp .eq. 'full' ) nl = dopsrc/nio
c               n in the p, part of ohmic contact
              if ( symtyp .eq. 'Mock' ) nl = (dopsrc/nio)*expm1
c               fixed n(y) for the center
              if ( symtyp .eq. 'center' ) nl = srccon(k,2)
            else
c               LOS boundary condition
              nl = n((i+1),k)
              LOSl = .true.
            endif
          endif

          nr = n((i+1),k)
          if ( i .eq. xnodes ) then
c
c           See if we are simulating a real contact, or a LOS
c
            if ( drncon(k,3) .eq. 1.0 ) then
c               n+ region ohmic contact
              if ( symtyp .eq. 'full' ) nr = dopdrn/nio
c               n in the p, part of ohmic contact
              if ( symtyp .eq. 'Mock' ) nr = (dopdrn/nio)*expm1
c               fixed n(y) for the center
              if ( symtyp .eq. 'center' ) nr = drncon(k,2)
            else
c               LOS boundary condition
              nr = n((i-1),k)
              LOSr = .true.
            endif
          endif

          nt = n(i,(k-1))
c               we can't have electrons above
c               the oxide
          if ( k .eq. 1 ) then
            nt = 0.0
          endif

          nb = n(i,(k+1))
c               LOS at the bottom - except for Mock
          if ( k .eq. ynodes ) then
            if ( symtyp .eq. 'Mock' ) then
              nb = (nio/dopsub)*expm1
            else
              nb = n(i,(k-1))
              LOSb = .true.
            endif
          endif

c               we are always where we are
          nij = n(i,k)
          nijold = nold(i,k)
c               check to see if we are in the oxide
          if ( k .le. (joxide-1) ) then
            nb = 0.0
          endif
          if ( k .le. joxide ) then
            nl = 0.0
            nr = 0.0
            nij = 0.0
            nijold = 0.0
          endif
          if ( k .le. (joxide+1) ) then
            nt = 0.0
          endif
c
c         now I need to compute the coefficients for the electric
c         field functions
c
          a = ( hl*vr - (hl+hr)*vij + hr*vl )/( hl*hr*(hl+hr) )
          b = ( ht*vb - (ht+hb)*vij + hb*vt )/( ht*hb*(ht+hb) )
          c = ( hl*hl*vr + (hr*hr-hl*hl)*vij - hr*hr*vl )/
     +        ( hl*hr*(hl+hr) )
          d = ( ht*ht*vb + (hb*hb-ht*ht)*vij - hb*hb*vt )/
     +        ( hb*ht*(hb+ht) )
c
c         now let's compute the actual normalized electric field
c         components for each of the four points of interest
c
          flx = a*hl - c
          fly = -d
          frx = -a*hr - c
          fry = -d
          ftx = -c
          fty = b*ht - d
          fbx = -c
          fby = -b*hb - d
c
c         we will also need the unit vector decomposition for
c         later decomposing of the velocity
c
          fl = sqrt( flx*flx + fly*fly )
          if ( fl .ne. 0.0 ) then
            axl = flx / fl
            ayl = fly / fl
c               un-normalize the electric field for computing
            fl = fl * kT / Ldi
          else
            axl = (hl - 1.0)/sqrt( (1.0 - hl)**2 + 1.0 )
            ayl = -1.0/sqrt( (1.0 - hl)**2 + 1.0 )
          endif

          fr = sqrt( frx*frx + fry*fry )
          if ( fr .ne. 0.0 ) then
            axr = frx / fr
            ayr = fry / fr
c               un-normalize the electric field for computing
            fr = fr * kT / Ldi
          else
            axr = (-hr - 1.0)/sqrt( (1.0 + hr)**2 + 1.0 )
            ayr = -1.0/sqrt( (1.0 + hr)**2 + 1.0 )
          endif

          ft = sqrt( ftx*ftx + fty*fty )
          if ( ft .ne. 0.0 ) then
            axt = ftx / ft
            ayt = fty / ft
c               un-normalize the electric field for computing
            ft = ft * kT / Ldi
          else
            axt = -1.0/sqrt( 1.0 + (1.0 - ht)**2 )
            ayt = (ht - 1.0)/sqrt( 1.0 + (1.0 - ht)**2 )
          endif

          fb = sqrt( fbx*fbx + fby*fby )
          if ( fb .ne. 0.0 ) then
            axb = fbx / fb
            ayb = fby / fb
c               un-normalize the electric field for computing
            fb = fb * kT / Ldi
          else
            axb = -1.0/sqrt( 1.0 + (1.0 + hb)**2 )
            ayb = (-hb - 1.0)/sqrt( 1.0 + (1.0 + hb)**2 )
          endif
c
c         Now I need to compute the normalized  vell, velr, velt, velb,
c         and the normalized diffl, diffr, difft, diffb (being careful
c         not to compute any impossible results)
c
          vsatl = Av*exp(-fl/Bv) + Cv*(1.0/(1.0 + (fl/Dv)**Fv ))
          ul = mobil*(1.0/sqrt(1.0+(mobil*fl/vsatl)**2))
          vell = ul*fl * Ldi
          ul = ul * kT
          if ( fl .ne. 0.0 ) then
            diffl = Do + dDo*exp( -((log(fl)-Ad)/Bd)**2 )
          else
            diffl = Do
          endif
          vell = vell * axl

          vsatr = Av*exp(-fr/Bv) + Cv*(1.0/(1.0 + (fr/Dv)**Fv ))
          ur = mobil*(1.0/sqrt(1.0+(mobil*fr/vsatr)**2))
          velr = ur*fr * Ldi
          ur = ur * kT
          if ( fr .ne. 0.0 ) then
            diffr = Do + dDo*exp( -((log(fr)-Ad)/Bd)**2 )
          else
            diffr = Do
          endif
          velr = velr * axr

          vsatt = Av*exp(-ft/Bv) + Cv*(1.0/(1.0 + (ft/Dv)**Fv ))
          ut = mobil*(1.0/sqrt(1.0+(mobil*ft/vsatt)**2))
          velt = ut*ft * Ldi
          ut = ut * kT
          if ( ft .ne. 0.0 ) then
            difft = Do + dDo*exp( -((log(ft)-Ad)/Bd)**2 )
          else
            difft = Do
          endif
          velt = velt * ayt

          vsatb = Av*exp(-fb/Bv) + Cv*(1.0/(1.0 + (fb/Dv)**Fv ))
          ub = mobil*(1.0/sqrt(1.0+(mobil*fb/vsatb)**2))
          velb = ub*fb * Ldi
          ub = ub * kT
          if ( fb .ne. 0.0 ) then
            diffb = Do + dDo*exp( -((log(fb)-Ad)/Bd)**2 )
          else
            diffb = Do
          endif
          velb = velb * ayb
c
c         Now compute some exponentials that we are going to be using
c         quite a lot.
c
          Lexpm1 = hl*vell/diffl
          Lexp = exp( Lexpm1 )
          if ( abs(Lexpm1) .gt. 100.0*vmax ) Lexpm1 = Lexp - 1.0
c               check for the divide by zero point
          if ( vell .eq. 0.0 ) then
            Lexp = 1.0
            Lexpm1 = 0.0
          endif

          Rexpm1 = hr*velr/diffr
          Rexp = exp( Rexpm1 )
          if ( abs(Rexpm1) .gt. 100.0*vmax ) Rexpm1 = Rexp - 1.0
c               check for the divide by zero point
          if ( velr .eq. 0.0 ) then
            Rexp = 1.0
            Rexpm1 = 0.0
          endif

          Texpm1 = ht*velt/difft
          Texp = exp( Texpm1 )
          if ( abs(Texpm1) .gt. 100.0*vmax ) Texpm1 = Texp - 1.0
c               check for the divide by zero point
          if ( velt .eq. 0.0 ) then
            Texp = 1.0
            Texpm1 = 0.0
          endif

          Bexpm1 = hb*velb/diffb
          Bexp = exp( Bexpm1 )
          if ( abs(Bexpm1) .gt. 100.0*vmax ) Bexpm1 = Bexp - 1.0
c               check for the divide by zero point
          if ( velb .eq. 0.0 ) then
            Bexp = 1.0
            Bexpm1 = 0.0
          endif
c
c         These are the partials of the Generation-Recombination
c         (SRH) with respect to the variables
c
          GR = 0.0
          dGRdV = 0.0
          dGRdn = 0.0
c
c         Now compute the currents for this node, depending on the
c         location of the node
c
          if ( k .le. joxide ) then
c               we are in the oxide
            Jl = 0.0
            Jr = 0.0
            Jt = 0.0
            Jb = 0.0
          else
c               we are NOT in the oxide
            if ( vell .eq. 0.0 ) then
              Jl = 0.0
            else
              Jl = -vell*( Lexp*nij - nl )/Lexpm1
            endif
            if ( velr .eq. 0.0 ) then
              Jr = 0.0
            else
              Jr = -velr*( Rexp*nr - nij )/Rexpm1
            endif
c               Jt into the oxide is zero
            if ( (k .eq. (joxide+1)) .or. (velt .eq. 0.0) ) then
              Jt = 0.0
            else
              Jt = -velt*( Texp*nij - nt )/Texpm1
            endif
            if ( velb .eq. 0.0 ) then
              Jb = 0.0
            else
              Jb = -velb*( Bexp*nb - nij )/Bexpm1
            endif
          endif
c
c         Now the tough diffusion partials, but be careful not to
c         try to compute log(0)!
c
          if ( fl .ne. 0.0 ) then
            dDldV = ( -(2.0/(Bd*Bd))*(log(fl)-Ad)*(diffl-Do)/fl )
     +              *(kT/Ldi)
c               normalize the field again
            fl = fl * Ldi / kT
            dDldVl = dDldV * ( -(c-a*hl)/(hl*fl) )
            dDldVr = 0.0
c               check for zero x-component
            if ( flx .eq. 0.0 ) then
              dDldVt = 0.0
              dDldVb = 0.0
            else
              dDldVt = dDldV * ( -hb*d/(fl*ht*(ht+hb)) )
              dDldVb = dDldV * ( ht*d/(fl*hb*(hb+ht)) )
            endif
            dDldV = dDldV * ( (c-a*hl)/hl + d*(hb-ht)/(ht*hb) )/fl
          else
            dDldV = 0.0
            dDldVl = 0.0
            dDldVr = 0.0
            dDldVt = 0.0
            dDldVb = 0.0
          endif
c               UN-normalize the field again
          fl = fl * kT / Ldi

          if ( fr .ne. 0.0 ) then
            dDrdV = ( -(2.0/(Bd*Bd))*(log(fr)-Ad)*(diffr-Do)/fr )
     +              *(kT/Ldi)
c               normalize the field again
            fr = fr * Ldi / kT
            dDrdVl = 0.0
            dDrdVr = dDrdV * (c+a*hr)/(fr*hr)
c               check for zero x-component
            if ( frx .eq. 0.0 ) then
              dDrdVt = 0.0
              dDrdVb = 0.0
            else
              dDrdVt = dDrdV * ( -hb*d/(fr*ht*(ht+hb)) )
              dDrdVb = dDrdV * ( ht*d/(fr*hb*(hb+ht)) )
            endif
            dDrdV = dDrdV * ( -(c+a*hr)/hr + d*(hb-ht)/(ht*hb) )/fr
          else
            dDrdV = 0.0
            dDrdVl = 0.0
            dDrdVr = 0.0
            dDrdVt = 0.0
            dDrdVb = 0.0
          endif
c               UN-normalize the field again
          fr = fr * kT / Ldi

          if ( ft .ne. 0.0 ) then
            dDtdV = ( -(2.0/(Bd*Bd))*(log(ft)-Ad)*(difft-Do)/ft )
     +              *(kT/Ldi)
c               normalize the field again
            ft = ft * Ldi / kT
c               check for zero y-component
            if ( fty .eq. 0.0 ) then
              dDtdVl = 0.0
              dDtdVr = 0.0
            else
              dDtdVl = dDtdV * ( -hr*c/(ft*hl*(hl+hr)) )
              dDtdVr = dDtdV * ( hl*c/(ft*hr*(hl+hr)) )
            endif
            dDtdVt = dDtdV * ( -(d-b*ht)/(ft*ht) )
            dDtdVb = 0.0
            dDtdV = dDtdV * ( (d-b*ht)/(ft*ht) )
          else
            dDtdV = 0.0
            dDtdVl = 0.0
            dDtdVr = 0.0
            dDtdVt = 0.0
            dDtdVb = 0.0
          endif
c               UN-normalize the field again
          ft = ft * kT / Ldi

          if ( fb .ne. 0.0 ) then
            dDbdV = ( -(2.0/(Bd*Bd))*(log(fb)-Ad)*(diffb-Do)/fb )
     +              *(kT/Ldi)
c               normalize the field again
            fb = fb * Ldi / kT
c               check for zero y-component
            if ( fby .eq. 0.0 ) then
              dDbdVl = 0.0
              dDbdVr = 0.0
            else
              dDbdVl = dDbdV * ( -hr*c/(fb*hl*(hl+hr)) )
              dDbdVr = dDbdV * ( hl*c/(fb*hr*(hl+hr)) )
            endif
            dDbdVt = 0.0
            dDbdVb = dDbdV * (d+b*hb)/(fb*hb)
            dDbdV = dDbdV * ( -(d+b*hb)/(fb*hb) )
          else
            dDbdV = 0.0
            dDbdVl = 0.0
            dDbdVr = 0.0
            dDbdVt = 0.0
            dDbdVb = 0.0
          endif
c               UN-normalize the field again
          fb = fb * kT / Ldi
c
c         now the tough velocity partials
c
          if ( fl .ne. 0.0 ) then
            Ql = sqrt( 1.0 + ( mobil*fl/vsatl )**2 )
            dvldV = (( (mobil/(Ql*Ql))*( Ql - (mobil*mobil*fl*fl/(Ql*
     +              vsatl*vsatl*vsatl))*( vsatl - fl*(-(Av/Bv)*exp(-fl/
     +              Bv) - (Cv*Fv/Dv)/((fl/Dv)**(1.0-Fv) * (1.0 + (fl/
     +              Dv)**Fv)**2 ) ) ) ) )*kT - ul )*axl
c               normalize the field again
            fl = fl * Ldi / kT
            dvldVl = dvldV * ( -(c-a*hl)/(fl*hl) ) + ul/hl
            dvldVr = 0.0
c               check for zero x-component
            if ( flx .eq. 0.0 ) then
              dvldVt = 0.0
              dvldVb = 0.0
            else
              dvldVt = dvldV * ( -d*hb/(fl*ht*(ht+hb)) )
              dvldVb = dvldV * ( d*ht/(fl*hb*(ht+hb)) )
            endif
            dvldV = dvldV * ( (-a*hl + c)/(fl*hl) + (d/fl)*
     +                      (hb - ht)/(ht*hb) ) - ul/hl
          else
            dvldV = ( mobil*kT - ul )*axl
            dvldVl = dvldV * (hl - 1.0)*1.41421356237/(hl*sqrt(
     +               (1.0 - hl)**2 + 1.0 )) + ul/hl
            dvldVr = 0.0
            dvldVt = 0.0
            dvldVb = 0.0
            dvldV = dvldV * ( (1.0 - hl)*1.41421356237/(hl*sqrt((1.0
     +              - hl)**2 + 1.0))+(hb - ht)*1.41421356237/(ht*hb*
     +              sqrt((1.0 - hl)**2 + 1.0)) ) - ul/hl
          endif

          if ( fr .ne. 0.0 ) then
            Qr = sqrt( 1.0 + ( mobil*fr/vsatr )**2 )
            dvrdV = (( (mobil/(Qr*Qr))*( Qr - (mobil*mobil*fr*fr/(Qr*
     +              vsatr*vsatr*vsatr))*( vsatr - fr*(-(Av/Bv)*exp(-fr/
     +              Bv) - (Cv*Fv/Dv)/((fr/Dv)**(1.0-Fv) * (1.0 + (fr/
     +              Dv)**Fv)**2 ) ) ) ) )*kT - ur )*axr
c               normalize the field again
            fr = fr * Ldi / kT
            dvrdVl = 0.0
            dvrdVr = dvrdV * (c+a*hr)/(fr*hr) - ur/hr
c               check for zero x-component
            if ( frx .eq. 0.0 ) then
              dvrdVt = 0.0
              dvrdVb = 0.0
            else
              dvrdVt = dvrdV * ( -hb*d/(fr*ht*(ht+hb)) )
              dvrdVb = dvrdV * ( ht*d/(fr*hb*(ht+hb)) )
            endif
            dvrdV = dvrdV * ( -(c+a*hr)/hr + d*(hb-ht)/(ht*hb) )/fr
     +                      + ur/hr
          else
            dvrdV = ( mobil*kT - ur )*axr
            dvrdVl = 0.0
            dvrdVr = dvrdV * (1.0 + hr)*1.41421356237/(hr*sqrt(1.0
     +               + (1.0 + hr)**2)) - ur/hr
            dvrdVt = 0.0
            dvrdVb = 0.0
            dvrdV = dvrdV * ( -(1.0+hr)/hr + (hb-ht)/(hb*ht) )*
     +              1.41421356237 / sqrt( 1.0 + (1.0+hr)**2 ) + ur/hr
          endif

          if ( ft .ne. 0.0 ) then
            Qt = sqrt( 1.0 + ( mobil*ft/vsatt )**2 )
            dvtdV = (( (mobil/(Qt*Qt))*( Qt - (mobil*mobil*ft*ft/(Qt*
     +              vsatt*vsatt*vsatt))*( vsatt - ft*(-(Av/Bv)*exp(-ft/
     +              Bv) - (Cv*Fv/Dv)/((ft/Dv)**(1.0-Fv) * (1.0 + (ft/
     +              Dv)**Fv)**2 ) ) ) ) )*kT - ut )*ayt
c               normalize the field again
            ft = ft * Ldi / kT
c               check for zero y-component
            if ( fty .eq. 0.0 ) then
              dvtdVl = 0.0
              dvtdVr = 0.0
            else
              dvtdVl = dvtdV * ( -hr*c/(ft*hl*(hl+hr)) )
              dvtdVr = dvtdV * ( hl*c/(ft*hr*(hl+hr)) )
            endif
            dvtdVt = dvtdV * ( -(d-b*ht)/(ft*ht) ) + ut/ht
            dvtdVb = 0.0
            dvtdV = dvtdV * ( (d-b*ht)/(ft*ht) ) - ut/ht
          else
            dvtdV = ( mobil*kT - ut )*ayt
            dvtdVl = 0.0
            dvtdVr = 0.0
            dvtdVt = dvtdV * ( -(d-b*ht)*1.41421356237/(ht*sqrt(
     +               1.0 + (1.0 - ht)**2 )) ) + ut/ht
            dvtdVb = 0.0
            dvtdV = dvtdV * ( (1.0-ht)/ht )*
     +              1.41421356237 / sqrt( 1.0 + (1.0-ht)**2 ) - ut/ht
          endif

          if ( fb .ne. 0.0 ) then
            Qb = sqrt( 1.0 + ( mobil*fb/vsatb )**2 )
            dvbdV = (( (mobil/(Qb*Qb))*( Qb - (mobil*mobil*fb*fb/(Qb*
     +              vsatb*vsatb*vsatb))*( vsatb - fb*(-(Av/Bv)*exp(-fb/
     +              Bv) - (Cv*Fv/Dv)/((fb/Dv)**(1.0-Fv) * (1.0 + (fb/
     +              Dv)**Fv)**2 ) ) ) ) )*kT - ub )*ayb
c               normalize the field again
            fb = fb * Ldi / kT
c               check for zero y-component
            if ( fby .eq. 0.0 ) then
              dvbdVl = 0.0
              dvbdVr = 0.0
            else
              dvbdVl = dvbdV * ( -hr*c/(fb*hl*(hl+hr)) )
              dvbdVr = dvbdV * ( hl*c/(fb*hr*(hl+hr)) )
            endif
            dvbdVt = 0.0
            dvbdVb = dvbdV * (d+b*hb)/(fb*hb) - ub/hb
            dvbdV = dvbdV * ( -(d+b*hb)/(fb*hb) ) + ub/hb
          else
            dvbdV = ( mobil*kT - ub )*ayb
            dvbdVl = 0.0
            dvbdVr = 0.0
            dvbdVt = 0.0
            dvbdVb = dvbdV * (1.0 + hb)*1.41421356237/(hb*sqrt(
     +               1.0 + (1.0 + hb)**2 ) ) - ub/hb
            dvbdV = dvbdV * ( -(1.0+hb)/hb )*
     +              1.41421356237 / sqrt( 1.0 + (1.0+hb)**2 ) + ub/hb
          endif
c
c         now the current potentials, being careful to use the proper
c         form for the region we are in
c
          if ( k .le. joxide ) then
c               we are in the insulator
            dJrdV = 0.0
            dJrdVl = 0.0
            dJrdVr = 0.0
            dJrdVt = 0.0
            dJrdVb = 0.0
            dJrdn = 0.0
            dJrdnl = 0.0
            dJrdnr = 0.0
            dJrdnt = 0.0
            dJrdnb = 0.0

            dJldV = 0.0
            dJldVl = 0.0
            dJldVr = 0.0
            dJldVt = 0.0
            dJldVb = 0.0
            dJldn = 0.0
            dJldnl = 0.0
            dJldnr = 0.0
            dJldnt = 0.0
            dJldnb = 0.0

            dJtdV = 0.0
            dJtdVl = 0.0
            dJtdVr = 0.0
            dJtdVt = 0.0
            dJtdVb = 0.0
            dJtdn = 0.0
            dJtdnl = 0.0
            dJtdnr = 0.0
            dJtdnt = 0.0
            dJtdnb = 0.0

            dJbdV = 0.0
            dJbdVl = 0.0
            dJbdVr = 0.0
            dJbdVt = 0.0
            dJbdVb = 0.0
            dJbdn = 0.0
            dJbdnl = 0.0
            dJbdnr = 0.0
            dJbdnt = 0.0
            dJbdnb = 0.0
          else
c
c               we are NOT in the insulator
c
            if ( velr .ne. 0.0 ) then
              dJrdV = ( velr*Rexp*hr/(diffr*diffr*Rexpm1) )*
     +                ( (Rexp*nr-nij)/Rexpm1 - nr )*
     +                ( diffr*dvrdV - velr*dDrdV ) -
     +                dvrdV*( Rexp*nr - nij )/Rexpm1
              dJrdVl = ( velr*Rexp*hr/(diffr*diffr*Rexpm1) )*
     +                ( (Rexp*nr-nij)/Rexpm1 - nr )*
     +                ( diffr*dvrdVl - velr*dDrdVl ) -
     +                dvrdVl*( Rexp*nr - nij )/Rexpm1
              dJrdVr = ( velr*Rexp*hr/(diffr*diffr*Rexpm1) )*
     +                ( (Rexp*nr-nij)/Rexpm1 - nr )*
     +                ( diffr*dvrdVr - velr*dDrdVr ) -
     +                dvrdVr*( Rexp*nr - nij )/Rexpm1
              dJrdVt = ( velr*Rexp*hr/(diffr*diffr*Rexpm1) )*
     +                ( (Rexp*nr-nij)/Rexpm1 - nr )*
     +                ( diffr*dvrdVt - velr*dDrdVt ) -
     +                dvrdVt*( Rexp*nr - nij )/Rexpm1
              dJrdVb = ( velr*Rexp*hr/(diffr*diffr*Rexpm1) )*
     +                ( (Rexp*nr-nij)/Rexpm1 - nr )*
     +                ( diffr*dvrdVb - velr*dDrdVb ) -
     +                dvrdVb*( Rexp*nr - nij )/Rexpm1
              dJrdn = velr/Rexpm1
              dJrdnl = 0.0
              dJrdnr = -velr*Rexp/Rexpm1
              dJrdnt = 0.0
              dJrdnb = 0.0
            else
              dJrdV = -dvrdV * nij
              dJrdVl = -dvrdVl * nij
              dJrdVr = -dvrdVr * nij
              dJrdVt = -dvrdVt * nij
              dJrdVb = -dvrdVb * nij
              dJrdn = diffr/hr
              dJrdnl = 0.0
              dJrdnr = -diffr/hr - velr
              dJrdnt = 0.0
              dJrdnb = 0.0
            endif

            if ( vell .ne. 0.0 ) then
              dJldV = ( vell*Lexp*hl/(diffl*diffl*Lexpm1) )*
     +                ( (Lexp*nij-nl)/Lexpm1 - nij )*
     +                ( diffl*dvldV - vell*dDldV ) -
     +                dvldV*( Lexp*nij - nl )/Lexpm1
              dJldVl = ( vell*Lexp*hl/(diffl*diffl*Lexpm1) )*
     +                ( (Lexp*nij-nl)/Lexpm1 - nij )*
     +                ( diffl*dvldVl - vell*dDldVl ) -
     +                dvldVl*( Lexp*nij - nl )/Lexpm1
              dJldVr = ( vell*Lexp*hl/(diffl*diffl*Lexpm1) )*
     +                ( (Lexp*nij-nl)/Lexpm1 - nij )*
     +                ( diffl*dvldVr - vell*dDldVr ) -
     +                dvldVr*( Lexp*nij - nl )/Lexpm1
              dJldVt = ( vell*Lexp*hl/(diffl*diffl*Lexpm1) )*
     +                ( (Lexp*nij-nl)/Lexpm1 - nij )*
     +                ( diffl*dvldVt - vell*dDldVt ) -
     +                dvldVt*( Lexp*nij - nl )/Lexpm1
              dJldVb = ( vell*Lexp*hl/(diffl*diffl*Lexpm1) )*
     +                ( (Lexp*nij-nl)/Lexpm1 - nij )*
     +                ( diffl*dvldVb - vell*dDldVb ) -
     +                dvldVb*( Lexp*nij - nl )/Lexpm1
              dJldn = -vell*Lexp/Lexpm1
              dJldnl = vell/Lexpm1
              dJldnr = 0.0
              dJldnt = 0.0
              dJldnb = 0.0
            else
              dJldV = -dvldV * nij
              dJldVl = -dvldVl * nij
              dJldVr = -dvldVr * nij
              dJldVt = -dvldVt * nij
              dJldVb = -dvldVb * nij
              dJldn = -diffl/hl - vell
              dJldnl = diffl/hl
              dJldnr = 0.0
              dJldnt = 0.0
              dJldnb = 0.0
            endif

            if ( velb .ne. 0.0 ) then
              dJbdV = ( velb*Bexp*hb/(diffb*diffb*Bexpm1) )*
     +                ( (Bexp*nb-nij)/Bexpm1 - nb )*
     +                ( diffb*dvbdV - velb*dDbdV ) -
     +                dvbdV*( Bexp*nb - nij )/Bexpm1
              dJbdVl = ( velb*Bexp*hb/(diffb*diffb*Bexpm1) )*
     +                ( (Bexp*nb-nij)/Bexpm1 - nb )*
     +                ( diffb*dvbdVl - velb*dDbdVl ) -
     +                dvbdVl*( Bexp*nb - nij )/Bexpm1
              dJbdVr = ( velb*Bexp*hb/(diffb*diffb*Bexpm1) )*
     +                ( (Bexp*nb-nij)/Bexpm1 - nb )*
     +                ( diffb*dvbdVr - velb*dDbdVr ) -
     +                dvbdVr*( Bexp*nb - nij )/Bexpm1
              dJbdVt = ( velb*Bexp*hb/(diffb*diffb*Bexpm1) )*
     +                ( (Bexp*nb-nij)/Bexpm1 - nb )*
     +                ( diffb*dvbdVt - velb*dDbdVt ) -
     +                dvbdVt*( Bexp*nb - nij )/Bexpm1
              dJbdVb = ( velb*Bexp*hb/(diffb*diffb*Bexpm1) )*
     +                ( (Bexp*nb-nij)/Bexpm1 - nb )*
     +                ( diffb*dvbdVb - velb*dDbdVb ) -
     +                dvbdVb*( Bexp*nb - nij )/Bexpm1
              dJbdn = velb/Bexpm1
              dJbdnl = 0.0
              dJbdnr = 0.0
              dJbdnt = 0.0
              dJbdnb = -velb*Bexp/Bexpm1
            else
              dJbdV = -dvbdV * nij
              dJbdVl = -dvbdVl * nij
              dJbdVr = -dvbdVr * nij
              dJbdVt = -dvbdVt * nij
              dJbdVb = -dvbdVb * nij
              dJbdn = diffb/hb
              dJbdnl = 0.0
              dJbdnr = 0.0
              dJbdnt = 0.0
              dJbdnb = -diffb/hb - velb
            endif
c
c           check to see if the center node is in the insulator,
c           so don't let current go streaming into the insulator
c
            if ( k .eq. (joxide+1) ) then
              dJtdV = 0.0
              dJtdVl = 0.0
              dJtdVr = 0.0
              dJtdVt = 0.0
              dJtdVb = 0.0
              dJtdn = 0.0
              dJtdnl = 0.0
              dJtdnr = 0.0
              dJtdnt = 0.0
              dJtdnb = 0.0
            else
              if ( velt .ne. 0.0 ) then
                dJtdV = ( velt*Texp*ht/(difft*difft*Texpm1) )*
     +                  ( (Texp*nij-nt)/Texpm1 - nij )*
     +                  ( difft*dvtdV - velt*dDtdV ) -
     +                  dvtdV*( Texp*nij - nt )/Texpm1
                dJtdVl = ( velt*Texp*ht/(difft*difft*Texpm1) )*
     +                  ( (Texp*nij-nt)/Texpm1 - nij )*
     +                  ( difft*dvtdVl - velt*dDtdVl ) -
     +                  dvtdVl*( Texp*nij - nt )/Texpm1
                dJtdVr = ( velt*Texp*ht/(difft*difft*Texpm1) )*
     +                  ( (Texp*nij-nt)/Texpm1 - nij )*
     +                  ( difft*dvtdVr - velt*dDtdVr ) -
     +                  dvtdVr*( Texp*nij - nt )/Texpm1
                dJtdVt = ( velt*Texp*ht/(difft*difft*Texpm1) )*
     +                  ( (Texp*nij-nt)/Texpm1 - nij )*
     +                  ( difft*dvtdVt - velt*dDtdVt ) -
     +                  dvtdVt*( Texp*nij - nt )/Texpm1
                dJtdVb = ( velt*Texp*ht/(difft*difft*Texpm1) )*
     +                  ( (Texp*nij-nt)/Texpm1 - nij )*
     +                  ( difft*dvtdVb - velt*dDtdVb ) -
     +                  dvtdVb*( Texp*nij - nt )/Texpm1
                dJtdn = -velt*Texp/Texpm1
                dJtdnl = 0.0
                dJtdnr = 0.0
                dJtdnt = velt/Texpm1
                dJtdnb = 0.0
              else
                dJtdV = -dvtdV * nij
                dJtdVl = -dvtdVl * nij
                dJtdVr = -dvtdVr * nij
                dJtdVt = -dvtdVt * nij
                dJtdVb = -dvtdVb * nij
                dJtdn = -difft/ht - velt
                dJtdnl = 0.0
                dJtdnr = 0.0
                dJtdnt = difft/ht
                dJtdnb = 0.0
              endif
            endif
          endif
c
c         now the main partials themselves
c
          dFndV = (2.0/(hr+hl))*(dJrdV - dJldV) + (2.0/(ht+hb))*
     +            (dJbdV - dJtdV) + dGRdV
          dFndVl = (2.0/(hr+hl))*(dJrdVl - dJldVl) + (2.0/(ht+hb))*
     +            (dJbdVl - dJtdVl)
          dFndVr = (2.0/(hr+hl))*(dJrdVr - dJldVr) + (2.0/(ht+hb))*
     +            (dJbdVr - dJtdVr)
          dFndVt = (2.0/(hr+hl))*(dJrdVt - dJldVt) + (2.0/(ht+hb))*
     +            (dJbdVt - dJtdVt)
          dFndVb = (2.0/(hr+hl))*(dJrdVb - dJldVb) + (2.0/(ht+hb))*
     +            (dJbdVb - dJtdVb)
          dFndn = (2.0/(hr+hl))*(dJrdn - dJldn) + (2.0/(ht+hb))*
     +            (dJbdn - dJtdn) + dGRdn
          dFndnl = (2.0/(hr+hl))*(dJrdnl - dJldnl) + (2.0/(ht+hb))*
     +            (dJbdnl - dJtdnl)
          dFndnr = (2.0/(hr+hl))*(dJrdnr - dJldnr) + (2.0/(ht+hb))*
     +            (dJbdnr - dJtdnr)
          dFndnt = (2.0/(hr+hl))*(dJrdnt - dJldnt) + (2.0/(ht+hb))*
     +            (dJbdnt - dJtdnt)
          dFndnb = (2.0/(hr+hl))*(dJrdnb - dJldnb) + (2.0/(ht+hb))*
     +            (dJbdnb - dJtdnb)
c
c         Now I can fill in the LPM sub-matrix of [J]
c
c               upper left corner element
          J( Fvl, LPMl ) = 1.0/(0.5*hl*(hl+hr))
c               upper right corner element
          J( Fvl, LPMr ) = 0.0
c               lower left corner element
          J( Fnl, LPMl ) = dFndVl
c               lower right corner element
          J( Fnl, LPMr ) = -nl*dFndnl
c
c         Now I can fill in the TPM sub-matrix of [J]
c
c               upper left corner element
          if ( k .eq. joxide ) then
            J( Fvl, TPMl ) = 2.0/(ht*hb*(1.0+(ksemi/kox)*(ht/hb)))
          else
            J( Fvl, TPMl ) = 2.0/(ht*(ht+hb))
          endif
c               upper right corner element
          J( Fvl, TPMr ) = 0.0
c               lower left corner element
          J( Fnl, TPMl ) = dFndVt
c               lower right corner element
          J( Fnl, TPMr ) = -nt*dFndnt
c
c         Now I can fill in the SPM sub-matrix of [J]
c
c               upper left corner element
          if ( k .lt. joxide ) then
            J( Fvl, SPMl ) = -2.0*(1.0/(hl*hr) + 1.0/(ht*hb))
          elseif ( k .eq. joxide ) then
            J( Fvl, SPMl ) = -2.0*(1.0/(hl*hr) + 1.0/(ht*hb))
          else
            if ( symtyp .eq. 'Mock' ) then
              J( Fvl, SPMl ) = -2.0*(1.0/(hl*hr) + 1.0/(ht*hb))
     +                         - exp(-vij-log(dopsrc/nio)) - nij
            else
              J( Fvl, SPMl ) = -2.0*(1.0/(hl*hr) + 1.0/(ht*hb))
     +                         - exp(-vij) - nij
            endif
          endif
c               upper right corner element
          J( Fvl, SPMr ) = nij
c               lower left corner element
          J( Fnl, SPMl ) = dFndV
c               lower right corner element
          if ( k .le. joxide ) then
            J( Fnl, SPMr ) = 1.0
          else
            J( Fnl, SPMr ) = -nij*(dFndn + Ldi*Ldi/dt)
          endif
c
c         Now I can fill in the BPM sub-matrix of [J]
c
c               upper left corner element
          if ( k .eq. joxide ) then
            J( Fvl, BPMl ) = 2.0/(hb*(ht+hb))*(1.0+(ht/hb))/
     +                       ((kox/ksemi) + (ht/hb))
          else
            J( Fvl, BPMl ) = 2.0/(hb*(ht+hb))
          endif
c               upper right corner element
          J( Fvl, BPMr ) = 0.0
c               lower left corner element
          J( Fnl, BPMl ) = dFndVb
c               lower right corner element
          J( Fnl, BPMr ) = -nb*dFndnb
c
c         Now I can fill in the RPM sub-matrix of [J]
c
c               upper left corner element
          J( Fvl, RPMl ) = 1.0/(0.5*hr*(hl+hr))
c               upper right corner element
          J( Fvl, RPMr ) = 0.0
c               lower left corner element
          J( Fnl, RPMl ) = dFndVr
c               lower right corner element
          J( Fnl, RPMr ) = -nr*dFndnr
c
c         Now that all the partials are calculated, we need to fix
c         them for the presence of L.O.S. boundary conditions, and
c         imaginary points from the gate and ohmic contact.
c
          if ( i .eq. 1 ) then
            if ( symtyp .ne. 'center' ) then
              if ( LOSl ) then
                J( Fvl, RPMl ) = J( Fvl, RPMl ) + J( Fvl, LPMl )
                J( Fvl, RPMr ) = J( Fvl, RPMr ) + J( Fvl, LPMr )
                J( Fnl, RPMl ) = J( Fnl, RPMl ) + J( Fnl, LPMl )
                J( Fnl, RPMr ) = J( Fnl, RPMr ) + J( Fnl, LPMr )
              endif
            endif
c               zero out LPM
            J( Fvl, LPMl ) = 0.0
            J( Fvl, LPMr ) = 0.0
            J( Fnl, LPMl ) = 0.0
            J( Fnl, LPMr ) = 0.0
          endif

          if ( k .eq. 1 ) then
c               zero out TPM
            J( Fvl, TPMl ) = 0.0
            J( Fvl, TPMr ) = 0.0
            J( Fnl, TPMl ) = 0.0
            J( Fnl, TPMr ) = 0.0
          endif

          if ( k .eq. ynodes ) then
c               add partials for LOS
            if ( LOSb ) then
              J( Fvl, TPMl ) = J( Fvl, TPMl ) + J( Fvl, BPMl )
              J( Fvl, TPMr ) = J( Fvl, TPMr ) + J( Fvl, BPMr )
              J( Fnl, TPMl ) = J( Fnl, TPMl ) + J( Fnl, BPMl )
              J( Fnl, TPMr ) = J( Fnl, TPMr ) + J( Fnl, BPMr )
            endif
c               zero out BPM
            J( Fvl, BPMl ) = 0.0
            J( Fvl, BPMr ) = 0.0
            J( Fnl, BPMl ) = 0.0
            J( Fnl, BPMr ) = 0.0
          endif

          if ( i .eq. xnodes ) then
            if ( symtyp .ne. 'center' ) then
              if ( LOSr ) then
                J( Fvl, LPMl ) = J( Fvl, LPMl ) + J( Fvl, RPMl )
                J( Fvl, LPMr ) = J( Fvl, LPMr ) + J( Fvl, RPMr )
                J( Fnl, LPMl ) = J( Fnl, LPMl ) + J( Fnl, RPMl )
                J( Fnl, LPMr ) = J( Fnl, LPMr ) + J( Fnl, RPMr )
              endif
            endif
c               zero out RPM
            J( Fvl, RPMl ) = 0.0
            J( Fvl, RPMr ) = 0.0
            J( Fnl, RPMl ) = 0.0
            J( Fnl, RPMr ) = 0.0
          endif
c
c         Now I can fill in the section of -F for this node
c
          if ( k .lt. joxide ) then
c               -Fv terms for this node in the oxide
            d2Vdx2 = (hr*vl-(hl+hr)*vij+hl*vr)/(0.5*hr*hl*(hl+hr))
            d2Vdy2 = (ht*vb-(ht+hb)*vij+hb*vt)/(0.5*ht*hb*(hb+ht))
            rho = 0.0
c               -Fn for this node in the oxide
            F( Fnl ) = 0.0
          elseif ( k .eq. joxide ) then
c               -Fv terms for this node at the interface
            d2Vdx2 = (hr*vl-(hl+hr)*vij+hl*vr)/(0.5*hr*hl*(hl+hr))
            d2Vdy2 = (2.0/(ht*hb*(ht+hb)))*( vt*(ht+hb)/(1.0+
     +               (ksemi/kox)*(ht/hb)) - vij*(ht+hb) + vb*ht*
     +               (1.0+ht/hb)/(kox/ksemi+ht/hb) )
            rho = 2.0 * Qss(i)/(nio*hb*Ldi)
c               -Fn for this node at the interface
            F( Fnl ) = 0.0
          elseif ( (k .eq. ynodes) .and. (symtyp .eq. 'Mock') ) then
c               -Fv terms for this node in the semi.
            d2Vdx2 = (hr*vl-(hl+hr)*vij+hl*vr)/(0.5*hr*hl*(hl+hr))
            d2Vdy2 = (ht*vb-(ht+hb)*vij+hb*vt)/(0.5*ht*hb*(hb+ht))
            rho = 0.0
c               -Fn for this node in the semi.
            F( Fnl ) = 0.0
          else
c               -Fv terms for this node in the semi.
            d2Vdx2 = (hr*vl-(hl+hr)*vij+hl*vr)/(0.5*hr*hl*(hl+hr))
            d2Vdy2 = (ht*vb-(ht+hb)*vij+hb*vt)/(0.5*ht*hb*(hb+ht))
            if ( symtyp .eq. 'Mock' ) then
              rho = exp(-vij-log(dopsrc/nio)) + dop(i,k) - nij
            else
              rho = exp(-vij) + dop(i,k) - nij
            endif
c               -Fn for this node in the semi.
            F( Fnl ) = -( 2.0*(Jr-Jl)/(hl+hr) + 2.0*(Jb-Jt)/(hb+ht) +
     +                  GR + (nij - nijold)*Ldi*Ldi/dt )
          endif
c
c         See if Fv is zero for our purposes
c
          if ( abs(rho) .gt. (1.0/nio) ) then
            if ( abs( ( d2Vdx2 + d2Vdy2 - rho )/rho ) .lt.
     +           (vmax/10.0) ) then
              F( Fvl ) = 0.0
            else
              F( Fvl ) = -( d2Vdx2 + d2Vdy2 + rho )
            endif
          else
            if ( abs(d2Vdx2 + d2Vdy2) .lt. (1.0/nio) ) then
              F( Fvl ) = 0.0
            else
              F( Fvl ) = -( d2Vdx2 + d2Vdy2 )
            endif
          endif
c
c         See if Fn is zero for our purposes
c
          if ( nij .lt. (1.0/nio) ) F( Fnl ) = 0.0
c
c         See if the entire row of J is zero due to a
c         depleation region of electrons.
c
          if ( ( J(Fnl, SPMr) .eq. 0.0 ) .and.
     +         ( F(Fnl) .eq. 0.0 ) ) then
            J(Fnl, SPMr) = 1.0
          endif
c
c         Gather the stats on the norm of F[]. We'll
c         use it later to see if we're "solved".
c
          Fnorm = max( Fnorm, abs(F(Fvl)), abs(F(Fvl)) )
c
c         Now normalize the rows of [J] w.r.t. the diagonal
c         elements and check for diagonal dominance, and l
c         infinity norm
c
          diagd = 0.0
          if ( J(Fvl, SPMl) .ne. 0.0 ) then
            J(Fvl, LPMl) = J(Fvl, LPMl) / J(Fvl, SPMl)
            J(Fvl, LPMr) = J(Fvl, LPMr) / J(Fvl, SPMl)
            J(Fvl, TPMl) = J(Fvl, TPMl) / J(Fvl, SPMl)
            J(Fvl, TPMr) = J(Fvl, TPMr) / J(Fvl, SPMl)
            J(Fvl, SPMr) = J(Fvl, SPMr) / J(Fvl, SPMl)
            J(Fvl, BPMl) = J(Fvl, BPMl) / J(Fvl, SPMl)
            J(Fvl, BPMr) = J(Fvl, BPMr) / J(Fvl, SPMl)
            J(Fvl, RPMl) = J(Fvl, RPMl) / J(Fvl, SPMl)
            J(Fvl, RPMr) = J(Fvl, RPMr) / J(Fvl, SPMl)
            F( Fvl ) = F( Fvl ) / J(Fvl, SPMl)
            J(Fvl, SPMl) = 1.0
            diagd = max( abs(J(Fvl, LPMl)), abs(J(Fvl, LPMr)),
     +                   abs(J(Fvl, TPMl)), abs(J(Fvl, TPMr)),
     +                   1.0,               abs(J(Fvl, SPMr)),
     +                   abs(J(Fvl, BPMl)), abs(J(Fvl, BPMr)),
     +                   abs(J(Fvl, RPMl)), abs(J(Fvl, RPMr)) )
            linf = max( linf, abs(diagd) )
          else
            write(6,'('' '',''F('',i5,'')=Fv('',i3,'','',i3,
     +                      '') has a zero diagonal.'')') Fvl, i, k
            if ( debug .ne. 'no' ) then
              write(6,'(10(g10.3,1x),3x,g10.3)') (J(Fvl,siter),
     +              siter=1,10), F(Fvl)
            endif
          endif
c
c         check for diagonal dominance
c
          if ( ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) .and.
     +         ( diagd .ge. 1.0 ) .and. ( soltyp .ne. 'linpk' ) ) then
            write(6,'('' '',''Trouble! The system matrix is non-'',
     +            ''diagonally dominant. An off-diagonal of Fv''/
     +                '' '',''is '',g9.3,'' times the diagonal '',
     +            ''at node '',i3,'','',i3,''.'')') diagd, i, k
          endif

          diagd = 0.0
          if ( J(Fnl, SPMr) .ne. 0.0 ) then
            J(Fnl, LPMl) = J(Fnl, LPMl) / J(Fnl, SPMr)
            J(Fnl, LPMr) = J(Fnl, LPMr) / J(Fnl, SPMr)
            J(Fnl, TPMl) = J(Fnl, TPMl) / J(Fnl, SPMr)
            J(Fnl, TPMr) = J(Fnl, TPMr) / J(Fnl, SPMr)
            J(Fnl, SPMl) = J(Fnl, SPMl) / J(Fnl, SPMr)
            J(Fnl, BPMl) = J(Fnl, BPMl) / J(Fnl, SPMr)
            J(Fnl, BPMr) = J(Fnl, BPMr) / J(Fnl, SPMr)
            J(Fnl, RPMl) = J(Fnl, RPMl) / J(Fnl, SPMr)
            J(Fnl, RPMr) = J(Fnl, RPMr) / J(Fnl, SPMr)
            F( Fnl ) = F( Fnl ) / J(Fnl, SPMr)
            J(Fnl, SPMr) = 1.0
            diagd = max( abs(J(Fnl, LPMl)), abs(J(Fnl, LPMr)),
     +                   abs(J(Fnl, TPMl)), abs(J(Fnl, TPMr)),
     +                   abs(J(Fnl, SPMl)), 1.0,
     +                   abs(J(Fnl, BPMl)), abs(J(Fnl, BPMr)),
     +                   abs(J(Fnl, RPMl)), abs(J(Fnl, RPMr)) )
            linf = max( linf, abs(diagd) )
          else
            write(6,'('' '',''F('',i5,'')=Fn('',i3,'','',i3,
     +                      '') has a zero diagonal.'')') Fnl, i, k
            if ( debug .ne. 'no' ) then
              write(6,'(10(g10.3,1x),3x,g10.3)') (J(Fnl,siter),
     +              siter=1,10), F(Fnl)
            endif
          endif
c
c         check for diagonal dominance
c
          if ( ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) .and.
     +         ( diagd .ge. 1.0 ) .and. ( soltyp .ne. 'linpk' ) ) then
            write(6,'('' '',''Trouble! The system matrix is non-'',
     +            ''diagonally dominant. An off-diagonal of Fn''/
     +                '' '',''is '',g9.3,'' times the diagonal '',
     +            ''at node '',i3,'','',i3,''.'')') diagd, i, k
          endif
c
c         All done filling out [J] and -F!
c
          if ( debug .eq. 'yes' ) then
            write(6,'(10(g10.3,1x),3x,g10.3)') (J(Fvl,siter),
     +            siter=1,10), F(Fvl)
            write(6,'(10(g10.3,1x),3x,g10.3)') (J(Fnl,siter),
     +            siter=1,10), F(Fnl)
#           ifndef Cyber
              call flush(6)
#           endif
          endif
  200 continue
c
c     See if the norm of F[] indicates that we're solved and
c     can bypass the solution of this Ju=-F
c
      if ( Fnorm .lt. (vmax/1.0e3) ) then
        if ( debug .ne. 'no' ) then
          write(6,'('' '',''The norm of F indicates that '',
     +                    ''the solution has been reached.''/
     +              '' '',''We are therefore skipping the '',
     +                    ''solution of this matrix.'')')
        endif
        goto 400
      endif
c
c     Print out the l infinity norm if debug is anything BUT 'no'
c
      if ( ( debug .ne. 'no' ) .and. ( linf .gt. 1.0 ) .and.
     +     ( soltyp .ne. 'linpk' ) ) then
        write(6,'('' '',''The l inf. norm of [J] is '',g10.3,
     +                  '' which may cause convergence '',
     +                  ''trouble.'')') linf
      endif
c
c     Get the setup time (we'll take out the bias later)
c
      sutime = dtime(tarray)
c
c     Now we call the solution system for this mess, depending on
c     how we want it solved
c
c     Is it to be Orthomin(k)?
c
      if ( oursol .eq. 'ortho' ) then
        call ortho( error, siter )
        if ( error ) then
          if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
            write(6,2)
    2       format(' ','*** Error! ***')
            write(6,3)
    3       format(' ','    While trying to solve the system of',
     +           ' equations, Shark has detected an error in')
            write(6,4)
    4       format(' ','    the solution itself. This has made it',
     +           ' impossible to solve the system.')
            write(6,'('' '',''    There were '',i4,'' iterations on'',
     +                      '' the Orthomin('',i1,'') algorithim '',
     +                      ''completed.'')') siter, keep
#           ifndef Cyber
              call flush(6)
#           endif
          endif
          error = .true.
          return
        endif
c
c     Is it to be the Gauss-Seidel method?
c
      elseif ( oursol .eq. 'gauss' ) then
        call gauss( error, siter )
        if ( error ) then
          if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
            write(6,2)
            write(6,3)
            write(6,4)
            write(6,'('' '',''    There were '',i4,'' iterations on'',
     +                      '' the Gauss-Seidel algorithim '',
     +                      ''completed.'')') siter
#           ifndef Cyber
              call flush(6)
#           endif
          endif
          error = .true.
          return
        endif
c
c     Is it to be Orthomin(0), the minimum residual form?
c
      elseif ( oursol .eq. 'minres' ) then
        call minres( error, siter )
        if ( error ) then
          if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
            write(6,2)
            write(6,3)
            write(6,4)
            write(6,'('' '',''    There were '',i4,'' iterations on'',
     +                      '' the Minimum Residual algorithim '',
     +                      ''completed.'')')
     +           siter
#           ifndef Cyber
              call flush(6)
#           endif
          endif
          error = .true.
          return
        endif
#     ifdef DIRECT
c
c     Is it to be a direct LINPK solution?
c
      elseif ( oursol .eq. 'linpk' ) then
#       ifdef Cyber
        call linpk( abd, fd, sd, error, info )
#       else
        call linpk( abd, fd, sd, ipvt, error, info )
#       endif
        if ( error ) then
          write(6,'('' '',''diagonal zero at '',i6,''.'')') info
          write(6,'(10(g10.3,1x),3x,g10.3)') (J(info,siter),
     +                                        siter=1,10), F(info)
          if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
            write(6,2)
            write(6,3)
            write(6,4)
#           ifndef Cyber
              call flush(6)
#           endif
          endif
          error = .true.
          return
        endif
#     endif
c
c     Flag an error for incorrect type of solution chosen
c
      else
        write(6,'('' '',''Sorry, you have chosen a solution type '',
     +                  ''that does not exist in''/
     +            '' '',''this executable version of Shark.'')')
        error = .true.
        return
      endif
c
c     Get the solution time (we'll take out the bias later)
c
      lptime = dtime(tarray)
c
c     Now I need to go through the system again, updating the variables
c     that I have the corrections for by using Brown and Lindsay's
c     scheme
c
c               assume the max change will be in voltage
      volt = .true.
      do 300 k=1,ynodes
        do 300 i=1,xnodes
c
c         calculate where I can find the changes I have just calculated
c
          if ( ynodes .lt. xnodes ) then
            l = (i - 1)*ynodes + k
          else
            l = (k - 1)*xnodes + i
          endif
          Fnl = 2*l
          Fvl = Fnl - 1
c
c         now change delu(i) to go along with Brown and Lindsay's
c         convergence criteria for the voltage.
c
          if ( debug .eq. 'yes' ) then
            write(6,*) 'delu(',Fvl,')=',delu(Fvl)
#           ifndef Cyber
              call flush(6)
#           endif
          endif
          if ( (abs(delu(Fvl)) .gt. 1.0) .and.
     +         (abs(delu(Fvl)) .le. 3.7) ) then
            delu(Fvl) = sign( abs(delu(Fvl))**0.2, delu(Fvl) )
          elseif ( abs(delu(Fvl)) .gt. 3.7 ) then
            delu(Fvl) = sign( log(abs(delu(Fvl))), delu(Fvl) )
          endif
c
c         update the voltage matrix
c
          v(i,k) = v(i,k) + delu(Fvl)
c
c         now change delu(i) to go along with Brown and Lindsay's
c         convergence criteria for the electron quasi-fermi
c         potential.
c
          if ( debug .eq. 'yes' ) then
            write(6,*) 'delu(',Fnl,')=',delu(Fnl)
#           ifndef Cyber
              call flush(6)
#           endif
          endif
          if ( (abs(delu(Fnl)) .gt. 1.0) .and.
     +         (abs(delu(Fnl)) .le. 3.7) ) then
            delu(Fnl) = sign( abs(delu(Fnl))**0.2, delu(Fnl) )
          elseif ( abs(delu(Fnl)) .gt. 3.7 ) then
            delu(Fnl) = sign( log(abs(delu(Fnl))), delu(Fnl) )
          endif
c
c         update the electron concentration matrix only if we aren't
c         in the oxide. (there it is always zero)
c
          if ( k .le. joxide ) then
            n(i,k) = 0.0
          else
            n(i,k) = n(i,k)*exp( delu(Fvl) - delu(Fnl) )
c               limit n to 10x less than rho=0 limit
            n(i,k) = max( n(i,k), (1.0/(10.0*nio)) )
          endif
c
c         save the maximum delu(i) for convergence tests, and the
c         location and type of the maximum update
c
          if ( abs( delu(Fvl) ) .gt. norm ) then
            norm = abs( delu(Fvl) )
            inorm = i
            jnorm = k
            volt = .true.
          endif
c
c         ...for the electron quasi-fermi level, we need to 
c         make sure that the change in phin (delu(Fnl)) is
c         significant. In some cases, the change is small and
c         the carrier concentration even smaller. Ignore these.
c
          deln = delu(Fvl) - delu(Fnl)
          if ( abs(deln) .gt. 100.0*vmax ) deln = exp(deln) - 1.0
          if ( (abs( delu(Fnl) ) .gt. norm) .and.
     +         ((deln * n(i,k)) .gt. (1.0/nio)) ) then
            norm = abs( delu(Fnl) )
            inorm = i
            jnorm = k
            volt = .false.
          endif
c
c     We are all done updating the system
c
  300 continue
c
c     well, we have done another pass, so count it
c
      iter = iter + 1
c
c     get the total elapsed time for this iteration
c
      ntime = dtime(tarray)
c
c     remove the bias in the times to create split times
c
      ittime = ntime - lptime
      lptime = lptime - sutime
      sutime = sutime - time0
      ntime = ntime - sitime
c
c     do the debugging stuff if need be
c
      if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
c
c       output, for now, the statistics for this pass
c
        if ( oursol .eq. 'ortho' ) then
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' and '',i3,
     +                      '' passes of ortho('',i1,
     +                      '') with norm = '',g9.3,'' at V('',i3,
     +                      '','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           siter, keep, norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' and '',i3,
     +                      '' passes of ortho('',i1,
     +                      '') with norm = '',g9.3,'' at n('',i3,
     +                      '','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           siter, keep, norm, inorm, jnorm
          endif
        elseif ( oursol .eq. 'minres' ) then
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' and '',i3,
     +                      '' passes of minres() with norm = '',
     +                      g9.3,'' at V('',i3,'','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           siter, norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' and '',i3,
     +                      '' passes of minres() with norm = '',
     +                      g9.3,'' at n('',i3,'','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           siter, norm, inorm, jnorm
          endif
        elseif ( oursol .eq. 'gauss' ) then
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' and '',i3,
     +                      '' passes of gauss() with norm = '',
     +                      g9.3,'' at V('',i3,'','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           siter, norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' and '',i3,
     +                      '' passes of gauss() with norm = '',
     +                      g9.3,'' at n('',i3,'','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           siter, norm, inorm, jnorm
          endif
        else
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' with norm = '',
     +                      g9.3,'' at V('',i3,'','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,
     +                      ''+'',i2,'':'',i2.2,'' with norm = '',
     +                      g9.3,'' at n('',i3,'','',i3,'').'')') iter,
     +           int(sutime/60.0), int(sutime-int(sutime/60.0)*60.0),
     +           int(lptime/60.0), int(lptime-int(lptime/60.0)*60.0),
     +           int(ittime/60.0), int(ittime-int(ittime/60.0)*60.0),
     +           norm, inorm, jnorm
          endif
        endif
#       ifndef Cyber
          call flush(6)
#       endif
c
c       Every so often we need to do some maintenance work on
c       the simulation. We need to let the user know what the
c       device is "looking" like, and check to see that we
c       don't need to add more nodes to the device due to high
c       band bending. Do this every 20 iterations - at
c       least, more often if the maxit is smaller that 100.
c
        if ( ((maxit .gt. 100) .and. 
     +        (iter .eq. (iter/20)*20)) .or.
     +       ((maxit .le. 100) .and. 
     +        (iter .eq. ((iter/(maxit/5))*(maxit/5)))) 
     +     ) then
c
c         Print out the depletion map finely, for debugging purposes,
c         stopping early if we can.
c
          call depmap( 1, 1 )
c
c         Now we need to see if we need to insert nodes into the
c         simulation because we have regions in the device that
c         have too quickly changing V(), and therefore don't have
c         as good an approximation to the differentials because
c         of the non-linear equations. But only do this is we
c         aren't already close to the solution anyway...
c
c          if ( norm .gt. 10.0*vmax ) then
          if ( norm .gt. 1.0 ) then
            call refine( .false., .false., .false., error )
            if ( error ) then
              write(6,'('' '',''**** Warning! ****''/
     +                  '' '',''     The periodic mesh refinement '',
     +                        ''failed and this does NOT bode''/
     +                  '' '',''     well for the stability of the '',
     +                        ''simulation. Please be careful of ''/
     +                  '' '',''     poor convergence, or none at '',
     +                        ''all.'')')
c               reset the error flag for now
              error = .false.
            endif
          else
            write(6,'('' '',''-- Note: --''/
     +                '' '',''   While the device warranted a mesh '',
     +                      ''refinement, the solution was''/
     +                '' '',''   close enough that it is believed '',
     +                      ''that this mesh can be used''/
     +                '' '',''   to obtain the solution.'')')
          endif
        endif
      endif
c
c     check and see if this next iteration could be solved with
c     a faster, iterative, solver
c
c      if ( (oursol .eq. 'minres') .or. (oursol .eq. 'ortho') ) then
c        if ( norm .lt. 0.1 ) then
c          oursol = 'gauss'
c        endif
c      endif
c      if ( oursol .eq. 'ortho' ) then
c        if ( norm .lt. 2.0 ) then
c          oursol = 'minres'
c        endif
c      elseif ( oursol .eq. 'minres' ) then
c        if ( norm .lt. 0.1 ) then
c          oursol = 'gauss'
c        endif
c      endif
c
c     check iteration count, and add an error message, if need be
c
      if ( iter .gt. maxit ) then
        if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
          write(6,'('' '')')
          write(6,'('' '',''*** Error! The solution failed to '',
     +                    ''converge in '',i4,
     +                    '' Newton iterations!'')') maxit
          write(6,'('' '')')
c               output the entire map for debugging
          call depmap( 1, 0 )
          write(6,'('' '',''This unsuccessful solution took a '',
     +                    ''total of '',i5,'':'',i2.2,''.'')')
     +              int(ntime/60.0), int(ntime-int(ntime/60.0)*60.0)
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
        if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
          write(6,'('' '')')
          write(6,'('' '',''*** Error! The solution has an '',
     +                    ''adjustment of '',f6.2,'' kT/q!''/
     +              '' '',''    This is divergent! So the solution '',
     +                    ''is being halted!'')') norm
          write(6,'('' '',''This unsuccessful solution took a '',
     +                    ''total of '',i5,'':'',i2.2,''.'')')
     +              int(ntime/60.0), int(ntime-int(ntime/60.0)*60.0)
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
 400  write(6,'('' '',''   The Newton solution converged in '',
     +                i4,'' iterations which took '',i5,'':'',
     +                i2.2,''.'')') iter, int(ntime/60.0), 
     +                              int(ntime-int(ntime/60.0)*60.0)
#     ifndef Cyber
        call flush(6)
#     endif
c
c     output the charge distribution for this equilibrium solution
c
      if ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) then
        write(6,'('' '')')
        write(6,'('' '')')
        write(6,'('' '',''==== Shark   v2.0                      '',
     +     ''                                             REB 1988'',
     +     ''         charge map ===='')')
        write(6,'('' '')')
        write(6,'('' '','' The following is the charge map for '',
     +                  ''equilibrium:''/
     +            '' '',''    g = resistive gate''/
     +            '' '',''    i = insulator, or oxide''/
     +            '' '',''    n = n-type''/
     +            '' '',''    p = p-type''/
     +            '' '',''    blank = depletion region'')')
        write(6,'('' '')')
        write(6,'('' '')')
c
c       call the depletion mapper with a fine set, and stop early
c
        call depmap( 1, 1 )
      endif
c
c     all done... so return to the caller
c
      error = .false.
      return
      end
