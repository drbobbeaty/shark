      subroutine full1D ( iline, error, verbose, iter )
c
c     This subroutine solves the system for a consistent V
c     assuming phin, and p held in equilibrium, using Newton
c     iteration. The Jacobian is set up, and RHS is set up, and
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
c     Now for some parameters that are constants at the machine's
c     precision
c
      PRECISION one
      parameter( one=1.0 )
c
c     now for the arguments to this routine
c
      integer iline, iter
      logical error, verbose
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
c     now for the function declaration I need
c
      real dtime
c
c     Now I need to see if the computer I am on can do a direct
c     maxtrix inversion
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
c     See if we are to override the default memory model
c
#     ifdef SMALL
#       undef DIRECT
#     endif
c
c     now for some local variables
c
#ifdef DIRECT
      integer info
#endif
      integer i, k, siter, inorm, jnorm, l,
     +        LPMl, LPMr, RPMl, RPMr, SPMl, SPMr,
     +        TPMl, TPMr, BPMl, BPMr, Fvl, Fnl
      real tarray(2), ittime, ntime
      PRECISION vt, vb, vij,
     +          ht, hb, vsatt, vsatb,
     +          nt, nb, nij, norm, b, d, Fnorm,
     +          ft, fty, fb, fby,
     +          Texp, Bexp, Texpm1, Bexpm1
      PRECISION GR, dGRdV, dGRdn, velt, velb,
     +          difft, diffb, dFndV, dFndVt, dFndVb,
     +          dDtdV, dDbdV, dvtdV, dvbdV,
     +          dDtdVt, dDbdVt, dDtdVb, dDbdVb,
     +          dvtdVt, dvbdVt, dvtdVb, dvbdVb,
     +          Qt, Qb, dFndn, dFndnt, dFndnb,
     +          Jt, Jb,
     +          dJtdV, dJbdV, dJtdn, dJbdn,
     +          dJtdVt, dJbdVt, dJtdVb, dJbdVb,
     +          dJtdnt, dJbdnt, dJtdnb, dJbdnb,
     +          ayt, ayb, ut, ub, d2Vdy2, rho, expm1
      PRECISION nold(Ny), nijold, Vob, Voc, diagd, linf
      logical volt

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
c              ml      mu
c              :       :
      fd = 2*( 3 ) + ( 3 ) + 1
c               keep the size of abd() the same, but don't
c				blow ipvt() out of the water
      sd = min( ( (6*Ny+4)*(2*Nx*Ny)/fd ), 2*Nx*Ny )
#     endif
c
c     Preset some local variables so the compiler won't complain
c
      volt = .false.
c
c     precalculate the bulk potential
c
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
c     reset the number of iterations
c
      iter = 0
c
c     Copy the vector n(iline,j) into nold(j) for the time evolution,
c     and reset the updating flag for quasi-fermi levels
c
      do 50 k=1,ynodes
   50   nold(k) = n(iline,k)
c
c     reset the timer for the start of the first pass
c
      ittime = dtime(tarray)
      ntime = 0.0
c
c     Get the number exp(-1) for speed's sake, we will use it
c     a lot.
c
      expm1 = exp( -one )
c
c     Set the locations of the sub-matricies in the full form of J
c
      if ( soltyp .eq. 'linpk' ) then
        LPMl = 9
        LPMr = 10
        TPMl = 1
        TPMr = 2
        SPMl = 3
        SPMr = 4
        BPMl = 5
        BPMr = 6
        RPMl = 7
        RPMr = 8
      else
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
      endif
c
c     start a pass, and reset the norm of delta V for this pass
c
  100 norm = -1.0
      Fnorm = -1.0
c
c     For every node in the system, collect the data needed and fill
c     in [J] and -F so that we can solve them to get the corrections
c     and, in the end, the result.
c
c     Fill it all in, for this line, and reset the l infinity norm of [J]
c
      linf = 0.0
      do 200 k=1,ynodes
c
c       We are on a fixed, inputted, line and the equation count is simply
c       the node number.
c
        i = iline
        l = k
c
c       Get the location of the Fv and Fn equations
c
        Fvl = 2*l - 1
        Fnl = 2*l
c
c       collect the variables I will need for the simplified equations
c
c       get ht, hb (use symmetry at the boundaries)
c
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
c       get  vt, vb, vij (using boundary condition vectors)
c
        if ( k .eq. 1 ) then
c				only specify the gate voltage where
c				it is there, otherwise it is the
c				same as where we are.
          if ( vgate(i) .eq. 1.0e20 ) then
            vt = v(i,k)
          else
            vt = vgate(i)
          endif
        else
          vt = v(i,(k-1))
        endif

c				LOS at the bottom
        if ( k .eq. ynodes ) then
          vb = v(i,(k-1))
        else
          vb = v(i,(k+1))
        endif

        vij = v(i,k)
c
c       get nt, nb, nij and nijold
c
c				we can't have electrons above
c				the oxide
        if ( k .eq. 1 ) then
          nt = 0.0
        else
          nt = n(i,(k-1))
        endif
c				LOS at the bottom
        if ( k .eq. ynodes ) then
          nb = n(i,(k-1))
        else
          nb = n(i,(k+1))
        endif
c				we are always where we are
        nij = n(i,k)
        nijold = nold(k)
c				check to see if we are in the oxide
        if ( k .le. (joxide-1) ) then
          nb = 0.0
        endif
        if ( k .le. joxide ) then
          nij = 0.0
          nijold = 0.0
        endif
        if ( k .le. (joxide+1) ) then
          nt = 0.0
        endif
c
c       now I need to compute the coefficients for the electric
c       field functions (but with no electric field components)
c
        b = ( ht*vb - (ht+hb)*vij + hb*vt )/( ht*hb*(ht+hb) )
        d = ( ht*ht*vb + (hb*hb-ht*ht)*vij - hb*hb*vt )/
     +      ( ht*hb*(ht+hb) )
c
c       now let's compute the actual normalized electric field
c       components for each of the four points of interest
c
        fty = b*ht - d
        fby = -b*hb - d
c
c       we will also need the unit vector decomposition for
c       later decomposing of the velocity (for 1D it will be +/-1)
c
        ft = abs( fty )
        ayt = sign( one, fty )

        fb = abs( fby )
        ayb = sign( one, fby )
c				un-normalize the electric field for computing
        ft = ft * (kT/Ldi)
        fb = fb * (kT/Ldi)
c
c       now I need to compute the normalized velt, velb,
c       and the normalized difft, diffb (being careful
c       not to compute impossible results)
c
        vsatt = Av*exp(-ft/Bv) + Cv*(1.0/(1.0 + (ft/Dv)**Fv ))
        velt = mobil*ft*(1.0/sqrt(1.0+(mobil*ft/vsatt)**2)) * Ldi
        ut = mobil*(1.0/sqrt(1.0+(mobil*ft/vsatt)**2)) * kT
        if ( (ft .ne. 0.0) .and. (dDo .ne. 0.0) ) then
          difft = Do + dDo*exp( -((log(ft)-Ad)/Bd)**2 )
        else
          difft = Do
        endif
        velt = velt*ayt

        vsatb = Av*exp(-fb/Bv) + Cv*(1.0/(1.0 + (fb/Dv)**Fv ))
        velb = mobil*fb*(1.0/sqrt(1.0+(mobil*fb/vsatb)**2)) * Ldi
        ub = mobil*(1.0/sqrt(1.0+(mobil*fb/vsatb)**2)) * kT
        if ( (fb .ne. 0.0) .and. (dDo .ne. 0.0) ) then
          diffb = Do + dDo*exp( -((log(fb)-Ad)/Bd)**2 )
        else
          diffb = Do
        endif
        velb = velb*ayb
c
c       Now compute some exponentials that we are going to be using
c       quite a lot
c
        Texpm1 = ht*velt/difft
        Texp = exp( Texpm1 )
        if ( abs(Texpm1) .gt. 100.0*vmax ) Texpm1 = Texp - 1.0
c				check for the divide by zero point
        if ( velt .eq. 0.0 ) then
          Texp = 1.0
          Texpm1 = 1.0
        endif

        Bexpm1 = hb*velb/diffb
        Bexp = exp( Bexpm1 )
        if ( abs(Bexpm1) .gt. 100.0*vmax ) Bexpm1 = Bexp - 1.0
c				check for the divide by zero point
        if ( velb .eq. 0.0 ) then
          Bexp = 1.0
          Bexpm1 = 1.0
        endif
c
c       These are the partials of the Generation-Recombination
c       (SRH) with respect to the variables
c
        GR = 0.0
        dGRdV = 0.0
        dGRdn = 0.0
c
c       now compute the currents for this node, depending on the
c       location of the node
c
        if ( k .le. joxide ) then
c				we are in the oxide
          Jt = 0.0
          Jb = 0.0
        else
c				we are NOT in the oxide
          Jb = -velb*( Bexp*nb - nij )/Bexpm1
c				Jt into the oxide is zero
          if ( k .eq. (joxide+1) ) then
            Jt = 0.0
          else
            Jt = -velt*( Texp*nij - nt )/Texpm1
          endif
        endif
c
c       now the tough diffusion partials, but be careful not to
c       try to compute log(0)!
c
        if ( (ft .ne. 0.0) .and. (abs(Do) .gt. 0.0) ) then
          dDtdV = ( -(2.0/(Bd*Bd))*(log(ft)-Ad)*(difft-Do)/ft )
     +            *(kT/Ldi)
          dDtdVt = dDtdV * ( ayt/ht )
          dDtdVb = 0.0
          dDtdV = -dDtdVt
        else
          dDtdV = 0.0
          dDtdVt = 0.0
          dDtdVb = 0.0
        endif

        if ( (fb .ne. 0.0) .and. (abs(Do) .gt. 0.0) ) then
          dDbdV = ( -(2.0/(Bd*Bd))*(log(fb)-Ad)*(diffb-Do)/fb )
     +            *(kT/Ldi)
          dDbdVt = 0.0
          dDbdVb = dDbdV * ( -ayb/hb )
          dDbdV = -dDbdVb
        else
          dDbdV = 0.0
          dDbdVt = 0.0
          dDbdVb = 0.0
        endif
c
c       now the tough velocity partials
c
        if ( ft .ne. 0.0 ) then
          Qt = sqrt( 1.0 + ( mobil*ft/vsatt )**2 )
          if ( ( (abs(Bv) .gt. 1.0e30) .and. (Cv .eq. 0.0) ) .or.
     +         ( (abs(Dv) .gt. 1.0e30) .and. (Av .eq. 0.0) ) ) then
            dvtdV = ( (mobil/Qt)*( 1.0 - (mobil*ft/(Qt*vsatt))**2 )
     +              )*kT*ayt
          else
            dvtdV = ( (mobil/(Qt*Qt))*( Qt - (mobil*mobil*ft*ft/(Qt*
     +              vsatt*vsatt*vsatt))*( vsatt - ft*(-(Av/Bv)*exp(-ft
     +              /Bv) - (Cv*Fv/Dv)/((ft/Dv)**(1.0-Fv) * (1.0 + (ft/
     +              Dv)**Fv)**2 ) ) ) ) )*kT*ayt
          endif
          dvtdVt = dvtdV * ( ayt/ht )
          dvtdVb = 0.0
          dvtdV = -dvtdVt
        else
          dvtdV = mobil*kT*ayt
          dvtdVt = dvtdV * ( ayt/ht )
c          dvtdVt = dvtdV * ( -(1.0-ht)/ht )*1.41421356237/
c     +             sqrt( 1.0 + (1.0-ht)**2 )
          dvtdVb = 0.0
          dvtdV = -dvtdVt
        endif

        if ( fb .ne. 0.0 ) then
          Qb = sqrt( 1.0 + ( mobil*fb/vsatb )**2 )
          if ( ( (abs(Bv) .gt. 1.0e30) .and. (Cv .eq. 0.0) ) .or.
     +         ( (abs(Dv) .gt. 1.0e30) .and. (Av .eq. 0.0) ) ) then
            dvbdV = ( (mobil/Qb)*( 1.0 - (mobil*fb/(Qb*vsatb))**2 )
     +              )*kT*ayb
          else
            dvbdV = ( (mobil/(Qb*Qb))*( Qb - (mobil*mobil*fb*fb/(Qb*
     +              vsatb*vsatb*vsatb))*( vsatb - fb*(-(Av/Bv)*exp(-fb
     +              /Bv) - (Cv*Fv/Dv)/((fb/Dv)**(1.0-Fv) * (1.0 + (fb/
     +              Dv)**Fv)**2 ) ) ) ) )*kT*ayb
          endif
          dvbdVt = 0.0
          dvbdVb = dvbdV * ( -ayb/hb )
          dvbdV = -dvbdVb
        else
          dvbdV = mobil*kT*ayb
          dvbdVt = 0.0
          dvbdVb = dvbdV * ( -ayb/hb )
c          dvbdVb = dvbdV * ( -(1.0-hb)/hb )*1.41421356237/
c     +             sqrt( 1.0 + (1.0-hb)**2 )
          dvbdV = -dvbdVb
        endif
c
c       now the current partials, being careful to use the proper
c       form for the region we are in
c
        if ( k .le. joxide ) then
c				we are in the insulator
          dJtdV = 0.0
          dJtdVt = 0.0
          dJtdVb = 0.0
          dJtdn = 0.0
          dJtdnt = 0.0
          dJtdnb = 0.0

          dJbdV = 0.0
          dJbdVt = 0.0
          dJbdVb = 0.0
          dJbdn = 0.0
          dJbdnt = 0.0
          dJbdnb = 0.0
        else
c				we are NOT in the insulator
          if ( velb .ne. 0.0 ) then
            dJbdV = (velb*hb/diffb)*(Bexp/Bexpm1**2)*(nb - nij)*
     +              ( dvbdV - (velb/diffb)*dDbdV ) +
     +              dvbdV*( nij - Bexp*nb )/Bexpm1
            dJbdVt = (velb*hb/diffb)*(Bexp/Bexpm1**2)*(nb - nij)*
     +               ( dvbdVt - (velb/diffb)*dDbdVt ) +
     +               dvbdVt*( nij - Bexp*nb )/Bexpm1
            dJbdVb = (velb*hb/diffb)*(Bexp/Bexpm1**2)*(nb - nij)*
     +               ( dvbdVb - (velb/diffb)*dDbdVb ) +
     +               dvbdVb*( nij - Bexp*nb )/Bexpm1
            dJbdn = velb/Bexpm1
            dJbdnt = 0.0
            dJbdnb = -velb*(Bexp/Bexpm1)
          else
            dJbdV = -dvbdV * nij
            dJbdVt = -dvbdVt * nij
            dJbdVb = -dvbdVb * nij
            dJbdn = diffb/hb
            dJbdnt = 0.0
            dJbdnb = -diffb/hb - velb
          endif
c
c         check to see if the center node is in the insulator,
c         so don't let current go streaming into the insulator
c
          if ( k .eq. (joxide+1) ) then
            dJtdV = 0.0
            dJtdVt = 0.0
            dJtdVb = 0.0
            if ( velt .ne. 0.0 ) then
              dJtdn = -velt*(Texp/Texpm1)
              dJtdnt = velt/Texpm1
            else
              dJtdn = -difft/ht - velt
              dJtdnt = difft/ht
            endif
            dJtdnb = 0.0
          else
            if ( velt .ne. 0.0 ) then
              dJtdV = (velt*ht/difft)*(Texp/Texpm1**2)*(nij - nt)*
     +                ( dvtdV - (velt/difft)*dDtdV ) +
     +                dvtdV*( nt - Texp*nij )/Texpm1
              dJtdVt = (velt*ht/difft)*(Texp/Texpm1**2)*(nij - nt)*
     +                 ( dvtdVt - (velt/difft)*dDtdVt ) +
     +                 dvtdVt*( nt - Texp*nij )/Texpm1
              dJtdVb = (velt*ht/difft)*(Texp/Texpm1**2)*(nij - nt)*
     +                 ( dvtdVb - (velt/difft)*dDtdVb ) +
     +                 dvtdVb*( nt - Texp*nij )/Texpm1
              dJtdn = -velt*(Texp/Texpm1)
              dJtdnt = velt/Texpm1
              dJtdnb = 0.0
            else
              dJtdV = -dvtdV * nij
              dJtdVt = -dvtdVt * nij
              dJtdVb = -dvtdVb * nij
              dJtdn = -difft/ht - velt
              dJtdnt = difft/ht
              dJtdnb = 0.0
            endif
          endif
        endif
c
c       now the main partials themselves
c
        dFndV = (2.0/(ht+hb))*(dJbdV - dJtdV) + dGRdV
        dFndVt = (2.0/(ht+hb))*(dJbdVt - dJtdVt)
        dFndVb = (2.0/(ht+hb))*(dJbdVb - dJtdVb)
        dFndn = (2.0/(ht+hb))*(dJbdn - dJtdn) + dGRdn
        dFndnt = (2.0/(ht+hb))*(dJbdnt - dJtdnt)
        dFndnb = (2.0/(ht+hb))*(dJbdnb - dJtdnb)
c
c       Now I can fill in the TPM sub-matrix of [J]
c
c				upper left corner element
        if ( k .eq. joxide ) then
          J( Fvl, TPMl ) = 2.0/(ht*hb*(1.0 + (ksemi/kox)*(ht/hb)))
        else
          J( Fvl, TPMl ) = 2.0/(ht*(ht+hb))
        endif
c				upper right corner element
        J( Fvl, TPMr ) = 0.0
c				lower left corner element
        J( Fnl, TPMl ) = dFndVt
c				lower right corner element
        J( Fnl, TPMr ) = -nt*dFndnt
c
c       Now I can fill in the SPM sub-matrix of [J]
c
c				upper left corner element
        if ( k .lt. joxide ) then
          J( Fvl, SPMl ) = -2.0/(ht*hb)
        elseif ( k .eq. joxide ) then
          J( Fvl, SPMl ) = -2.0/(ht*hb)
        else
          if ( symtyp .eq. 'Mock' ) then
            J( Fvl, SPMl ) = -2.0/(ht*hb)
     +                       - exp(-vij - log(dopsrc/nio)) - nij
          else
            J( Fvl, SPMl ) = -2.0/(ht*hb) - exp(-vij) - nij
          endif
        endif
c				upper right corner element
        J( Fvl, SPMr ) = nij
c				lower left corner element
        J( Fnl, SPMl ) = dFndV
c				lower right corner element
        if ( k .le. joxide ) then
          J( Fnl, SPMr ) = 1.0
        else
c          J( Fnl, SPMr ) = -nij*(dFndn + Ldi*Ldi/dt)
          J( Fnl, SPMr ) = -nij*dFndn
        endif
c
c       Now I can fill in the BPM sub-matrix of [J]
c
c				upper left corner element
        if ( k .eq. joxide ) then
          J( Fvl, BPMl ) = 2.0/(hb*(ht + hb))*(1.0+(ht/hb))/
     +                     ((kox/ksemi) + (ht/hb))
        else
          J( Fvl, BPMl ) = 2.0/(hb*(ht + hb))
        endif
c				upper right corner element
        J( Fvl, BPMr ) = 0.0
c				lower left corner element
        J( Fnl, BPMl ) = dFndVb
c				lower right corner element
        J( Fnl, BPMr ) = -nb*dFndnb
c
c       Zero out the right and left partials because they are
c       not going to be used
c
        J( Fvl, RPMl ) = 0.0
        J( Fvl, RPMr ) = 0.0
        J( Fnl, RPMl ) = 0.0
        J( Fnl, RPMr ) = 0.0

        J( Fvl, LPMl ) = 0.0
        J( Fvl, LPMr ) = 0.0
        J( Fnl, LPMl ) = 0.0
        J( Fnl, LPMr ) = 0.0
c
c       Now that all the partials are calculated we need to fix
c       them for the presence of L.O.S. boundary conditions, and
c       imaginary points from the gate and ohmic contact
c
        if ( k .eq. 1 ) then
c				zero out TPM
          J( Fvl, TPMl ) = 0.0
          J( Fvl, TPMr ) = 0.0
          J( Fnl, TPMl ) = 0.0
          J( Fnl, TPMr ) = 0.0
        endif
        if ( k .eq. ynodes ) then
c				add the partials for LOS
          J( Fvl, TPMl ) = J( Fvl, TPMl ) + J( Fvl, BPMl )
          J( Fvl, TPMr ) = J( Fvl, TPMr ) + J( Fvl, BPMr )
          J( Fnl, TPMl ) = J( Fnl, TPMl ) + J( Fnl, BPMl )
          J( Fnl, TPMr ) = J( Fnl, TPMr ) + J( Fnl, BPMr )
c				zero out BPM
          J( Fvl, BPMl ) = 0.0
          J( Fvl, BPMr ) = 0.0
          J( Fnl, BPMl ) = 0.0
          J( Fnl, BPMr ) = 0.0
        endif
c
c       Now I can fill in the section of -F for this node
c
        if ( k .lt. joxide ) then
c				-Fv terms for this node in the oxide
          d2Vdy2 = (ht*vb - (ht+hb)*vij + hb*vt)/(0.5*ht*hb*(ht+hb))
          rho = 0.0
c				-Fn terms for this node in the oxide
          F( Fnl ) = 0.0
        elseif ( k .eq. joxide ) then
c				-Fv terms for this node at the interface
          d2Vdy2 = (2.0/(ht*hb*(ht+hb)))*( vt*(ht+hb)/(1.0 + 
     +             (ksemi/kox)*(ht/hb)) - vij*(ht+hb) + vb*ht*
     +             (1.0+ht/hb)/(kox/ksemi + ht/hb) )
          rho = 2.0 * Qss(i)/(nio * hb * Ldi)
c				-Fn terms for this node at the interface
          F( Fnl ) = 0.0
        else
c				-Fv terms for this node in the semi.
          d2Vdy2 = (ht*vb - (ht+hb)*vij + hb*vt)/(0.5*ht*hb*(ht+hb))
          if ( symtyp .eq. 'Mock' ) then
            rho = exp(-vij - log(dopsrc/nio)) + dop(i,k) - nij
          else
            rho = exp(-vij) + dop(i,k) - nij
          endif
c				-Fn for this node in the semi
c          F( Fnl ) = -( 2.0*(Jb-Jt)/(hb+ht) + 
c     +               GR + (nij - nijold)*Ldi*Ldi/dt )
          F( Fnl ) = -( 2.0*(Jb-Jt)/(hb+ht) + GR )
        endif
c
c       See if Fv is zero for our purposes
c
        if ( abs(rho) .gt. (1.0/nio) ) then
          if ( abs( ( d2Vdy2 - rho )/rho ) .lt. (vmax/10.0) ) then
            F( Fvl ) = 0.0
          else
            F( Fvl ) = -( d2Vdy2 + rho )
          endif
        else
          if ( abs( d2Vdy2 ) .lt. (1.0/nio) ) then
            F( Fvl ) = 0.0
          else
            F( Fvl ) = -( d2Vdy2 )
          endif
        endif
c
c       see if the Fn's are small enough to be considered solved
c
        if ( nij .lt. (1.0/nio) ) F( Fnl ) = 0.0
c
c       see if the entire row of Fn is zero due to a depletion
c       region of electrons
c
        if ( ( J( Fnl, SPMr ) .eq. 0.0 ) .and.
     +       ( F( Fnl ) .eq. 0.0 ) ) then
          J( Fnl, SPMr ) = 1.0
        endif
c
c       Check the norm of F to see if it is small enough to be
c       considered solved without having to go through the solver.
c
        Fnorm = max( Fnorm, abs(F(Fvl)), abs(F(Fnl)) )
c
c       see if the Fn diagonal is so small as to be considered
c       ineffective in the solution
c
c        if ( abs( J( Fnl, SPMr ) ) .lt. (vmax*1.0e-6) ) then
c          J( Fnl, TPMl ) = 0.0
c          J( Fnl, TPMr ) = 0.0
c          J( Fnl, SPMl ) = 0.0
c          J( Fnl, SPMr ) = 0.0
c          J( Fnl, BPMl ) = 0.0
c          J( Fnl, BPMr ) = 0.0
c          F( Fnl ) = 0.0
c        endif
c
c       now normalize the rows of [J] w.r.t. the diagonal
c       elements and check for diagonal dominance, and l infinity
c
        diagd = 0.0
        if ( J( Fvl, SPMl ) .ne. 0.0 ) then
          b = 1.0/J( Fvl, SPMl )
          J( Fvl, TPMl ) = J( Fvl, TPMl ) * b
          J( Fvl, TPMr ) = J( Fvl, TPMr ) * b
          J( Fvl, SPMl ) = 1.0
          J( Fvl, SPMr ) = J( Fvl, SPMr ) * b
          J( Fvl, BPMl ) = J( Fvl, BPMl ) * b
          J( Fvl, BPMr ) = J( Fvl, BPMr ) * b
          F( Fvl ) = F( Fvl ) * b
          diagd = max( abs(J(Fvl,TPMl)), abs(J(Fvl,TPMr)),
     +                 abs(J(Fvl,SPMr)), abs(J(Fvl,BPMl)),
     +                 abs(J(Fvl,BPMr)) )
          if ( diagd .lt. 1.0 ) diagd = 0.0
          linf = max ( linf, diagd )
        else
          write(6,'('' '',''F('',i5,'')=Fv('',i3,'','',i3,
     +          '') has a zero diagonal.'')') Fvl, i, k
          if ( ( debug .eq. 'no' ) .and. ( verbose ) ) then
            write(6,'(10(g10.3,1x),3x,g10.3)') (J(Fvl,siter),
     +          siter=1,10), F(Fvl)
          endif
        endif
c
c       check for diagonal dominance
c
        if ( ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) .and.
     +       ( diagd .ne. 0.0 ) .and. ( soltyp .ne. 'linpk' ) .and.
     +       ( verbose ) ) then
          write(6,'('' '',''Trouble! The system matrix is non-'',
     +         ''diagonally dominant. An off diagonal of Fv''/
     +         '' '',''is '',g9.3,'' times the diagonal'',
     +         '' at node '',i3,'','',i3,''.'')') diagd, i, k
        endif

        diagd = 0.0
        if ( J( Fnl, SPMr ) .ne. 0.0 ) then
          b = 1.0/J( Fnl, SPMr )
          J( Fnl, TPMl ) = J( Fnl, TPMl ) * b
          J( Fnl, TPMr ) = J( Fnl, TPMr ) * b
          J( Fnl, SPMl ) = J( Fnl, SPMl ) * b
          J( Fnl, SPMr ) = 1.0
          J( Fnl, BPMl ) = J( Fnl, BPMl ) * b
          J( Fnl, BPMr ) = J( Fnl, BPMr ) * b
          F( Fnl ) = F( Fnl ) * b
          diagd = max( abs(J(Fnl,TPMl)), abs(J(Fnl,TPMr)),
     +                 abs(J(Fnl,SPMl)), abs(J(Fnl,BPMl)),
     +                 abs(J(Fnl,BPMr)) )
          if ( diagd .lt. 1.0 ) diagd = 0.0
          linf = max ( linf, diagd )
        else
          write(6,'('' '',''F('',i5,'')=Fn('',i3,'','',i3,
     +          '') has a zero diagonal.'')') Fnl, i, k
          if ( ( debug .eq. 'no' ) .and. ( verbose ) ) then
            write(6,'(10(g10.3,1x),3x,g10.3)') (J(Fnl,siter),
     +          siter=1,10), F(Fnl)
          endif
        endif
c
c       check for diagonal dominance
c
        if ( ( (debug .eq. 'yes') .or. (debug .eq. 'sortof') ) .and.
     +       ( diagd .ne. 0.0 ) .and. ( soltyp .ne. 'linpk' ) .and.
     +       ( verbose ) ) then
          write(6,'('' '',''Trouble! The system matrix is non-'',
     +         ''diagonally dominant. An off diagonal of Fn''/
     +         '' '',''is '',g9.3,'' times the diagonal'',
     +         '' at node '',i3,'','',i3,''.'')') diagd, i, k
        endif
c
c       All done filling out [J] and -F!
c
        if ( ( debug .eq. 'yes' ) .and. ( verbose ) ) then
          if ( soltyp .eq. 'linpk' ) then
            write(6,'(6(g10.3,1x),3x,g10.3)') (J(Fvl,siter),
     +          siter=1,6), F(Fvl)
            write(6,'(6(g10.3,1x),3x,g10.3)') (J(Fnl,siter),
     +          siter=1,6), F(Fnl)
          else
            write(6,'(6(g10.3,1x),3x,g10.3)') (J(Fvl,siter),
     +          siter=3,8), F(Fvl)
            write(6,'(6(g10.3,1x),3x,g10.3)') (J(Fnl,siter),
     +          siter=3,8), F(Fnl)
          endif
#         ifndef Cyber
            call flush(6)
#         endif
        endif
  200 continue
c
c     See if we even need to call the solver based on the L-infinity
c     norm of the F vector
c
      if ( Fnorm .lt. (vmax/1.0e3) ) then
        if ( ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) .and.
     +     ( verbose ) ) then
          write(6,'('' '',''   The L-inf norm of F is '',g9.3,
     +         '' and is therefore considered solved.'')') Fnorm
        endif
        goto 400
      endif
c
c     Print out the l infinity norm if debug is anything BUT 'no'
c
      if ( ( debug .ne. 'no' ) .and. ( linf .gt. 0.0 ) .and.
     +     ( soltyp .ne. 'linpk' ) .and. ( verbose ) ) then
        write(6,'('' '',''The L-inf norm of [J] is '',g10.3,
     +        '' which may cause convergence trouble.'')') linf
      endif
c
c     Now we can call the solution system for this mess, depending on
c     how we want it solved
c
c     Is it to be Orthomin(k)?
c
      if ( soltyp .eq. 'ortho' ) then
        call ortho( error, siter )
        if ( error ) then
          if ( ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) .and.
     +         ( verbose ) ) then
            write(6,'('' '',''*** Error! ***''/
     +                '' '',''    While trying to solve the system of '',
     +                      ''equations, Shark has detected an error''/
     +                '' '',''    in the solution itself. This has '',
     +                      ''made it impossible to solve the system.''/
     +                '' '',''    There were '',i4,'' iterations on '',
     +                      ''the Orthomin('',i1,'') algorithim ''/
     +                '' '',''    completed.'')') siter, keep
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
      elseif ( soltyp .eq. 'gauss' ) then
        call gauss( error, siter )
        if ( error ) then
          if ( ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) .and.
     +         ( verbose ) ) then
            write(6,'('' '',''*** Error! ***''/
     +                '' '',''    While trying to solve the system of'',
     +                      '' equations, Shark has detected an error''/
     +                '' '',''    in the solution itself. This has '',
     +                      '' made it impossible to solve the''/
     +                '' '',''    system. There were '',i4,
     +                      '' iterations on the Gauss-Seidel '',
     +                      ''algorithim completed.'')') siter
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
      elseif ( soltyp .eq. 'minres' ) then
        call minres( error, siter )
        if ( error ) then
          if ( ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) .and.
     +         ( verbose ) ) then
            write(6,'('' '',''*** Error! ***''/
     +                '' '',''    While trying to solve the system of'',
     +                      '' equations, Shark has detected an error''/
     +                '' '',''    in the solution itself. This has '',
     +                      '' made it impossible to solve the''/
     +                '' '',''    system. There were '',i4,
     +                      '' iterations on the Minimum Residual '',
     +                      ''algorithim completed.'')') siter
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
      elseif ( soltyp .eq. 'linpk' ) then
#       ifdef Cyber
        call Lpk1D( abd, fd, sd, error, info )
#       else
        call Lpk1D( abd, fd, sd, ipvt, error, info )
#       endif
        if ( error ) then
          write(6,'('' '',''diagonal zero at '',i6,''.'')') info
          if ( verbose ) then
            write(6,'(10(g10.3,1x),3x,g10.3)') (J(info,siter),
     +             siter=1,10), F(info)
          endif
          if ( ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) .and.
     +         ( verbose ) ) then
            write(6,'('' '',''*** Error! ***''/
     +                '' '',''    While trying to solve the system of '',
     +                      ''equations, Shark has detected an error''/
     +                '' '',''    in the solution itself. This has '',
     +                      ''made it impossible to solve the '',
     +                      ''system.'')')
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
        if ( verbose ) then
          write(6,'('' '',''Sorry, you have chosen a solution type'',
     +                    '' that does not exist in''/
     +              '' '',''this executable version of Shark.'')')
        endif
        error = .true.
        return
      endif
c
c     Now I need to go through the system again, updating the variables
c     that I have the corrections for by using Brown and Lindsay's
c     scheme
c
      do 300 k=1,ynodes
c
c       We are on a fixed, inputted, line and the equation count is simply
c       the node number
c
        i = iline
        l = k
c
c       Get the location of the Fv and Fn equations
c
        Fvl = 2*l - 1
        Fnl = 2*l
c
c       now change delu(i) to go along with Brown and Lindsay's
c       convergence criteria for the voltage
c
        if ( (debug .eq. 'yes') .and. ( verbose ) ) then
          write(6,*) 'delu(',Fvl,')=',delu(Fvl)
#         ifndef Cyber
            call flush(6)
#         endif
        endif
        if ( (abs(delu(Fvl)) .gt. 1.0) .and.
     +       (abs(delu(Fvl)) .le. 3.7) ) then
          delu(Fvl) = sign( abs(delu(Fvl))**0.2, delu(Fvl) )
        elseif ( abs(delu(Fvl)) .gt. 3.7 ) then
          delu(Fvl) = sign( log(abs(delu(Fvl))), delu(Fvl) )
        endif
c
c       update the voltage matrix
c
        v(i,k) = v(i,k) + delu(Fvl)
c
c       now change delu(i) to go along with Brown and Lindsay's
c       convergence criteria for the electron quasi-fermi potential.
c
        if ( (debug .eq. 'yes') .and. ( verbose ) ) then
          write(6,*) 'delu(',Fnl,')=',delu(Fnl)
#         ifndef Cyber
            call flush(6)
#         endif
        endif
        if ( (abs(delu(Fnl)) .gt. 1.0) .and.
     +       (abs(delu(Fnl)) .le. 3.7) ) then
          delu(Fnl) = sign( abs(delu(Fnl))**0.2, delu(Fnl) )
        elseif ( abs(delu(Fnl)) .gt. 3.7 ) then
          delu(Fnl) = sign( log(abs(delu(Fnl))), delu(Fnl) )
        endif
c
c       update the electron concentration matrix only if we aren't
c       in the oxide. (there it is always zero)
c
        if ( k .le. joxide ) then
          n(i,k) = 0.0
        else
          n(i,k) = n(i,k)*exp( delu(Fvl) - delu(Fnl) )
c				limit n to 10x less than rho=0 limit
          n(i,k) = max( n(i,k), (1.0/(10.0*nio)) )
        endif
c
c       save the maximum delu(i) for convergence tests, and the
c       location and type of the maximum update
c
        if ( abs( delu(Fvl) ) .gt. norm ) then
          norm = abs( delu(Fvl) )
          inorm = i
          jnorm = k
          volt = .true.
        endif
        if ( abs( delu(Fnl) ) .gt. norm ) then
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
c     get the time for this iteration
c
      ittime = dtime(tarray)
      ntime = ntime + ittime
c
c     do the debugging stuff if need be
c
      if ( ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) .and.
     +     ( verbose ) ) then
c
c       output, for now, the statistics for this pass
c
        if ( soltyp .eq. 'ortho' ) then
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +         '' and '',i3,
     +         '' passes of ortho('',i1,'') with norm = '',g9.3,
     +         '' at V('',i3,'','',i3,'').'')') iter, int(ittime/
     +         60.0), int(ittime-int(ittime/60.0)*60.0), siter, keep,
     +         norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +         '' and '',i3,
     +         '' passes of ortho('',i1,'') with norm = '',g9.3,
     +         '' at n('',i3,'','',i3,'').'')') iter, int(ittime/
     +         60.0), int(ittime-int(ittime/60.0)*60.0), siter, keep,
     +         norm, inorm, jnorm
          endif
        elseif ( soltyp .eq. 'minres' ) then
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +         '' and '',i3,
     +         '' passes of minres() with norm = '',g9.3,
     +         '' at V('',i3,'','',i3,'').'')') iter, int(ittime/
     +         60.0), int(ittime-int(ittime/60.0)*60.0), siter,
     +         norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +         '' and '',i3,
     +         '' passes of minres() with norm = '',g9.3,
     +         '' at n('',i3,'','',i3,'').'')') iter, int(ittime/
     +         60.0), int(ittime-int(ittime/60.0)*60.0), siter,
     +         norm, inorm, jnorm
          endif
        elseif ( soltyp .eq. 'gauss' ) then
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +         '' and '',i3,
     +         '' passes of gauss() with norm = '',g9.3,
     +         '' at V('',i3,'','',i3,'').'')') iter, int(ittime/
     +         60.0), int(ittime-int(ittime/60.0)*60.0), siter,
     +         norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',i2,'':'',i2.2,
     +         '' and '',i3,
     +         '' passes of gauss() with norm = '',g9.3,
     +         '' at n('',i3,'','',i3,'').'')') iter, int(ittime/
     +         60.0), int(ittime-int(ittime/60.0)*60.0), siter,
     +         norm, inorm, jnorm
          endif
        else
          if ( volt ) then
            write(6,'('' '',''Iter. #'',i3,'' took '',f6.1,
     +         '' msec. with norm = '',g9.3,
     +         '' at V('',i3,'','',i3,'').'')') iter, (ittime*1000.0),
     +         norm, inorm, jnorm
          else
            write(6,'('' '',''Iter. #'',i3,'' took '',f6.1,
     +         '' msec. with norm = '',g9.3,
     +         '' at n('',i3,'','',i3,'').'')') iter, (ittime*1000.0),
     +         norm, inorm, jnorm
          endif
        endif
#       ifndef Cyber
          call flush(6)
#       endif
      endif
c
c     check and see if this next iteration could be solved with 
c     a faster, iterative, solver
c
c      if ( (soltyp .eq. 'minres') .or. (soltyp .eq. 'ortho') ) then
c        if ( norm .lt. 0.1 ) then
c          soltyp = 'gauss'
c        endif
c      endif
c      if ( soltyp .eq. 'ortho' ) then
c        if ( norm .lt. 2.0 ) then
c          soltyp = 'minres'
c        endif
c      endif
c      if ( soltyp .eq. 'minres' ) then
c        if ( norm .lt. 0.1 ) then
c          soltyp = 'gauss'
c        endif
c      endif
c
c     check iteration count and add an error message, if need be
c
      if ( iter .gt. (2*maxit) ) then
        if ( ( (debug .eq. 'sortof') .or. (debug .eq. 'yes') ) .and.
     +       ( verbose ) ) then
          write(6,'('' '')')
          write(6,'('' '',''*** Error! The solution failed to '',
     +       '' converge in '',i4,'' Newton iterations!'')') maxit
          write(6,'('' '')')
          write(6,'('' '',''This unsuccessful solution took a '',
     +       ''total of '',f5.1,'' sec.'')') ntime
#         ifndef Cyber
            call flush(6)
#         endif
        endif
        error = .true.
        return
      endif
      if ( norm .gt. vmax ) goto 100
c
c     write a success message
c
  400 if ( verbose ) then
        write(6,'('' '',''   The Newton solution converged in'',
     +        i4,'' iterations which took '',f5.1,'' sec.'')')
     +        iter, ntime
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
