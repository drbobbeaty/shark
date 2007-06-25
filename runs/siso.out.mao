 *******************************************************************************************************************
 *                                                                                                                 *
 *         SSSSSSSSS    h                                k                        2222222        000000000         *
 *        S             h           aaaaaa               k                       2       2      0         0        *
 *        S             h                 a   r rrrrr    k   kkk                         2      0         0        *
 *         SSSSSSSSS    h hhhhh     aaaaa a   rr     r   k kkk          v     v   2222222       0         0        *
 *                  S   hh     h   a     aa   r          kkk             v   v   2              0         0        *
 *                  S   h      h   a      a   r          k  kkk           v v    2          ..  0         0        *
 *         SSSSSSSSS    h      h    aaaaaaa   r          k    kkk          v     222222222  ..   000000000         *
 *                                                                                                                 *
 *                                                                                                                 *
 *            T W O   D I M E N S I O N A L   T R A N S I E N T   G a A s   F E T   S I M U L A T O R              *
 *                                                                                                                 *
 *                  Written by Robert E. Beaty as a part of a Ph.D. in Electrical Engineering at                   *
 *                           Purdue University, West Lafayette, Indiana 47906  June 1988                           *
 *                                                                                                                 *
 *******************************************************************************************************************
                                                                                                                    
                                                                                                                    
                                                                                                                    
 
 
 ==== Shark   v2.0                                                                  REB 1988         input deck ====
 
 *
 * This will be a comparative test of Sutherland's work for a Silicon
 * MOSFET with SiO2, biased into strong inversion, saturation mode
 *
 source thick=0.35
 drain thick=0.35
 bulk doping=2.0e16
 oxide type=SiO2 thick=1000
 channel type=MODFET length=3.0 material=Si
 * setting dDo=0 will make D constant
 DvsE Do=18.13 dDo=0.0
 * setting A=C=inf. will make vsat=inf. and then v=u*E
 VvsE mobility=700.0 A=2.0e38 C=2.0e38
 * bias the device in the triode region
 bias vsource=0.00 vdrain=4.00 eavg=0.0 vg1=4.0 delV=1.0
 * do the steady state operation: delt -> infinity
 time tout=10.0 tstop=100.0 delt=1.0e38
 solution type=linpk
 start with=DC maxit=300 debug=sortof guess=MOS
 
 
 ==== Shark   v2.0                                                                  REB 1988     run parameters ====
 
    The length of the simulated channel is   3.00 (um).
    The thickness of the simulated device is   5.10 (um).
    The doping density of the substrate, given as (Na-Nd), is 0.20E+17 (/cm**3).
 
    The insulator is SiO2   with a thickness of 1000.0 (Angst.).
 
    The dielectric constants for the materials are:
         semiconductor = 11.70
         insulator =  3.90
 
    The channel is based on the MODFET channel.
    The channel material type is Si   .
 
    The biasing on the terminals of the device is:
         source voltage =  0.000 (V).
         drain voltage =  4.000 (V).
         source end of the gate =  4.000 (V).
         drain end of the gate =  4.000 (V).
         field in the gate =     0.0 (V/cm).
 
    The parameters for the velocity versus field equation are:
          low-field mobility =   700.0 (cm**2/V-s).
                           A =  0.2000E+39
                           B =  0.1000E+39
                           C =  0.2000E+39
                           D =   1.000    
                           F =   1.000    
 
    The parameters for the diffusivity versus field equation are:
       low-field diffusivity =   18.13 (cm**2/sec).
            Do-to-peak value =    0.00 (cm**2/sec).
                           A =  8.130
                           B =  0.5960
 
    Starting the channel with the steady state electron concentration in the channel.
 
    The Newton iterations will halt if the largest change in the voltage is less than 0.10E-05 kT/q,
    or  300 iterations have been completed.
    The systems of equations will be solved by the small LINPACK algorithm.
 
    The starting time of the simulation is   0.0 (psec).
    The stopping time of the simulation is 100.0 (psec).
    Output will be generated every  10.0 (psec) of the simulation.
 
 
 
 
 ==== Shark   v2.0                                                                  REB 1988     mesh selection ====
 
 The simple MOS electrostatics for the channel region are:
       Vgmax =  4.000 V
  Vg+Qss/Cox =  4.000 V
         Cox = 0.345E-07 F/cm**2
          Vo = 2.7840 V
          V1 = 0.0000 V
        phis =  1.305 V
          xd =  0.290 um
 
 Because Shark is selecting the mesh, Shark is setting the thickness of the
 substrate to  5.124 um
 
    automatically selecting the y mesh for a MODFET channel.
 
    placing y nodes uniformly in tox...
    placing y nodes quadratically in xd...
    placing y nodes uniformly past xd to tsub...
 
    placing x nodes non-uniformly in the channel...
    setting xnodes =  32
    setting ynodes =  66
    maximum hx ratio =  1.111 from:
             hx max. = 0.1962 (um)
             hx min. = 0.0404 (um)
    maximum hy ratio =  4.654 from:
             hy max. = 0.0529 (um)
             hy min. = 0.0036 (um)
 
 
 ==== Shark   v2.0                                                                  REB 1988  setting up device ====
 
   setting up the oxide...
   setting up the channel...
   setting the substrate...
 
 
 ==== Shark   v2.0                                                                  REB 1988    applying biases ====
 
 
 
 ==== Shark   v2.0                                                                  REB 1988      initial guess ====
 
     Computing the initial guess based on MOS electrostatics for the
     center section device simulation.
     Using Pao-Sah to calculate the initial quasi-fermi levels...
 
 Solving the electrostatics for the left boundary...
 Iter. #  1 took   0.00 msec. with norm =  1.00     at V(  1,  9).
 Iter. #  2 took   0.00 msec. with norm =  1.00     at V(  1,  9).
 Iter. #  3 took   0.00 msec. with norm =  1.00     at V(  1,  9).
 Iter. #  4 took  10.00 msec. with norm =  1.00     at V(  1,  9).
 Iter. #  5 took   0.00 msec. with norm =  1.00     at V(  1,  9).
 Iter. #  6 took   0.00 msec. with norm = 1.000     at V(  1, 10).
 Iter. #  7 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. #  8 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. #  9 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. # 10 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. # 11 took   0.00 msec. with norm = 1.000     at V(  1,  8).
 Iter. # 12 took   0.00 msec. with norm = 1.000     at V(  1,  8).
 Iter. # 13 took   0.00 msec. with norm = 1.000     at V(  1,  8).
 Iter. # 14 took   0.00 msec. with norm = 0.999     at V(  1,  8).
 Iter. # 15 took   0.00 msec. with norm = 0.989     at V(  1,  7).
 Iter. # 16 took   0.00 msec. with norm = 0.964     at V(  1,  7).
 Iter. # 17 took   0.00 msec. with norm = 0.882     at V(  1,  7).
 Iter. # 18 took   0.00 msec. with norm = 0.674     at V(  1,  7).
 Iter. # 19 took   0.00 msec. with norm = 0.333     at V(  1,  7).
 Iter. # 20 took  10.00 msec. with norm = 0.645E-01 at V(  1,  7).
 Iter. # 21 took   0.00 msec. with norm = 0.203E-02 at V(  1,  7).
 Iter. # 22 took   0.00 msec. with norm = 0.193E-05 at V(  1,  7).
 Iter. # 23 took   0.00 msec. with norm = 0.176E-11 at V(  1,  7).
    The Newton solution converged in   23 iterations which took  0.02 sec.
 
 Solving for the exact n and V for the left boundary...
 Iter. #  1 took    0.0 msec. with norm = 0.233E-05 at n(  1, 42).
 Iter. #  2 took    0.0 msec. with norm = 0.280E-07 at n(  1, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 
 Solving the electrostatics for the right boundary...
 Iter. #  1 took   0.00 msec. with norm = 0.665     at V( 32, 21).
 Iter. #  2 took   0.00 msec. with norm = 0.165     at V( 32, 22).
 Iter. #  3 took   0.00 msec. with norm = 0.696E-02 at V( 32, 22).
 Iter. #  4 took  10.00 msec. with norm = 0.119E-04 at V( 32, 22).
 Iter. #  5 took   0.00 msec. with norm = 0.372E-10 at V( 32, 22).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 
 Solving for the exact n and V for the right boundary...
 Iter. #  1 took    0.0 msec. with norm = 0.147E-12 at n( 32, 65).
    The Newton solution converged in   1 iterations which took   0.0 sec.
 
 Solving the electrostatic problem for slice i=  2.
 Iter. #  1 took   0.00 msec. with norm = 0.608E-07 at V(  2, 12).
    The Newton solution converged in    1 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  2.
 Iter. #  1 took    0.0 msec. with norm = 0.233E-05 at n(  2, 42).
 Iter. #  2 took    0.0 msec. with norm = 0.280E-07 at n(  2, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  3.
 Iter. #  1 took   0.00 msec. with norm = 0.368E-01 at V(  3,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.623E-03 at V(  3,  7).
 Iter. #  3 took  10.00 msec. with norm = 0.183E-06 at V(  3,  7).
    The Newton solution converged in    3 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i=  3.
 Iter. #  1 took    0.0 msec. with norm = 0.233E-05 at n(  3, 42).
 Iter. #  2 took    0.0 msec. with norm = 0.280E-07 at n(  3, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  4.
 Iter. #  1 took   0.00 msec. with norm = 0.180     at V(  4,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.135E-01 at V(  4,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.863E-04 at V(  4,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.349E-08 at V(  4,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  4.
 Iter. #  1 took   10.0 msec. with norm = 0.234E-05 at n(  4, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.280E-07 at n(  4, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  5.
 Iter. #  1 took   0.00 msec. with norm = 0.204     at V(  5,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.171E-01 at V(  5,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.139E-03 at V(  5,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.901E-08 at V(  5,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  5.
 Iter. #  1 took    0.0 msec. with norm = 0.236E-05 at n(  5, 41).
 Iter. #  2 took    0.0 msec. with norm = 0.281E-07 at n(  5, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  6.
 Iter. #  1 took   0.00 msec. with norm = 0.231     at V(  6,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.216E-01 at V(  6,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.222E-03 at V(  6,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.231E-07 at V(  6,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  6.
 Iter. #  1 took    0.0 msec. with norm = 0.237E-05 at n(  6, 42).
 Iter. #  2 took    0.0 msec. with norm = 0.282E-07 at n(  6, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  7.
 Iter. #  1 took   0.00 msec. with norm = 0.263     at V(  7,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.274E-01 at V(  7,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.359E-03 at V(  7,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.603E-07 at V(  7,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  7.
 Iter. #  1 took    0.0 msec. with norm = 0.239E-05 at n(  7, 49).
 Iter. #  2 took    0.0 msec. with norm = 0.282E-07 at n(  7, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  8.
 Iter. #  1 took   0.00 msec. with norm = 0.301     at V(  8,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.349E-01 at V(  8,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.585E-03 at V(  8,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.160E-06 at V(  8,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  8.
 Iter. #  1 took    0.0 msec. with norm = 0.241E-05 at n(  8, 51).
 Iter. #  2 took   10.0 msec. with norm = 0.283E-07 at n(  8, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  9.
 Iter. #  1 took   0.00 msec. with norm = 0.345     at V(  9,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.447E-01 at V(  9,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.962E-03 at V(  9,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.433E-06 at V(  9,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  9.
 Iter. #  1 took   10.0 msec. with norm = 0.243E-05 at n(  9, 42).
 Iter. #  2 took    0.0 msec. with norm = 0.284E-07 at n(  9, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 10.
 Iter. #  1 took   0.00 msec. with norm = 0.398     at V( 10,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.574E-01 at V( 10,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.160E-02 at V( 10,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.120E-05 at V( 10,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.670E-12 at V( 10,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 10.
 Iter. #  1 took    0.0 msec. with norm = 0.246E-05 at n( 10, 51).
 Iter. #  2 took   10.0 msec. with norm = 0.285E-07 at n( 10, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 11.
 Iter. #  1 took   0.00 msec. with norm = 0.462     at V( 11,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.742E-01 at V( 11,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.271E-02 at V( 11,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.342E-05 at V( 11,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.547E-11 at V( 11,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 11.
 Iter. #  1 took    0.0 msec. with norm = 0.249E-05 at n( 11, 50).
 Iter. #  2 took    0.0 msec. with norm = 0.286E-07 at n( 11, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 12.
 Iter. #  1 took   0.00 msec. with norm = 0.540     at V( 12,  7).
 Iter. #  2 took  10.00 msec. with norm = 0.964E-01 at V( 12,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.462E-02 at V( 12,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.100E-04 at V( 12,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.466E-10 at V( 12,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 12.
 Iter. #  1 took    0.0 msec. with norm = 0.252E-05 at n( 12, 42).
 Iter. #  2 took    0.0 msec. with norm = 0.287E-07 at n( 12, 32).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 13.
 Iter. #  1 took   0.00 msec. with norm = 0.636     at V( 13,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.126     at V( 13,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.803E-02 at V( 13,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.302E-04 at V( 13,  7).
 Iter. #  5 took  10.00 msec. with norm = 0.423E-09 at V( 13,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 13.
 Iter. #  1 took    0.0 msec. with norm = 0.255E-05 at n( 13, 35).
 Iter. #  2 took    0.0 msec. with norm = 0.289E-07 at n( 13, 35).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 14.
 Iter. #  1 took   0.00 msec. with norm = 0.756     at V( 14,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.165     at V( 14,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.142E-01 at V( 14,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.941E-04 at V( 14,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.411E-08 at V( 14,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 14.
 Iter. #  1 took   10.0 msec. with norm = 0.260E-05 at n( 14, 51).
 Iter. #  2 took    0.0 msec. with norm = 0.290E-07 at n( 14, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 15.
 Iter. #  1 took   0.00 msec. with norm = 0.913     at V( 15,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.219     at V( 15,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.258E-01 at V( 15,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.314E-03 at V( 15,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.457E-07 at V( 15,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 15.
 Iter. #  1 took    0.0 msec. with norm = 0.265E-05 at n( 15, 51).
 Iter. #  2 took   10.0 msec. with norm = 0.292E-07 at n( 15, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 16.
 Iter. #  1 took   0.00 msec. with norm =  1.02     at V( 16,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.273     at V( 16, 14).
 Iter. #  3 took   0.00 msec. with norm = 0.261E-01 at V( 16,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.322E-03 at V( 16,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.480E-07 at V( 16,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 16.
 Iter. #  1 took    0.0 msec. with norm = 0.270E-05 at n( 16, 37).
 Iter. #  2 took    0.0 msec. with norm = 0.294E-07 at n( 16, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 17.
 Iter. #  1 took   0.00 msec. with norm =  1.07     at V( 17,  7).
 Iter. #  2 took  10.00 msec. with norm = 0.312     at V( 17,  4).
 Iter. #  3 took   0.00 msec. with norm = 0.103E-01 at V( 17,  8).
 Iter. #  4 took   0.00 msec. with norm = 0.496E-04 at V( 17,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.113E-08 at V( 17,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 17.
 Iter. #  1 took    0.0 msec. with norm = 0.277E-05 at n( 17, 37).
 Iter. #  2 took    0.0 msec. with norm = 0.296E-07 at n( 17, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 18.
 Iter. #  1 took   0.00 msec. with norm =  1.06     at V( 18,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.310     at V( 18, 18).
 Iter. #  3 took   0.00 msec. with norm = 0.133E-01 at V( 18,  8).
 Iter. #  4 took   0.00 msec. with norm = 0.819E-04 at V( 18,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.307E-08 at V( 18,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 18.
 Iter. #  1 took    0.0 msec. with norm = 0.284E-05 at n( 18, 34).
 Iter. #  2 took    0.0 msec. with norm = 0.298E-07 at n( 18, 34).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 19.
 Iter. #  1 took   0.00 msec. with norm =  1.05     at V( 19,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.298     at V( 19, 17).
 Iter. #  3 took   0.00 msec. with norm = 0.163E-01 at V( 19,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.123E-03 at V( 19,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.687E-08 at V( 19,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 19.
 Iter. #  1 took   10.0 msec. with norm = 0.290E-05 at n( 19, 37).
 Iter. #  2 took    0.0 msec. with norm = 0.300E-07 at n( 19, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 20.
 Iter. #  1 took   0.00 msec. with norm =  1.04     at V( 20,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.289     at V( 20,  5).
 Iter. #  3 took   0.00 msec. with norm = 0.191E-01 at V( 20,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.168E-03 at V( 20,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.127E-07 at V( 20,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 20.
 Iter. #  1 took    0.0 msec. with norm = 0.296E-05 at n( 20, 36).
 Iter. #  2 took   10.0 msec. with norm = 0.301E-07 at n( 20, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 21.
 Iter. #  1 took   0.00 msec. with norm =  1.03     at V( 21,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.276     at V( 21, 15).
 Iter. #  3 took   0.00 msec. with norm = 0.218E-01 at V( 21,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.217E-03 at V( 21,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.211E-07 at V( 21,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 21.
 Iter. #  1 took    0.0 msec. with norm = 0.294E-05 at n( 21, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.256E-07 at n( 21, 22).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 22.
 Iter. #  1 took   0.00 msec. with norm =  1.02     at V( 22,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.264     at V( 22, 14).
 Iter. #  3 took  10.00 msec. with norm = 0.241E-01 at V( 22,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.263E-03 at V( 22,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.307E-07 at V( 22,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 22.
 Iter. #  1 took    0.0 msec. with norm = 0.278E-05 at n( 22, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.182E-07 at n( 22, 21).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 23.
 Iter. #  1 took   0.00 msec. with norm =  1.01     at V( 23,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.251     at V( 23, 12).
 Iter. #  3 took   0.00 msec. with norm = 0.261E-01 at V( 23,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.305E-03 at V( 23,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.409E-07 at V( 23,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 23.
 Iter. #  1 took    0.0 msec. with norm = 0.282E-05 at n( 23, 34).
 Iter. #  2 took    0.0 msec. with norm = 0.182E-07 at n( 23, 21).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 24.
 Iter. #  1 took   0.00 msec. with norm =  1.00     at V( 24,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.240     at V( 24, 10).
 Iter. #  3 took   0.00 msec. with norm = 0.274E-01 at V( 24,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.333E-03 at V( 24,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.480E-07 at V( 24,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 24.
 Iter. #  1 took   10.0 msec. with norm = 0.286E-05 at n( 24, 37).
 Iter. #  2 took    0.0 msec. with norm = 0.182E-07 at n( 24, 21).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 25.
 Iter. #  1 took   0.00 msec. with norm =  1.00     at V( 25,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.234     at V( 25,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.275E-01 at V( 25,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.330E-03 at V( 25,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.467E-07 at V( 25,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 25.
 Iter. #  1 took    0.0 msec. with norm = 0.290E-05 at n( 25, 33).
 Iter. #  2 took   10.0 msec. with norm = 0.148E-01 at n( 25, 21).
 Iter. #  3 took    0.0 msec. with norm = 0.110E-03 at n( 25, 21).
 Iter. #  4 took    0.0 msec. with norm = 0.605E-08 at n( 25, 21).
    The Newton solution converged in   4 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 26.
 Iter. #  1 took   0.00 msec. with norm = 0.985     at V( 26,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.225     at V( 26,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.249E-01 at V( 26,  7).
 Iter. #  4 took  10.00 msec. with norm = 0.266E-03 at V( 26,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.297E-07 at V( 26,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 26.
 Iter. #  1 took    0.0 msec. with norm = 0.243E-05 at n( 26, 34).
 Iter. #  2 took    0.0 msec. with norm = 0.109E-07 at n( 26, 20).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 27.
 Iter. #  1 took   0.00 msec. with norm = 0.980     at V( 27,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.220     at V( 27,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.230E-01 at V( 27,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.221E-03 at V( 27,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.201E-07 at V( 27,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 27.
 Iter. #  1 took   10.0 msec. with norm = 0.246E-05 at n( 27, 36).
 Iter. #  2 took    0.0 msec. with norm = 0.111E-07 at n( 27, 20).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 28.
 Iter. #  1 took   0.00 msec. with norm = 0.993     at V( 28,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.218     at V( 28,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.219E-01 at V( 28,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.193E-03 at V( 28,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.148E-07 at V( 28,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 28.
 Iter. #  1 took   10.0 msec. with norm = 0.249E-05 at n( 28, 37).
 Iter. #  2 took    0.0 msec. with norm = 0.111E-07 at n( 28, 20).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 29.
 Iter. #  1 took   0.00 msec. with norm =  1.01     at V( 29,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.215     at V( 29, 10).
 Iter. #  3 took   0.00 msec. with norm = 0.181E-01 at V( 29,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.126E-03 at V( 29,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.605E-08 at V( 29,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 29.
 Iter. #  1 took    0.0 msec. with norm = 0.251E-05 at n( 29, 37).
 Iter. #  2 took   10.0 msec. with norm = 0.110E-07 at n( 29, 20).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 30.
 Iter. #  1 took   0.00 msec. with norm =  1.02     at V( 30,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.212     at V( 30, 14).
 Iter. #  3 took   0.00 msec. with norm = 0.102E-01 at V( 30,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.370E-04 at V( 30,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.479E-09 at V( 30,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 30.
 Iter. #  1 took    0.0 msec. with norm = 0.254E-05 at n( 30, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.534E-03 at n( 30, 20).
 Iter. #  3 took   10.0 msec. with norm = 0.143E-06 at n( 30, 20).
    The Newton solution converged in   3 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 31.
 Iter. #  1 took   0.00 msec. with norm =  1.04     at V( 31,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.190     at V( 31,  5).
 Iter. #  3 took   0.00 msec. with norm = 0.985E-03 at V( 31,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.283E-06 at V( 31,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 31.
 Iter. #  1 took    0.0 msec. with norm = 0.257E-05 at n( 31, 34).
 Iter. #  2 took    0.0 msec. with norm = 0.225E-02 at n( 31, 20).
 Iter. #  3 took   10.0 msec. with norm = 0.253E-05 at n( 31, 20).
 Iter. #  4 took    0.0 msec. with norm = 0.276E-11 at n( 31, 20).
    The Newton solution converged in   4 iterations which took   0.0 sec.
 
  +------> 1 x nodes/character         maximum voltage =   16.71 kT/q =   0.433 V
  !                    at  32, 18 going v the local dV =    8.63 kT/q =   0.224 V
  !                                        elec. conc. = 0.28E+08 nio = 0.41E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN d
 sNNNNNNNNNNNNNNNNNNNNNn--------- d
 s---------------------           d
 s                                d
 s                                d
 s                                d
 s                                d
 s                                d
 s9999999999999888888888          d
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
 
 
 
 
 Now attempting to use DC() to solve for the biases:
   Vs =   0.000   Vd =   2.000   Vg1 =   4.000   Eavg =     0.0
 Solving for t=  0.000     sec. with dt= 0.1000E+33 sec.
 Iter. #  1 took  0:00+ 0:01 with norm =  33.8     at n( 32,  9).
 Iter. #  2 took  0:00+ 0:01 with norm =  4.32     at V( 32, 21).
 Iter. #  3 took  0:00+ 0:01 with norm =  4.27     at V( 32, 22).
 Iter. #  4 took  0:00+ 0:01 with norm =  4.21     at V( 32, 22).
 Iter. #  5 took  0:00+ 0:01 with norm =  5.66     at n( 31, 22).
 Iter. #  6 took  0:00+ 0:01 with norm =  4.06     at V( 32, 22).
 Iter. #  7 took  0:00+ 0:01 with norm =  3.96     at V( 32, 22).
 Iter. #  8 took  0:00+ 0:01 with norm =  3.82     at V( 32, 22).
 Iter. #  9 took  0:00+ 0:01 with norm =  8.89     at n(  4, 19).
 Iter. # 10 took  0:00+ 0:01 with norm =  23.8     at n( 32, 24).
 Iter. # 11 took  0:00+ 0:01 with norm =  2.62     at V( 32, 23).
 Iter. # 12 took  0:00+ 0:01 with norm =  14.6     at n(  1, 24).
 Iter. # 13 took  0:00+ 0:01 with norm =  13.0     at n(  2, 24).
 Iter. # 14 took  0:00+ 0:01 with norm =  9.14     at n(  3, 25).
 Iter. # 15 took  0:00+ 0:01 with norm =  9.59     at n(  4, 25).
 Iter. # 16 took  0:00+ 0:01 with norm =  7.45     at n(  5, 23).
 Iter. # 17 took  0:00+ 0:01 with norm =  20.7     at n(  2, 20).
 Iter. # 18 took  0:00+ 0:01 with norm =  25.4     at n(  3, 20).
 Iter. # 19 took  0:00+ 0:01 with norm =  4.64     at n(  6, 20).
 Iter. # 20 took  0:00+ 0:01 with norm =  1.22     at n( 32, 25).
 
  +------> 1 x nodes/character         maximum voltage =   42.64 kT/q =   1.104 V
  !                    at  32, 23 going v the local dV =   29.42 kT/q =   0.762 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  d
 sNNNNNNNNNNNNNNNNNNn------------Nd
 sNNN----------------           -Nd
 sN---                          NNd
 sN-                            NNd
 s-                             NNd
 s-                             -Nd
 s                               Nd
 s     88999998888888            Nd
     8ppppppppppppppppppppp8       
  ppppppppppppppppppppppppppp98    
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
 
 
 
 Adding a node to the left of ( 32,  7) because:
   Er=0.22710E+06 V/cm  vr=-.15374E+09 cm/s  v/E= -677.0 cm^2/V-s  hr=0.40E-01 um
   Dr=  18.1 cm^2/sec  2D/v=0.23585E-02 um  h*E=  35.412 kT  dV=  34.248 kT
 A node was added at the right side, so now xnodes= 33 and ynodes= 66.
 Iter. # 21 took  0:00+ 0:01 with norm =  39.0     at n( 33,  7).
 Iter. # 22 took  0:00+ 0:01 with norm =  1.29     at V( 32,  7).
 Iter. # 23 took  0:00+ 0:01 with norm =  1.21     at V( 32,  7).
 Iter. # 24 took  0:00+ 0:01 with norm =  1.13     at V( 32,  7).
 Iter. # 25 took  0:00+ 0:01 with norm =  1.07     at V( 33, 21).
 Iter. # 26 took  0:00+ 0:01 with norm =  1.04     at V( 33, 21).
 Iter. # 27 took  0:00+ 0:01 with norm =  1.02     at V( 33, 21).
 Iter. # 28 took  0:00+ 0:01 with norm =  1.01     at V( 33, 21).
 Iter. # 29 took  0:00+ 0:01 with norm = 0.994     at V( 33, 21).
 Iter. # 30 took  0:00+ 0:01 with norm = 0.955     at V( 33, 21).
 Iter. # 31 took  0:00+ 0:01 with norm = 0.954     at V( 33, 21).
 Iter. # 32 took  0:00+ 0:01 with norm = 0.953     at V( 33, 21).
 Iter. # 33 took  0:00+ 0:01 with norm = 0.952     at V( 33, 21).
 Iter. # 34 took  0:00+ 0:01 with norm = 0.951     at V( 33, 21).
 Iter. # 35 took  0:00+ 0:01 with norm = 0.950     at V( 33, 21).
 Iter. # 36 took  0:00+ 0:01 with norm = 0.948     at V( 33, 21).
 Iter. # 37 took  0:00+ 0:01 with norm = 0.947     at V( 33, 21).
 Iter. # 38 took  0:00+ 0:01 with norm = 0.946     at V( 33, 21).
 Iter. # 39 took  0:00+ 0:01 with norm = 0.944     at V( 33, 21).
 Iter. # 40 took  0:00+ 0:01 with norm = 0.943     at V( 33, 21).
 
  +------> 1 x nodes/character         maximum voltage =   59.93 kT/q =   1.552 V
  !                    at  33, 23 going v the local dV =   39.91 kT/q =   1.034 V
  !                                        elec. conc. = 0.83E+08 nio = 0.12E+19 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  ggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  Nd
 sNNNNNNNNNNNNNNNNNNn-----------  Nd
 sNNN---------------              Nd
 sN---                           -Nd
 sN-                            -NNd
 s-                             -NNd
 s-                              -Nd
 s                                Nd
 s     8899999888888              Nd
     8pppppppppppppppppppp8         
  pppppppppppppppppppppppppp9       
  pppppppppppppppppppppppppppppp999 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 41 took  0:00+ 0:01 with norm = 0.941     at V( 33, 21).
 Iter. # 42 took  0:00+ 0:01 with norm = 0.939     at V( 33, 21).
 Iter. # 43 took  0:00+ 0:01 with norm = 0.937     at V( 33, 21).
 Iter. # 44 took  0:00+ 0:01 with norm = 0.935     at V( 33, 21).
 Iter. # 45 took  0:00+ 0:01 with norm = 0.933     at V( 33, 21).
 Iter. # 46 took  0:00+ 0:01 with norm = 0.931     at V( 33, 21).
 Iter. # 47 took  0:00+ 0:01 with norm = 0.929     at V( 33, 21).
 Iter. # 48 took  0:00+ 0:01 with norm = 0.926     at V( 33, 21).
 Iter. # 49 took  0:00+ 0:01 with norm = 0.923     at V( 33, 21).
 Iter. # 50 took  0:00+ 0:01 with norm = 0.920     at V( 33, 21).
 Iter. # 51 took  0:00+ 0:01 with norm = 0.917     at V( 33, 21).
 Iter. # 52 took  0:00+ 0:01 with norm = 0.914     at V( 33, 21).
 Iter. # 53 took  0:00+ 0:01 with norm = 0.910     at V( 33, 21).
 Iter. # 54 took  0:00+ 0:01 with norm = 0.906     at V( 33, 21).
 Iter. # 55 took  0:00+ 0:01 with norm = 0.901     at V( 33, 21).
 Iter. # 56 took  0:00+ 0:01 with norm = 0.896     at V( 33, 21).
 Iter. # 57 took  0:00+ 0:01 with norm = 0.890     at V( 33, 21).
 Iter. # 58 took  0:00+ 0:01 with norm = 0.884     at V( 33, 21).
 Iter. # 59 took  0:00+ 0:01 with norm = 0.877     at V( 33, 21).
 Iter. # 60 took  0:00+ 0:01 with norm = 0.869     at V( 33, 21).
 
  +------> 1 x nodes/character         maximum voltage =   70.60 kT/q =   1.829 V
  !                    at  33, 23 going v the local dV =   46.44 kT/q =   1.203 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  ggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNn---  -d
 sNNNNNNNNNNNNNNNNNN------------  -d
 sNNN---------------              -d
 sN---                            nd
 sN-                              Nd
 s-                               Nd
 s-                               Nd
 s                                Nd
 s     8899999888888              Nd
     8ppppppppppppppppppp98         
  ppppppppppppppppppppppppp98       
  pppppppppppppppppppppppppppp99888 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 61 took  0:00+ 0:01 with norm = 0.859     at V( 33, 21).
 Iter. # 62 took  0:00+ 0:01 with norm = 0.848     at V( 33, 21).
 Iter. # 63 took  0:00+ 0:01 with norm = 0.835     at V( 33, 21).
 Iter. # 64 took  0:00+ 0:01 with norm = 0.820     at V( 33, 21).
 Iter. # 65 took  0:00+ 0:01 with norm = 0.801     at V( 33, 21).
 Iter. # 66 took  0:00+ 0:01 with norm =  7.95     at n( 32, 14).
 Iter. # 67 took  0:00+ 0:01 with norm =  1.23     at n( 32, 14).
 Iter. # 68 took  0:00+ 0:01 with norm = 0.710     at V( 33, 21).
 Iter. # 69 took  0:00+ 0:01 with norm = 0.657     at V( 33, 21).
 Iter. # 70 took  0:00+ 0:01 with norm = 0.583     at V( 33, 21).
 Iter. # 71 took  0:00+ 0:01 with norm = 0.479     at V( 33, 21).
 Iter. # 72 took  0:00+ 0:01 with norm =  4.67     at n( 31, 17).
 Iter. # 73 took  0:00+ 0:01 with norm =  1.08     at n( 31, 17).
 Iter. # 74 took  0:00+ 0:01 with norm = 0.152     at n( 31, 17).
 Iter. # 75 took  0:00+ 0:01 with norm = 0.210     at n( 31, 17).
 Iter. # 76 took  0:00+ 0:01 with norm =  28.2     at n( 32, 13).
 Iter. # 77 took  0:00+ 0:01 with norm =  34.9     at n( 32,  7).
 Iter. # 78 took  0:00+ 0:01 with norm =  14.4     at n( 29, 19).
 Iter. # 79 took  0:00+ 0:01 with norm =  11.2     at n( 27, 19).
 Iter. # 80 took  0:00+ 0:01 with norm =  4.60     at n( 27, 19).
 
  +------> 1 x nodes/character         maximum voltage =   71.79 kT/q =   1.859 V
  !                    at  33, 24 going ^ the local dV =   46.43 kT/q =   1.203 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  ggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNn------Nd
 sNNNNNNNNNNNNNNNNNn--------    --Nd
 sNNN---------------              Nd
 sN---                            Nd
 sN-                              -d
 s-                               -d
 s-                               -d
 s                                 d
 s     8899999888888               d
     8ppppppppppppppppppp9          
  ppppppppppppppppppppppppp9        
  pppppppppppppppppppppppppppp98888 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
 
 
 
 Adding a node to the right of ( 30, 21) because:
   Er=0.10416E+06 V/cm  vr=-.67930E+08 cm/s  v/E= -652.2 cm^2/V-s  hr=0.45E-01 um
   Dr=  18.1 cm^2/sec  2D/v=0.53378E-02 um  h*E=  18.047 kT  dV=  16.814 kT
 Adding a node to the right of ( 29, 22) because:
   Er= 90890.     V/cm  vr=-.54870E+08 cm/s  v/E= -603.7 cm^2/V-s  hr=0.50E-01 um
   Dr=  18.1 cm^2/sec  2D/v=0.66083E-02 um  h*E=  17.497 kT  dV=  15.090 kT
 Adding a node below ( 31, 23) because:
   Eb=0.11492E+06 V/cm  vb=0.52932E+08 cm/s  v/E=  460.6 cm^2/V-s  hb=0.53E-01 um
   Db=  18.1 cm^2/sec  2D/v=0.68503E-02 um  h*E=  23.475 kT  dV=  15.446 kT
 Adding a node to the right of ( 34, 23) because:
   Er=0.23144E+06 V/cm  vr=-.13589E+09 cm/s  v/E= -587.1 cm^2/V-s  hr=0.20E-01 um
   Dr=  18.1 cm^2/sec  2D/v=0.26684E-02 um  h*E=  18.045 kT  dV=  15.135 kT
 A node was added at the right side, so now xnodes= 36 and ynodes= 67.
 Iter. # 81 took  0:00+ 0:01 with norm =  10.5     at n( 34, 23).
 Iter. # 82 took  0:00+ 0:01 with norm =  1.71     at V( 36, 23).
 Iter. # 83 took  0:00+ 0:01 with norm =  3.13     at n( 15, 21).
 Iter. # 84 took  0:00+ 0:01 with norm =  5.27     at n( 15, 21).
 Iter. # 85 took  0:00+ 0:01 with norm =  1.05     at V( 36, 23).
 Iter. # 86 took  0:00+ 0:01 with norm =  9.18     at n( 34, 23).
 Iter. # 87 took  0:00+ 0:01 with norm = 0.496     at n( 35, 23).
 Iter. # 88 took  0:00+ 0:01 with norm = 0.366E-01 at n( 26, 50).
 Iter. # 89 took  0:00+ 0:01 with norm = 0.583E-01 at n( 17, 56).
 Iter. # 90 took  0:00+ 0:01 with norm = 0.135E-01 at n( 36, 56).
 Iter. # 91 took  0:00+ 0:01 with norm = 0.551E-01 at n( 18, 52).
 Iter. # 92 took  0:00+ 0:01 with norm = 0.301E-02 at V( 18,  7).
 Iter. # 93 took  0:00+ 0:01 with norm = 0.690E-01 at n( 36, 49).
 Iter. # 94 took  0:00+ 0:01 with norm = 0.503E-01 at n( 18, 50).
 Iter. # 95 took  0:00+ 0:01 with norm = 0.128E-02 at V( 16,  7).
 Iter. # 96 took  0:00+ 0:01 with norm = 0.349E-01 at n( 17, 47).
 Iter. # 97 took  0:00+ 0:01 with norm = 0.654E-03 at V( 15,  7).
 Iter. # 98 took  0:00+ 0:01 with norm = 0.770E-01 at n( 23, 50).
 Iter. # 99 took  0:00+ 0:01 with norm = 0.214E-02 at n(  1, 54).
 Iter. #100 took  0:00+ 0:01 with norm = 0.426E-01 at n( 17, 46).
 
  +------> 1 x nodes/character         maximum voltage =   73.20 kT/q =   1.896 V
  !                    at  36, 23 going v the local dV =   35.02 kT/q =   0.907 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.13E-10 nio = 0.19     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNn-------NNNd
 sNNNNNNNNNNNNNNNNNn--------     ---NNd
 sNNN---------------               -NNd
 sN---                             --Nd
 sN-                                -Nd
 s-                                 -Nd
 s-                                  Nd
 s                                   -d
 s     889998888888                   d
      8ppppppppppppppp998              
  889ppppppppppppppppppppp9            
  pppppppppppppppppppppppppp98         
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. #101 took  0:00+ 0:01 with norm = 0.436E-02 at n( 18, 46).
 Iter. #102 took  0:00+ 0:01 with norm = 0.900E-01 at n( 25, 50).
 Iter. #103 took  0:00+ 0:01 with norm = 0.376E-02 at n(  1, 54).
 Iter. #104 took  0:00+ 0:01 with norm = 0.415E-01 at n( 17, 46).
 Iter. #105 took  0:00+ 0:01 with norm = 0.116E-01 at n( 20, 47).
 Iter. #106 took  0:00+ 0:01 with norm = 0.769E-01 at n( 25, 50).
 Iter. #107 took  0:00+ 0:01 with norm = 0.276E-02 at n(  1, 49).
 Iter. #108 took  0:00+ 0:01 with norm = 0.437E-01 at n( 23, 49).
 Iter. #109 took  0:00+ 0:01 with norm = 0.302E-01 at n( 20, 47).
 Iter. #110 took  0:00+ 0:01 with norm = 0.412E-01 at n( 17, 47).
 Iter. #111 took  0:00+ 0:01 with norm = 0.270E-04 at n( 14,  7).
 Iter. #112 took  0:00+ 0:01 with norm = 0.385E-01 at n( 24, 49).
 Iter. #113 took  0:00+ 0:01 with norm = 0.207E-01 at n( 15, 46).
 Iter. #114 took  0:00+ 0:01 with norm = 0.812E-01 at n( 18, 48).
 Iter. #115 took  0:00+ 0:01 with norm = 0.167E-01 at n( 16, 44).
 Iter. #116 took  0:00+ 0:01 with norm = 0.949E-05 at V( 19,  7).
 Iter. #117 took  0:00+ 0:01 with norm = 0.603E-01 at n( 24, 49).
 Iter. #118 took  0:00+ 0:01 with norm = 0.813E-02 at n( 36, 62).
 Iter. #119 took  0:00+ 0:01 with norm = 0.568E-01 at n( 36, 61).
 Iter. #120 took  0:00+ 0:01 with norm = 0.285E-01 at n( 15, 46).
 
  +------> 1 x nodes/character         maximum voltage =   73.20 kT/q =   1.896 V
  !                    at  36, 23 going v the local dV =   35.02 kT/q =   0.907 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.17E-10 nio = 0.26     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNn-------NNNd
 sNNNNNNNNNNNNNNNNNn--------     ---NNd
 sNNN---------------               -NNd
 sN---                             --Nd
 sN-                                -Nd
 s-                                 -Nd
 s-                                  Nd
 s                                   -d
 s     889998888888                   d
      8ppppppppppppppp998              
  889ppppppppppppppppppppp9            
  pppppppppppppppppppppppppp98         
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. #121 took  0:00+ 0:01 with norm = 0.350E-05 at V( 16,  7).
 Iter. #122 took  0:00+ 0:01 with norm = 0.660E-01 at n( 24, 49).
 Iter. #123 took  0:00+ 0:01 with norm = 0.384E-01 at n( 17, 47).
 Iter. #124 took  0:00+ 0:01 with norm = 0.135E-05 at V( 15,  7).
 Iter. #125 took  0:00+ 0:01 with norm = 0.341E-01 at n( 20, 46).
 Iter. #126 took  0:00+ 0:01 with norm = 0.328E-01 at n( 23, 50).
 Iter. #127 took  0:00+ 0:01 with norm = 0.149E-01 at n( 17, 46).
 Iter. #128 took  0:00+ 0:01 with norm = 0.390E-01 at n( 20, 47).
 Iter. #129 took  0:00+ 0:01 with norm = 0.457E-01 at n( 17, 46).
 Iter. #130 took  0:00+ 0:01 with norm = 0.473E-06 at n( 19,  7).
    The Newton solution converged in  130 iterations which took   3:30.
 
 
 ==== Shark   v2.0                                                                   REB 1988         charge map ====
 
  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region
 
 
 
  +------> 1 x nodes/character         maximum voltage =   73.20 kT/q =   1.896 V
  !                    at  36, 23 going v the local dV =   35.02 kT/q =   0.907 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.12E-10 nio = 0.19     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNn-------NNNd
 sNNNNNNNNNNNNNNNNNn--------     ---NNd
 sNNN---------------               -NNd
 sN---                             --Nd
 sN-                                -Nd
 s-                                 -Nd
 s-                                  Nd
 s                                   -d
 s     889998888888                   d
      8ppppppppppppppp998              
  889ppppppppppppppppppppp9            
  pppppppppppppppppppppppppp98         
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
 
 
 
  ...scaling n and v by   0.000
 
 Now attempting to use DC() to solve for the biases:
   Vs =   0.000   Vd =   3.000   Vg1 =   4.000   Eavg =     0.0
 Solving for t=  0.000     sec. with dt= 0.1000E+33 sec.
 Iter. #  1 took  0:00+ 0:01 with norm =  4.21     at V( 36, 24).
 Iter. #  2 took  0:00+ 0:01 with norm =  4.13     at V( 36, 24).
 Iter. #  3 took  0:00+ 0:01 with norm =  4.05     at V( 36, 24).
 Iter. #  4 took  0:00+ 0:01 with norm =  5.49     at n( 34, 23).
 Iter. #  5 took  0:00+ 0:01 with norm =  5.06     at n(  6, 25).
 Iter. #  6 took  0:00+ 0:01 with norm =  15.8     at n(  1, 25).
 Iter. #  7 took  0:00+ 0:01 with norm =  3.64     at V( 36, 24).
 Iter. #  8 took  0:00+ 0:01 with norm =  3.54     at V( 36, 24).
 Iter. #  9 took  0:00+ 0:01 with norm =  3.44     at V( 36, 24).
 Iter. # 10 took  0:00+ 0:01 with norm =  3.33     at V( 36, 24).
 Iter. # 11 took  0:00+ 0:01 with norm =  3.20     at V( 36, 24).
 Iter. # 12 took  0:00+ 0:01 with norm =  2.89     at V( 36, 24).
 Iter. # 13 took  0:00+ 0:01 with norm =  1.78     at V( 36, 24).
 Iter. # 14 took  0:00+ 0:01 with norm =  1.35     at n( 28, 22).
 Iter. # 15 took  0:00+ 0:01 with norm =  1.14     at V( 36, 24).
 Iter. # 16 took  0:00+ 0:01 with norm =  1.09     at V( 36, 24).
 Iter. # 17 took  0:00+ 0:01 with norm =  1.05     at V( 36, 24).
 Iter. # 18 took  0:00+ 0:01 with norm =  1.03     at V( 36, 24).
 Iter. # 19 took  0:00+ 0:01 with norm =  1.00     at V( 36, 24).
 Iter. # 20 took  0:00+ 0:01 with norm = 0.933     at V( 36, 24).
 
  +------> 1 x nodes/character         maximum voltage =   94.74 kT/q =   2.454 V
  !                    at  36, 24 going v the local dV =   40.30 kT/q =   1.044 V
  !                                        elec. conc. = 0.20E+09 nio = 0.30E+19 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNn--------nNd
 sNNNNNNNNNNNNNNNNNn----------------NNd
 sNNN---------------              --NNd
 sN---                             -NNd
 sN-                                -Nd
 s-                                 -Nd
 s-                                 -Nd
 s                                  -Nd
 s      89998888888                  Nd
      89ppppppppppppp998               
    89ppppppppppppppppppp9             
  ppppppppppppppppppppppppp9           
  pppppppppppppppppppppppppppp99888888 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPPpppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPpppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPpppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPpPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 21 took  0:00+ 0:01 with norm = 0.915     at V( 36, 24).
 Iter. # 22 took  0:00+ 0:01 with norm = 0.910     at V( 36, 24).
 Iter. # 23 took  0:00+ 0:01 with norm = 0.906     at V( 36, 24).
 Iter. # 24 took  0:00+ 0:01 with norm = 0.901     at V( 36, 24).
 Iter. # 25 took  0:00+ 0:01 with norm = 0.895     at V( 36, 24).
 Iter. # 26 took  0:00+ 0:01 with norm = 0.888     at V( 36, 24).
 Iter. # 27 took  0:00+ 0:01 with norm = 0.880     at V( 36, 24).
 Iter. # 28 took  0:00+ 0:01 with norm = 0.872     at V( 36, 24).
 Iter. # 29 took  0:00+ 0:01 with norm = 0.861     at V( 36, 24).
 Iter. # 30 took  0:00+ 0:01 with norm = 0.849     at V( 36, 24).
 Iter. # 31 took  0:00+ 0:01 with norm = 0.833     at V( 36, 24).
 Iter. # 32 took  0:00+ 0:01 with norm = 0.814     at V( 36, 24).
 Iter. # 33 took  0:00+ 0:01 with norm = 0.790     at V( 36, 24).
 Iter. # 34 took  0:00+ 0:01 with norm = 0.757     at V( 36, 24).
 Iter. # 35 took  0:00+ 0:01 with norm = 0.720     at V( 36, 16).
 Iter. # 36 took  0:00+ 0:01 with norm = 0.682     at V( 36, 16).
 Iter. # 37 took  0:00+ 0:01 with norm = 0.631     at V( 36, 16).
 Iter. # 38 took  0:00+ 0:01 with norm = 0.563     at V( 36, 16).
 Iter. # 39 took  0:00+ 0:01 with norm = 0.469     at V( 36, 16).
 Iter. # 40 took  0:00+ 0:01 with norm = 0.345     at V( 36, 16).
 
  +------> 1 x nodes/character         maximum voltage =  109.85 kT/q =   2.845 V
  !                    at  36, 24 going v the local dV =   45.05 kT/q =   1.167 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNn--------- -d
 sNNNNNNNNNNNNNNNNNn-------         -Nd
 sNNN--------------                  Nd
 sN---                               Nd
 sN-                                 Nd
 s-                                  Nd
 s-                                  Nd
 s                                   Nd
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp9             
  ppppppppppppppppppppppppp8           
  ppppppppppppppppppppppppppp9888      
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPPpppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 41 took  0:00+ 0:01 with norm = 0.204     at V( 36, 16).
 Iter. # 42 took  0:00+ 0:01 with norm = 0.850E-01 at V( 36, 16).
 Iter. # 43 took  0:00+ 0:01 with norm = 0.227E-01 at V( 36, 16).
 Iter. # 44 took  0:00+ 0:01 with norm = 0.449E-01 at n( 17, 47).
 Iter. # 45 took  0:00+ 0:01 with norm = 0.426E-01 at n( 36, 58).
 Iter. # 46 took  0:00+ 0:01 with norm = 0.313E-01 at n( 23, 51).
 Iter. # 47 took  0:00+ 0:01 with norm = 0.176E-01 at n( 15, 46).
 Iter. # 48 took  0:00+ 0:01 with norm = 0.280E-01 at n( 22, 55).
 Iter. # 49 took  0:00+ 0:01 with norm = 0.540E-01 at n( 32, 54).
 Iter. # 50 took  0:00+ 0:01 with norm = 0.123E-02 at V( 16,  7).
 Iter. # 51 took  0:00+ 0:01 with norm = 0.184E-01 at n( 23, 50).
 Iter. # 52 took  0:00+ 0:01 with norm = 0.499E-01 at n( 18, 58).
 Iter. # 53 took  0:00+ 0:01 with norm = 0.478E-01 at n( 36, 47).
 Iter. # 54 took  0:00+ 0:01 with norm = 0.305E-01 at n( 16, 45).
 Iter. # 55 took  0:00+ 0:01 with norm = 0.174E-01 at n( 24, 51).
 Iter. # 56 took  0:00+ 0:01 with norm = 0.223E-01 at n( 15, 46).
 Iter. # 57 took  0:00+ 0:01 with norm = 0.880E-01 at n( 18, 47).
 Iter. # 58 took  0:00+ 0:01 with norm = 0.202E-03 at n( 19,  7).
 Iter. # 59 took  0:00+ 0:01 with norm = 0.793E-02 at n(  5, 64).
 Iter. # 60 took  0:00+ 0:01 with norm = 0.616E-01 at n( 19, 48).
 
  +------> 1 x nodes/character         maximum voltage =  110.16 kT/q =   2.853 V
  !                    at  36, 24 going v the local dV =   45.06 kT/q =   1.167 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.86E-11 nio = 0.13     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNn--------  -d
 sNNNNNNNNNNNNNNNNNn-------          -d
 sNNN--------------                  Nd
 sN---                               Nd
 sN-                                 Nd
 s-                                  Nd
 s-                                  Nd
 s                                   -d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp9             
  ppppppppppppppppppppppppp8           
  ppppppppppppppppppppppppppp9888      
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPPpppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 61 took  0:00+ 0:01 with norm = 0.658E-01 at n( 36, 52).
 Iter. # 62 took  0:00+ 0:01 with norm = 0.970E-01 at n( 36, 60).
 Iter. # 63 took  0:00+ 0:01 with norm = 0.667E-04 at n( 16,  7).
 Iter. # 64 took  0:00+ 0:01 with norm = 0.361E-01 at n( 36, 56).
 Iter. # 65 took  0:00+ 0:01 with norm = 0.697E-01 at n( 36, 48).
 Iter. # 66 took  0:00+ 0:01 with norm = 0.192E-01 at n( 17, 45).
 Iter. # 67 took  0:00+ 0:01 with norm = 0.133E-01 at n( 24, 51).
 Iter. # 68 took  0:00+ 0:01 with norm = 0.182E-01 at n( 15, 45).
 Iter. # 69 took  0:00+ 0:01 with norm = 0.809E-01 at n( 18, 47).
 Iter. # 70 took  0:00+ 0:01 with norm = 0.921E-05 at V( 19,  7).
 Iter. # 71 took  0:00+ 0:01 with norm = 0.367E-01 at n( 18, 46).
 Iter. # 72 took  0:00+ 0:01 with norm = 0.377E-01 at n( 35, 51).
 Iter. # 73 took  0:00+ 0:01 with norm = 0.618E-05 at V( 18,  7).
 Iter. # 74 took  0:00+ 0:01 with norm = 0.511E-01 at n( 18, 47).
 Iter. # 75 took  0:00+ 0:01 with norm = 0.281E-01 at n( 15, 46).
 Iter. # 76 took  0:00+ 0:01 with norm = 0.314E-05 at V( 16,  7).
 Iter. # 77 took  0:00+ 0:01 with norm = 0.667E-01 at n( 18, 48).
 Iter. # 78 took  0:00+ 0:01 with norm = 0.173E-05 at V( 15,  7).
 Iter. # 79 took  0:00+ 0:01 with norm = 0.357E-01 at n( 17, 47).
 Iter. # 80 took  0:00+ 0:01 with norm = 0.424E-01 at n( 35, 54).
 
  +------> 1 x nodes/character         maximum voltage =  110.16 kT/q =   2.853 V
  !                    at  36, 24 going v the local dV =   45.06 kT/q =   1.167 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.16E-10 nio = 0.24     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNn--------  -d
 sNNNNNNNNNNNNNNNNNn-------          -d
 sNNN--------------                  Nd
 sN---                               Nd
 sN-                                 Nd
 s-                                  Nd
 s-                                  Nd
 s                                   -d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp9             
  ppppppppppppppppppppppppp8           
  ppppppppppppppppppppppppppp9888      
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPPpppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPpppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 81 took  0:00+ 0:01 with norm = 0.358E-01 at n( 15, 46).
 Iter. # 82 took  0:00+ 0:01 with norm = 0.576E-06 at n( 21,  7).
    The Newton solution converged in   82 iterations which took   2:34.
 
 
 ==== Shark   v2.0                                                                   REB 1988         charge map ====
 
  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region
 
 
 
  +------> 1 x nodes/character         maximum voltage =  110.16 kT/q =   2.853 V
  !                    at  36, 24 going v the local dV =   45.06 kT/q =   1.167 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.15E-10 nio = 0.22     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNn--------  -d
 sNNNNNNNNNNNNNNNNNn-------          -d
 sNNN--------------                  Nd
 sN---                               Nd
 sN-                                 Nd
 s-                                  Nd
 s-                                  Nd
 s                                   -d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp9             
  ppppppppppppppppppppppppp8           
  ppppppppppppppppppppppppppp9888      
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPPpppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPpppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
  ...scaling n and v by   0.500
 
 Now attempting to use DC() to solve for the biases:
   Vs =   0.000   Vd =   3.500   Vg1 =   4.000   Eavg =     0.0
 Solving for t=  0.000     sec. with dt= 0.1000E+33 sec.
 Iter. #  1 took  0:00+ 0:01 with norm =  21.6     at n( 35, 24).
 Iter. #  2 took  0:00+ 0:01 with norm =  33.7     at n( 34, 11).
 Iter. #  3 took  0:00+ 0:01 with norm =  35.5     at n( 33,  7).
 Iter. #  4 took  0:00+ 0:01 with norm =  37.8     at n( 36, 16).
 Iter. #  5 took  0:00+ 0:01 with norm =  1.91     at V( 36, 24).
 Iter. #  6 took  0:00+ 0:01 with norm =  1.58     at V( 36, 24).
 Iter. #  7 took  0:00+ 0:01 with norm =  1.27     at V( 36, 24).
 Iter. #  8 took  0:00+ 0:01 with norm =  1.15     at V( 36, 24).
 Iter. #  9 took  0:00+ 0:01 with norm =  1.09     at n( 34, 24).
 Iter. # 10 took  0:00+ 0:01 with norm =  10.0     at n( 36, 24).
 Iter. # 11 took  0:00+ 0:01 with norm =  14.0     at n( 36, 24).
 Iter. # 12 took  0:00+ 0:01 with norm = 0.152E-02 at n( 16,  7).
 Iter. # 13 took  0:00+ 0:01 with norm = 0.490E-01 at n( 23, 51).
 Iter. # 14 took  0:00+ 0:01 with norm = 0.867E-01 at n( 36, 54).
 Iter. # 15 took  0:00+ 0:01 with norm = 0.518E-03 at n( 14,  7).
 Iter. # 16 took  0:00+ 0:01 with norm = 0.268E-01 at n( 36, 61).
 Iter. # 17 took  0:00+ 0:01 with norm = 0.481E-01 at n( 36, 51).
 Iter. # 18 took  0:00+ 0:01 with norm = 0.253E-03 at V( 20,  7).
 Iter. # 19 took  0:00+ 0:01 with norm = 0.604E-01 at n( 18, 47).
 Iter. # 20 took  0:00+ 0:01 with norm = 0.482E-01 at n( 36, 53).
 
  +------> 1 x nodes/character         maximum voltage =  128.67 kT/q =   3.332 V
  !                    at  36, 24 going v the local dV =   49.62 kT/q =   1.285 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.48E-10 nio = 0.71     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN--------    d
 sNNNNNNNNNNNNNNNNN-------           -d
 sNNN--------------                  -d
 sN---                               -d
 sN-                                 -d
 s-                                  -d
 s-                                  -d
 s                                   -d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp8             
  pppppppppppppppppppppppp9            
  pppppppppppppppppppppppppp98         
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPpppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 21 took  0:00+ 0:01 with norm = 0.215E-01 at n( 23, 50).
 Iter. # 22 took  0:00+ 0:01 with norm = 0.399E-01 at n( 15, 46).
 Iter. # 23 took  0:00+ 0:01 with norm = 0.888E-02 at n( 23, 50).
 Iter. # 24 took  0:00+ 0:01 with norm = 0.416E-01 at n( 17, 47).
 Iter. # 25 took  0:00+ 0:01 with norm = 0.307E-01 at n( 24, 51).
 Iter. # 26 took  0:00+ 0:01 with norm = 0.259E-01 at n( 17, 46).
 Iter. # 27 took  0:00+ 0:01 with norm = 0.293E-01 at n( 23, 50).
 Iter. # 28 took  0:00+ 0:01 with norm = 0.247E-01 at n( 17, 46).
 Iter. # 29 took  0:00+ 0:01 with norm = 0.376E-01 at n( 36, 58).
 Iter. # 30 took  0:00+ 0:01 with norm = 0.845E-02 at n( 17, 46).
 Iter. # 31 took  0:00+ 0:01 with norm = 0.570E-01 at n( 31, 53).
 Iter. # 32 took  0:00+ 0:01 with norm = 0.128E-04 at n( 20,  7).
 Iter. # 33 took  0:00+ 0:01 with norm = 0.489E-02 at n(  1, 55).
 Iter. # 34 took  0:00+ 0:01 with norm = 0.439E-01 at n( 16, 46).
 Iter. # 35 took  0:00+ 0:01 with norm = 0.103E-01 at n( 16, 45).
 Iter. # 36 took  0:00+ 0:01 with norm = 0.806E-01 at n( 18, 47).
 Iter. # 37 took  0:00+ 0:01 with norm = 0.521E-05 at n( 17,  7).
 Iter. # 38 took  0:00+ 0:01 with norm = 0.156E-01 at n( 17, 46).
 Iter. # 39 took  0:00+ 0:01 with norm = 0.494E-01 at n( 16, 46).
 Iter. # 40 took  0:00+ 0:01 with norm = 0.345E-01 at n( 21, 52).
 
  +------> 1 x nodes/character         maximum voltage =  128.67 kT/q =   3.332 V
  !                    at  36, 24 going v the local dV =   49.62 kT/q =   1.285 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.47E-10 nio = 0.71     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN--------    d
 sNNNNNNNNNNNNNNNNN-------           -d
 sNNN--------------                  -d
 sN---                               -d
 sN-                                 -d
 s-                                  -d
 s-                                  -d
 s                                   -d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp8             
  pppppppppppppppppppppppp9            
  pppppppppppppppppppppppppp98         
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPpppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 41 took  0:00+ 0:01 with norm = 0.191E-01 at n( 32, 52).
 Iter. # 42 took  0:00+ 0:01 with norm = 0.448E-01 at n( 18, 56).
 Iter. # 43 took  0:00+ 0:01 with norm = 0.675E-06 at n( 13,  7).
    The Newton solution converged in   43 iterations which took   1:21.
 
 
 ==== Shark   v2.0                                                                   REB 1988         charge map ====
 
  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region
 
 
 
  +------> 1 x nodes/character         maximum voltage =  128.67 kT/q =   3.332 V
  !                    at  36, 24 going v the local dV =   49.62 kT/q =   1.285 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.46E-10 nio = 0.69     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN--------    d
 sNNNNNNNNNNNNNNNNN-------           -d
 sNNN--------------                  -d
 sN---                               -d
 sN-                                 -d
 s-                                  -d
 s-                                  -d
 s                                   -d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp8             
  pppppppppppppppppppppppp9            
  pppppppppppppppppppppppppp98         
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPpppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
  ...scaling n and v by   1.000
 
 Now attempting to use DC() to solve for the biases:
   Vs =   0.000   Vd =   4.000   Vg1 =   4.000   Eavg =     0.0
 Solving for t=  0.000     sec. with dt= 0.1000E+33 sec.
 Iter. #  1 took  0:00+ 0:01 with norm =  7.39     at n( 35, 24).
 Iter. #  2 took  0:00+ 0:01 with norm =  29.2     at n( 36, 23).
 Iter. #  3 took  0:00+ 0:01 with norm =  3.28     at V( 36, 24).
 Iter. #  4 took  0:00+ 0:01 with norm =  13.4     at n(  1, 26).
 Iter. #  5 took  0:00+ 0:01 with norm =  3.00     at V( 36, 24).
 Iter. #  6 took  0:00+ 0:01 with norm =  2.83     at V( 36, 24).
 Iter. #  7 took  0:00+ 0:01 with norm =  2.65     at V( 36, 24).
 Iter. #  8 took  0:00+ 0:01 with norm =  2.44     at V( 36, 24).
 Iter. #  9 took  0:00+ 0:01 with norm =  2.21     at V( 36, 24).
 Iter. # 10 took  0:00+ 0:01 with norm =  1.96     at V( 36, 24).
 Iter. # 11 took  0:00+ 0:01 with norm =  1.66     at V( 36, 24).
 Iter. # 12 took  0:00+ 0:01 with norm =  1.31     at V( 36, 24).
 Iter. # 13 took  0:00+ 0:01 with norm =  1.20     at V( 36, 24).
 Iter. # 14 took  0:00+ 0:01 with norm =  1.06     at V( 36, 24).
 Iter. # 15 took  0:00+ 0:01 with norm = 0.716     at V( 36, 12).
 Iter. # 16 took  0:00+ 0:01 with norm = 0.673     at V( 36, 12).
 Iter. # 17 took  0:00+ 0:01 with norm = 0.616     at V( 36, 13).
 Iter. # 18 took  0:00+ 0:01 with norm = 0.539     at V( 36, 14).
 Iter. # 19 took  0:00+ 0:01 with norm = 0.434     at V( 36, 14).
 Iter. # 20 took  0:00+ 0:01 with norm = 0.299     at V( 36, 15).
 
  +------> 1 x nodes/character         maximum voltage =  147.02 kT/q =   3.808 V
  !                    at  36, 24 going v the local dV =   54.10 kT/q =   1.401 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN-------    -d
 sNNNNNNNNNNNNNNNNN-------           Nd
 sNNN--------------                  Nd
 sN---                               Nd
 sN-                                 Nd
 s-                                  Nd
 s-                                  Nd
 s                                   nd
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp              
  pppppppppppppppppppppppp8            
  pppppppppppppppppppppppppp8          
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPpppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 21 took  0:00+ 0:01 with norm = 0.155     at V( 36, 15).
 Iter. # 22 took  0:00+ 0:01 with norm = 0.497E-01 at V( 36, 15).
 Iter. # 23 took  0:00+ 0:01 with norm = 0.243E-01 at n( 15, 46).
 Iter. # 24 took  0:00+ 0:01 with norm =  1.43     at n( 34, 23).
 Iter. # 25 took  0:00+ 0:01 with norm =  1.15     at n( 34, 23).
 Iter. # 26 took  0:00+ 0:01 with norm = 0.153     at n( 34, 23).
 Iter. # 27 took  0:00+ 0:01 with norm = 0.105     at n( 36, 53).
 Iter. # 28 took  0:00+ 0:01 with norm = 0.207E-01 at n( 34, 23).
 Iter. # 29 took  0:00+ 0:01 with norm = 0.242E-01 at n( 23, 50).
 Iter. # 30 took  0:00+ 0:01 with norm = 0.702E-01 at n( 36, 62).
 Iter. # 31 took  0:00+ 0:01 with norm = 0.149E-01 at n( 16, 45).
 Iter. # 32 took  0:00+ 0:01 with norm = 0.353E-01 at n( 18, 56).
 Iter. # 33 took  0:00+ 0:01 with norm = 0.241E-03 at V( 13,  7).
 Iter. # 34 took  0:00+ 0:01 with norm = 0.736E-01 at n( 35, 65).
 Iter. # 35 took  0:00+ 0:01 with norm = 0.580E-01 at n( 17, 47).
 Iter. # 36 took  0:00+ 0:01 with norm = 0.163E-03 at n( 20,  7).
 Iter. # 37 took  0:00+ 0:01 with norm = 0.227E-01 at n( 15, 45).
 Iter. # 38 took  0:00+ 0:01 with norm = 0.827E-01 at n( 36, 54).
 Iter. # 39 took  0:00+ 0:01 with norm = 0.103E-03 at n( 18,  7).
 Iter. # 40 took  0:00+ 0:01 with norm = 0.435E-01 at n( 18, 48).
 
  +------> 1 x nodes/character         maximum voltage =  147.22 kT/q =   3.813 V
  !                    at  36, 24 going v the local dV =   54.10 kT/q =   1.401 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN-------     d
 sNNNNNNNNNNNNNNNNN-------            d
 sNNN--------------                  -d
 sN---                               -d
 sN-                                 -d
 s-                                  -d
 s-                                  -d
 s                                    d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp              
  pppppppppppppppppppppppp8            
  pppppppppppppppppppppppppp8          
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPpppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 41 took  0:00+ 0:01 with norm = 0.334E-01 at n( 22, 50).
 Iter. # 42 took  0:00+ 0:01 with norm =  2.41     at n( 36, 24).
 Iter. # 43 took  0:00+ 0:01 with norm =  19.7     at n( 36, 24).
 Iter. # 44 took  0:00+ 0:01 with norm = 0.224     at n( 35, 24).
 Iter. # 45 took  0:00+ 0:01 with norm = 0.209E-04 at n( 15,  7).
 Iter. # 46 took  0:00+ 0:01 with norm = 0.210E-01 at n( 24, 52).
 Iter. # 47 took  0:00+ 0:01 with norm = 0.827E-01 at n( 18, 47).
 Iter. # 48 took  0:00+ 0:01 with norm = 0.625E-05 at V( 20,  7).
 Iter. # 49 took  0:00+ 0:01 with norm = 0.250E-01 at n( 24, 51).
 Iter. # 50 took  0:00+ 0:01 with norm = 0.525E-01 at n( 17, 57).
 Iter. # 51 took  0:00+ 0:01 with norm = 0.384E-01 at n( 36, 61).
 Iter. # 52 took  0:00+ 0:01 with norm = 0.287E-01 at n( 16, 46).
 Iter. # 53 took  0:00+ 0:01 with norm = 0.485E-02 at n( 23, 51).
 Iter. # 54 took  0:00+ 0:01 with norm = 0.190E-01 at n( 15, 45).
 Iter. # 55 took  0:00+ 0:01 with norm = 0.800E-01 at n( 18, 47).
 Iter. # 56 took  0:00+ 0:01 with norm = 0.187E-05 at V( 16,  7).
 Iter. # 57 took  0:00+ 0:01 with norm = 0.169E-01 at n( 15, 45).
 Iter. # 58 took  0:00+ 0:01 with norm = 0.606E-01 at n( 18, 47).
 Iter. # 59 took  0:00+ 0:01 with norm = 0.648E-02 at n( 23, 50).
 Iter. # 60 took  0:00+ 0:01 with norm = 0.699E-01 at n( 18, 47).
 
  +------> 1 x nodes/character         maximum voltage =  147.22 kT/q =   3.813 V
  !                    at  36, 24 going v the local dV =   54.10 kT/q =   1.401 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.22E-10 nio = 0.34     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN-------     d
 sNNNNNNNNNNNNNNNNN-------            d
 sNNN--------------                  -d
 sN---                               -d
 sN-                                 -d
 s-                                  -d
 s-                                  -d
 s                                    d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp              
  pppppppppppppppppppppppp8            
  pppppppppppppppppppppppppp8          
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPpppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 -- Note: --
    While the device warranted a mesh refinement, the solution was
    close enough that it is believed that this mesh can be used
    to obtain the solution.
 Iter. # 61 took  0:00+ 0:01 with norm = 0.368E-06 at n( 21,  7).
    The Newton solution converged in   61 iterations which took   1:55.
 
 
 ==== Shark   v2.0                                                                   REB 1988         charge map ====
 
  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region
 
 
 
  +------> 1 x nodes/character         maximum voltage =  147.22 kT/q =   3.813 V
  !                    at  36, 24 going v the local dV =   54.10 kT/q =   1.401 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.19E-10 nio = 0.28     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN-------     d
 sNNNNNNNNNNNNNNNNN-------            d
 sNNN--------------                  -d
 sN---                               -d
 sN-                                 -d
 s-                                  -d
 s-                                  -d
 s                                    d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp              
  pppppppppppppppppppppppp8            
  pppppppppppppppppppppppppp8          
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPpppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPpp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 
 
 ==== Shark   v2.0                                                                  REB 1988 starting transient ====
 
 Outputting the simulation results for t =   0.00 (psec).
 Iter. #  1 took  0:00+ 0:01 with norm = 0.136E-01 at n( 15, 46).
 Iter. #  2 took  0:00+ 0:01 with norm = 0.780E-01 at n( 36, 57).
 Iter. #  3 took  0:00+ 0:01 with norm = 0.757E-02 at n( 36, 60).
 Iter. #  4 took  0:00+ 0:01 with norm = 0.270E-01 at n( 15, 46).
 Iter. #  5 took  0:00+ 0:01 with norm = 0.496E-01 at n( 18, 47).
 Iter. #  6 took  0:00+ 0:01 with norm = 0.168E-01 at n( 15, 46).
 Iter. #  7 took  0:00+ 0:01 with norm = 0.179E-01 at n( 23, 50).
 Iter. #  8 took  0:00+ 0:01 with norm = 0.478E-01 at n( 35, 60).
 Iter. #  9 took  0:00+ 0:01 with norm = 0.785E-07 at n( 16,  7).
    The Newton solution converged in    9 iterations which took   0:16.
 
 
 ==== Shark   v2.0                                                                   REB 1988         charge map ====
 
  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region
 
 
 
  +------> 1 x nodes/character         maximum voltage =  147.22 kT/q =   3.813 V
  !                    at  36, 24 going v the local dV =   54.10 kT/q =   1.401 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.19E-10 nio = 0.29     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNN-------     d
 sNNNNNNNNNNNNNNNNN-------            d
 sNNN--------------                  -d
 sN---                               -d
 sN-                                 -d
 s-                                  -d
 s-                                  -d
 s                                    d
 s      89998888888                   d
      89ppppppppppppp998               
    89ppppppppppppppppppp              
  pppppppppppppppppppppppp8            
  pppppppppppppppppppppppppp8          
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppppp 
  pppppppppPPPpppppppppppppppppppppppp 
  ppppppPPPPPPPPPPPPPppppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPppppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPpppppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPppppppppppp 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
  PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP 
 
 
 
 Outputting the simulation results for t = ****** (psec).
  Total CPU time used:  0 days, 00:09:40
