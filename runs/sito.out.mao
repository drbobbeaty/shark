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
 * MOSFET with SiO2, biased into strong inversion, triode mode
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
 bias vsource=0.00 vdrain=1.00 eavg=0.0 vg1=4.0
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
         drain voltage =  1.000 (V).
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
 Iter. #  4 took   0.00 msec. with norm =  1.00     at V(  1,  9).
 Iter. #  5 took   0.00 msec. with norm =  1.00     at V(  1,  9).
 Iter. #  6 took   0.00 msec. with norm = 1.000     at V(  1, 10).
 Iter. #  7 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. #  8 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. #  9 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. # 10 took   0.00 msec. with norm = 1.000     at V(  1,  9).
 Iter. # 11 took   0.00 msec. with norm = 1.000     at V(  1,  8).
 Iter. # 12 took   0.00 msec. with norm = 1.000     at V(  1,  8).
 Iter. # 13 took   0.00 msec. with norm = 1.000     at V(  1,  8).
 Iter. # 14 took  10.00 msec. with norm = 0.999     at V(  1,  8).
 Iter. # 15 took   0.00 msec. with norm = 0.989     at V(  1,  7).
 Iter. # 16 took   0.00 msec. with norm = 0.964     at V(  1,  7).
 Iter. # 17 took   0.00 msec. with norm = 0.882     at V(  1,  7).
 Iter. # 18 took   0.00 msec. with norm = 0.674     at V(  1,  7).
 Iter. # 19 took   0.00 msec. with norm = 0.333     at V(  1,  7).
 Iter. # 20 took   0.00 msec. with norm = 0.645E-01 at V(  1,  7).
 Iter. # 21 took   0.00 msec. with norm = 0.203E-02 at V(  1,  7).
 Iter. # 22 took   0.00 msec. with norm = 0.193E-05 at V(  1,  7).
 Iter. # 23 took   0.00 msec. with norm = 0.176E-11 at V(  1,  7).
    The Newton solution converged in   23 iterations which took  0.01 sec.
 
 Solving for the exact n and V for the left boundary...
 Iter. #  1 took    0.0 msec. with norm = 0.233E-05 at n(  1, 42).
 Iter. #  2 took   10.0 msec. with norm = 0.280E-07 at n(  1, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 
 Solving the electrostatics for the right boundary...
 Iter. #  1 took   0.00 msec. with norm = 0.665     at V( 32, 21).
 Iter. #  2 took   0.00 msec. with norm = 0.165     at V( 32, 22).
 Iter. #  3 took   0.00 msec. with norm = 0.696E-02 at V( 32, 22).
 Iter. #  4 took   0.00 msec. with norm = 0.119E-04 at V( 32, 22).
 Iter. #  5 took   0.00 msec. with norm = 0.372E-10 at V( 32, 22).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 
 Solving for the exact n and V for the right boundary...
 Iter. #  1 took    0.0 msec. with norm = 0.135E-05 at n( 32, 18).
 Iter. #  2 took    0.0 msec. with norm = 0.894E-02 at n( 32, 18).
 Iter. #  3 took   10.0 msec. with norm = 0.398E-04 at n( 32, 18).
 Iter. #  4 took    0.0 msec. with norm = 0.735E-09 at n( 32, 18).
    The Newton solution converged in   4 iterations which took   0.0 sec.
 
 Solving the electrostatic problem for slice i=  2.
 Iter. #  1 took   0.00 msec. with norm = 0.977E-02 at V(  2,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.446E-04 at V(  2,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.936E-09 at V(  2,  7).
    The Newton solution converged in    3 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  2.
 Iter. #  1 took    0.0 msec. with norm = 0.233E-05 at n(  2, 41).
 Iter. #  2 took   10.0 msec. with norm = 0.280E-07 at n(  2, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  3.
 Iter. #  1 took   0.00 msec. with norm = 0.158     at V(  3,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.106E-01 at V(  3,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.536E-04 at V(  3,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.135E-08 at V(  3,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  3.
 Iter. #  1 took    0.0 msec. with norm = 0.234E-05 at n(  3, 42).
 Iter. #  2 took    0.0 msec. with norm = 0.280E-07 at n(  3, 40).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  4.
 Iter. #  1 took   0.00 msec. with norm = 0.179     at V(  4,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.134E-01 at V(  4,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.850E-04 at V(  4,  7).
 Iter. #  4 took  10.00 msec. with norm = 0.339E-08 at V(  4,  7).
    The Newton solution converged in    4 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i=  4.
 Iter. #  1 took    0.0 msec. with norm = 0.235E-05 at n(  4, 41).
 Iter. #  2 took    0.0 msec. with norm = 0.281E-07 at n(  4, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  5.
 Iter. #  1 took   0.00 msec. with norm = 0.203     at V(  5,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.169E-01 at V(  5,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.136E-03 at V(  5,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.867E-08 at V(  5,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  5.
 Iter. #  1 took   10.0 msec. with norm = 0.237E-05 at n(  5, 40).
 Iter. #  2 took    0.0 msec. with norm = 0.281E-07 at n(  5, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  6.
 Iter. #  1 took   0.00 msec. with norm = 0.230     at V(  6,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.214E-01 at V(  6,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.219E-03 at V(  6,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.224E-07 at V(  6,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  6.
 Iter. #  1 took    0.0 msec. with norm = 0.238E-05 at n(  6, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.282E-07 at n(  6, 42).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  7.
 Iter. #  1 took  10.00 msec. with norm = 0.262     at V(  7,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.272E-01 at V(  7,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.354E-03 at V(  7,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.586E-07 at V(  7,  7).
    The Newton solution converged in    4 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i=  7.
 Iter. #  1 took    0.0 msec. with norm = 0.240E-05 at n(  7, 34).
 Iter. #  2 took    0.0 msec. with norm = 0.283E-07 at n(  7, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  8.
 Iter. #  1 took   0.00 msec. with norm = 0.300     at V(  8,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.347E-01 at V(  8,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.577E-03 at V(  8,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.156E-06 at V(  8,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  8.
 Iter. #  1 took   10.0 msec. with norm = 0.242E-05 at n(  8, 44).
 Iter. #  2 took    0.0 msec. with norm = 0.283E-07 at n(  8, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  9.
 Iter. #  1 took   0.00 msec. with norm = 0.344     at V(  9,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.444E-01 at V(  9,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.952E-03 at V(  9,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.424E-06 at V(  9,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  9.
 Iter. #  1 took    0.0 msec. with norm = 0.244E-05 at n(  9, 36).
 Iter. #  2 took   10.0 msec. with norm = 0.284E-07 at n(  9, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 10.
 Iter. #  1 took   0.00 msec. with norm = 0.397     at V( 10,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.570E-01 at V( 10,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.158E-02 at V( 10,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.117E-05 at V( 10,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.636E-12 at V( 10,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 10.
 Iter. #  1 took    0.0 msec. with norm = 0.247E-05 at n( 10, 38).
 Iter. #  2 took    0.0 msec. with norm = 0.285E-07 at n( 10, 42).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 11.
 Iter. #  1 took  10.00 msec. with norm = 0.460     at V( 11,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.737E-01 at V( 11,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.266E-02 at V( 11,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.331E-05 at V( 11,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.512E-11 at V( 11,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 11.
 Iter. #  1 took    0.0 msec. with norm = 0.249E-05 at n( 11, 63).
 Iter. #  2 took    0.0 msec. with norm = 0.287E-07 at n( 11, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 12.
 Iter. #  1 took   0.00 msec. with norm = 0.537     at V( 12,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.956E-01 at V( 12,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.454E-02 at V( 12,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.965E-05 at V( 12,  7).
 Iter. #  5 took  10.00 msec. with norm = 0.434E-10 at V( 12,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 12.
 Iter. #  1 took    0.0 msec. with norm = 0.253E-05 at n( 12, 55).
 Iter. #  2 took    0.0 msec. with norm = 0.288E-07 at n( 12, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 13.
 Iter. #  1 took   0.00 msec. with norm = 0.633     at V( 13,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.125     at V( 13,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.791E-02 at V( 13,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.293E-04 at V( 13,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.398E-09 at V( 13,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 13.
 Iter. #  1 took   10.0 msec. with norm = 0.256E-05 at n( 13, 65).
 Iter. #  2 took    0.0 msec. with norm = 0.289E-07 at n( 13, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 14.
 Iter. #  1 took   0.00 msec. with norm = 0.753     at V( 14,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.164     at V( 14,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.140E-01 at V( 14,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.918E-04 at V( 14,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.391E-08 at V( 14,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 14.
 Iter. #  1 took    0.0 msec. with norm = 0.261E-05 at n( 14, 51).
 Iter. #  2 took   10.0 msec. with norm = 0.291E-07 at n( 14, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 15.
 Iter. #  1 took   0.00 msec. with norm = 0.908     at V( 15,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.217     at V( 15,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.254E-01 at V( 15,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.303E-03 at V( 15,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.426E-07 at V( 15,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 15.
 Iter. #  1 took    0.0 msec. with norm = 0.265E-05 at n( 15, 64).
 Iter. #  2 took    0.0 msec. with norm = 0.293E-07 at n( 15, 31).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 16.
 Iter. #  1 took   0.00 msec. with norm =  1.02     at V( 16,  7).
 Iter. #  2 took  10.00 msec. with norm = 0.274     at V( 16, 14).
 Iter. #  3 took   0.00 msec. with norm = 0.264E-01 at V( 16,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.329E-03 at V( 16,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.499E-07 at V( 16,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 16.
 Iter. #  1 took    0.0 msec. with norm = 0.271E-05 at n( 16, 36).
 Iter. #  2 took    0.0 msec. with norm = 0.295E-07 at n( 16, 36).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 17.
 Iter. #  1 took   0.00 msec. with norm =  1.07     at V( 17,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.309     at V( 17,  4).
 Iter. #  3 took   0.00 msec. with norm = 0.106E-01 at V( 17,  8).
 Iter. #  4 took   0.00 msec. with norm = 0.524E-04 at V( 17,  7).
 Iter. #  5 took  10.00 msec. with norm = 0.126E-08 at V( 17,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 17.
 Iter. #  1 took    0.0 msec. with norm = 0.278E-05 at n( 17, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.297E-07 at n( 17, 33).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 18.
 Iter. #  1 took   0.00 msec. with norm =  1.06     at V( 18,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.311     at V( 18, 18).
 Iter. #  3 took   0.00 msec. with norm = 0.136E-01 at V( 18,  8).
 Iter. #  4 took   0.00 msec. with norm = 0.857E-04 at V( 18,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.336E-08 at V( 18,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 18.
 Iter. #  1 took   10.0 msec. with norm = 0.284E-05 at n( 18, 34).
 Iter. #  2 took    0.0 msec. with norm = 0.298E-07 at n( 18, 34).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 19.
 Iter. #  1 took   0.00 msec. with norm =  1.05     at V( 19,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.299     at V( 19, 17).
 Iter. #  3 took   0.00 msec. with norm = 0.166E-01 at V( 19,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.127E-03 at V( 19,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.730E-08 at V( 19,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 19.
 Iter. #  1 took    0.0 msec. with norm = 0.291E-05 at n( 19, 36).
 Iter. #  2 took   10.0 msec. with norm = 0.300E-07 at n( 19, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 20.
 Iter. #  1 took   0.00 msec. with norm =  1.04     at V( 20,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.287     at V( 20, 16).
 Iter. #  3 took   0.00 msec. with norm = 0.194E-01 at V( 20,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.173E-03 at V( 20,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.136E-07 at V( 20,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 20.
 Iter. #  1 took    0.0 msec. with norm = 0.297E-05 at n( 20, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.301E-07 at n( 20, 30).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 21.
 Iter. #  1 took   0.00 msec. with norm =  1.03     at V( 21,  7).
 Iter. #  2 took  10.00 msec. with norm = 0.275     at V( 21, 15).
 Iter. #  3 took   0.00 msec. with norm = 0.221E-01 at V( 21,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.223E-03 at V( 21,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.222E-07 at V( 21,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 21.
 Iter. #  1 took    0.0 msec. with norm = 0.294E-05 at n( 21, 35).
 Iter. #  2 took    0.0 msec. with norm = 0.256E-07 at n( 21, 22).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 22.
 Iter. #  1 took   0.00 msec. with norm =  1.02     at V( 22,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.261     at V( 22, 14).
 Iter. #  3 took   0.00 msec. with norm = 0.244E-01 at V( 22,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.270E-03 at V( 22,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.324E-07 at V( 22,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 22.
 Iter. #  1 took    0.0 msec. with norm = 0.278E-05 at n( 22, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.182E-07 at n( 22, 21).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 23.
 Iter. #  1 took   0.00 msec. with norm =  1.01     at V( 23,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.249     at V( 23, 11).
 Iter. #  3 took   0.00 msec. with norm = 0.264E-01 at V( 23,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.311E-03 at V( 23,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.426E-07 at V( 23,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 23.
 Iter. #  1 took   10.0 msec. with norm = 0.282E-05 at n( 23, 37).
 Iter. #  2 took    0.0 msec. with norm = 0.182E-07 at n( 23, 21).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 24.
 Iter. #  1 took   0.00 msec. with norm =  1.00     at V( 24,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.240     at V( 24,  9).
 Iter. #  3 took   0.00 msec. with norm = 0.277E-01 at V( 24,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.340E-03 at V( 24,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.500E-07 at V( 24,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 24.
 Iter. #  1 took    0.0 msec. with norm = 0.286E-05 at n( 24, 33).
 Iter. #  2 took   10.0 msec. with norm = 0.182E-07 at n( 24, 21).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 25.
 Iter. #  1 took   0.00 msec. with norm = 0.994     at V( 25,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.232     at V( 25,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.269E-01 at V( 25,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.315E-03 at V( 25,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.424E-07 at V( 25,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 25.
 Iter. #  1 took    0.0 msec. with norm = 0.290E-05 at n( 25, 36).
 Iter. #  2 took    0.0 msec. with norm = 0.174E-01 at n( 25, 21).
 Iter. #  3 took   10.0 msec. with norm = 0.150E-03 at n( 25, 21).
 Iter. #  4 took    0.0 msec. with norm = 0.114E-07 at n( 25, 21).
    The Newton solution converged in   4 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 26.
 Iter. #  1 took   0.00 msec. with norm = 0.978     at V( 26,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.223     at V( 26,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.243E-01 at V( 26,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.251E-03 at V( 26,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.265E-07 at V( 26,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 26.
 Iter. #  1 took    0.0 msec. with norm = 0.244E-05 at n( 26, 42).
 Iter. #  2 took   10.0 msec. with norm = 0.109E-07 at n( 26, 20).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 27.
 Iter. #  1 took   0.00 msec. with norm = 0.973     at V( 27,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.217     at V( 27,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.224E-01 at V( 27,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.208E-03 at V( 27,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.177E-07 at V( 27,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 27.
 Iter. #  1 took    0.0 msec. with norm = 0.246E-05 at n( 27, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.111E-07 at n( 27, 20).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 28.
 Iter. #  1 took   0.00 msec. with norm = 0.984     at V( 28,  7).
 Iter. #  2 took  10.00 msec. with norm = 0.214     at V( 28,  7).
 Iter. #  3 took   0.00 msec. with norm = 0.211E-01 at V( 28,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.179E-03 at V( 28,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.127E-07 at V( 28,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 28.
 Iter. #  1 took    0.0 msec. with norm = 0.249E-05 at n( 28, 33).
 Iter. #  2 took    0.0 msec. with norm = 0.111E-07 at n( 28, 20).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 29.
 Iter. #  1 took   0.00 msec. with norm =  1.00     at V( 29,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.213     at V( 29,  9).
 Iter. #  3 took   0.00 msec. with norm = 0.186E-01 at V( 29,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.132E-03 at V( 29,  7).
 Iter. #  5 took   0.00 msec. with norm = 0.654E-08 at V( 29,  7).
    The Newton solution converged in    5 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 29.
 Iter. #  1 took    0.0 msec. with norm = 0.252E-05 at n( 29, 34).
 Iter. #  2 took    0.0 msec. with norm = 0.532E-05 at n( 29, 20).
 Iter. #  3 took    0.0 msec. with norm = 0.182E-09 at n( 29, 20).
    The Newton solution converged in   3 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 30.
 Iter. #  1 took   0.00 msec. with norm =  1.02     at V( 30,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.208     at V( 30, 13).
 Iter. #  3 took   0.00 msec. with norm = 0.106E-01 at V( 30,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.396E-04 at V( 30,  7).
 Iter. #  5 took  10.00 msec. with norm = 0.544E-09 at V( 30,  7).
    The Newton solution converged in    5 iterations which took  0.01 sec.
 Solving the exact n and V problem for slice i= 30.
 Iter. #  1 took    0.0 msec. with norm = 0.254E-05 at n( 30, 37).
 Iter. #  2 took    0.0 msec. with norm = 0.587E-03 at n( 30, 20).
 Iter. #  3 took    0.0 msec. with norm = 0.172E-06 at n( 30, 20).
    The Newton solution converged in   3 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 31.
 Iter. #  1 took   0.00 msec. with norm =  1.04     at V( 31,  7).
 Iter. #  2 took   0.00 msec. with norm = 0.184     at V( 31,  5).
 Iter. #  3 took   0.00 msec. with norm = 0.123E-02 at V( 31,  7).
 Iter. #  4 took   0.00 msec. with norm = 0.433E-06 at V( 31,  7).
    The Newton solution converged in    4 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 31.
 Iter. #  1 took   10.0 msec. with norm = 0.257E-05 at n( 31, 32).
 Iter. #  2 took    0.0 msec. with norm = 0.240E-02 at n( 31, 20).
 Iter. #  3 took    0.0 msec. with norm = 0.287E-05 at n( 31, 20).
 Iter. #  4 took    0.0 msec. with norm = 0.374E-11 at n( 31, 20).
    The Newton solution converged in   4 iterations which took   0.0 sec.
 
  +------> 1 x nodes/character         maximum voltage =   16.71 kT/q =   0.433 V
  !                    at  32, 19 going ^ the local dV =    8.63 kT/q =   0.224 V
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
   Vs =   0.000   Vd =   1.000   Vg1 =   4.000   Eavg =     0.0
 Solving for t=  0.000     sec. with dt= 0.1000E+33 sec.
 Iter. #  1 took  0:00+ 0:01 with norm =  16.7     at n( 32, 18).
 Iter. #  2 took  0:00+ 0:01 with norm =  3.80     at V( 32, 21).
 Iter. #  3 took  0:00+ 0:01 with norm =  3.75     at V( 32, 22).
 Iter. #  4 took  0:00+ 0:01 with norm =  3.67     at V( 32, 22).
 Iter. #  5 took  0:00+ 0:01 with norm =  3.57     at V( 32, 22).
 Iter. #  6 took  0:00+ 0:01 with norm =  5.50     at n( 31, 22).
 Iter. #  7 took  0:00+ 0:01 with norm =  3.34     at V( 32, 22).
 Iter. #  8 took  0:00+ 0:01 with norm =  3.21     at V( 32, 22).
 Iter. #  9 took  0:00+ 0:01 with norm =  9.82     at n(  4, 19).
 Iter. # 10 took  0:00+ 0:01 with norm =  7.70     at n(  3, 19).
 Iter. # 11 took  0:00+ 0:01 with norm =  2.63     at V( 32, 22).
 Iter. # 12 took  0:00+ 0:01 with norm =  14.6     at n(  1, 24).
 Iter. # 13 took  0:00+ 0:01 with norm =  13.0     at n(  2, 24).
 Iter. # 14 took  0:00+ 0:01 with norm =  9.14     at n(  3, 25).
 Iter. # 15 took  0:00+ 0:01 with norm =  9.59     at n(  4, 25).
 Iter. # 16 took  0:00+ 0:01 with norm =  8.20     at n(  8, 22).
 Iter. # 17 took  0:00+ 0:01 with norm =  21.1     at n(  2, 20).
 Iter. # 18 took  0:00+ 0:01 with norm =  25.4     at n(  3, 20).
 Iter. # 19 took  0:00+ 0:01 with norm =  5.27     at n(  8, 21).
 Iter. # 20 took  0:00+ 0:01 with norm = 0.851     at V( 32, 20).
 
  +------> 1 x nodes/character         maximum voltage =   35.10 kT/q =   0.909 V
  !                    at  32, 23 going v the local dV =   25.65 kT/q =   0.664 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNd
 sNNNNNNNNNNNNNNNNNNn-----------NNd
 sNNN----------------           -Nd
 sN---                          -Nd
 sN-                            -Nd
 s-                             -Nd
 s-                              Nd
 s                               Nd
 s     88999998888888             d
     8ppppppppppppppppppppp9       
  ppppppppppppppppppppppppppp98    
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppp 
 
 
 
 Adding a node to the left of ( 32, 19) because:
   Er= 98972.     V/cm  vr=-.67735E+08 cm/s  v/E= -684.4 cm^2/V-s  hr=0.40E-01 um
   Dr=  18.1 cm^2/sec  2D/v=0.53532E-02 um  h*E=15.433 kT  dV=15.089 kT
 A node was added at the right side, so now xnodes= 33 and ynodes= 66.
 Iter. # 21 took  0:00+ 0:01 with norm =  2.44     at V( 33, 23).
 Iter. # 22 took  0:00+ 0:01 with norm =  2.17     at V( 33, 23).
 Iter. # 23 took  0:00+ 0:01 with norm =  1.87     at V( 33, 23).
 Iter. # 24 took  0:00+ 0:01 with norm =  1.51     at V( 33, 23).
 Iter. # 25 took  0:00+ 0:01 with norm =  1.25     at V( 33, 23).
 Iter. # 26 took  0:00+ 0:01 with norm =  1.14     at V( 33, 23).
 Iter. # 27 took  0:00+ 0:01 with norm =  1.02     at V( 33, 23).
 Iter. # 28 took  0:00+ 0:01 with norm = 0.795     at V( 33, 21).
 Iter. # 29 took  0:00+ 0:01 with norm = 0.772     at V( 33, 21).
 Iter. # 30 took  0:00+ 0:01 with norm = 0.742     at V( 33, 21).
 Iter. # 31 took  0:00+ 0:01 with norm = 0.705     at V( 33, 21).
 Iter. # 32 took  0:00+ 0:01 with norm = 0.654     at V( 33, 21).
 Iter. # 33 took  0:00+ 0:01 with norm = 0.584     at V( 33, 21).
 Iter. # 34 took  0:00+ 0:01 with norm = 0.486     at V( 33, 21).
 Iter. # 35 took  0:00+ 0:01 with norm = 0.355     at V( 33, 21).
 Iter. # 36 took  0:00+ 0:01 with norm = 0.202     at V( 33, 21).
 Iter. # 37 took  0:00+ 0:01 with norm = 0.759E-01 at V( 33, 21).
 Iter. # 38 took  0:00+ 0:01 with norm = 0.696E-01 at n( 26, 48).
 Iter. # 39 took  0:00+ 0:01 with norm =  12.4     at n( 30, 20).
 Iter. # 40 took  0:00+ 0:01 with norm =  16.3     at n( 30, 19).
 
  +------> 1 x nodes/character         maximum voltage =   36.55 kT/q =   0.947 V
  !                    at  33, 24 going ^ the local dV =   32.74 kT/q =   0.848 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V 
 
  ggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNd
 sNNNNNNNNNNNNNNNNNNn-----------NNNd
 sNNN---------------           -nNNd
 sN---                          -NNd
 sN-                             -Nd
 s-                              -Nd
 s-                               nd
 s                                -d
 s     88999998888888              d
     8pppppppppppppppppppp98        
  ppppppppppppppppppppppppppp8      
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
  ppppppppppppppppppppppppppppppppp 
 
 
 
 Adding a node below ( 31, 23) because:
   Eb=0.12579E+06 V/cm  vb=0.57849E+08 cm/s  v/E=  459.9 cm^2/V-s  hb=0.53E-01 um
   Db=  18.1 cm^2/sec  2D/v=0.62680E-02 um  h*E=25.694 kT  dV=16.881 kT
 Adding a node to the left of ( 33, 23) because:
   Er=0.22824E+06 V/cm  vr=-.13752E+09 cm/s  v/E= -602.5 cm^2/V-s  hr=0.20E-01 um
   Dr=  18.1 cm^2/sec  2D/v=0.26368E-02 um  h*E=17.795 kT  dV=15.317 kT
 A node was added at the right side, so now xnodes= 34 and ynodes= 67.
 Iter. # 41 took  0:00+ 0:01 with norm =  32.1     at n( 32, 21).
 Iter. # 42 took  0:00+ 0:01 with norm =  1.56     at n( 26, 20).
 Iter. # 43 took  0:00+ 0:01 with norm =  2.64     at n( 15, 21).
 Iter. # 44 took  0:00+ 0:01 with norm =  6.05     at n( 15, 21).
 Iter. # 45 took  0:00+ 0:01 with norm =  2.60     at n( 32, 24).
 Iter. # 46 took  0:00+ 0:01 with norm =  4.30     at n( 33, 24).
 Iter. # 47 took  0:00+ 0:01 with norm = 0.182     at n( 33, 24).
 Iter. # 48 took  0:00+ 0:01 with norm = 0.134E-01 at n( 33, 24).
 Iter. # 49 took  0:00+ 0:01 with norm =  1.90     at n( 34, 24).
 Iter. # 50 took  0:00+ 0:01 with norm = 0.416E-02 at V( 20,  7).
 Iter. # 51 took  0:00+ 0:01 with norm = 0.434E-01 at n( 34, 53).
 Iter. # 52 took  0:00+ 0:01 with norm = 0.661E-01 at n( 18, 47).
 Iter. # 53 took  0:00+ 0:01 with norm = 0.204E-01 at n( 16, 45).
 Iter. # 54 took  0:00+ 0:01 with norm = 0.314E-01 at n( 16, 47).
 Iter. # 55 took  0:00+ 0:01 with norm = 0.101E-02 at V( 16,  7).
 Iter. # 56 took  0:00+ 0:01 with norm = 0.612E-01 at n( 20, 48).
 Iter. # 57 took  0:00+ 0:01 with norm = 0.543E-01 at n( 17, 46).
 Iter. # 58 took  0:00+ 0:01 with norm = 0.335E-03 at V( 14,  7).
 Iter. # 59 took  0:00+ 0:01 with norm = 0.556E-01 at n( 32, 51).
 Iter. # 60 took  0:00+ 0:01 with norm = 0.427E-01 at n( 16, 48).
 
  +------> 1 x nodes/character         maximum voltage =   36.79 kT/q =   0.953 V
  !                    at  34, 23 going v the local dV =   24.97 kT/q =   0.647 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.14E-10 nio = 0.22     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNd
 sNNNNNNNNNNNNNNNNNNn----------NNNNNd
 sNNN---------------           -NNNNd
 sN---                          -NNNd
 sN-                            --NNd
 s-                              -NNd
 s-                               -Nd
 s                                 Nd
 s     8899988888888                d
      8ppppppppppppppp9998           
  889pppppppppppppppppppppp98        
  pppppppppppppppppppppppppppp988888 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
 
 
 
 Adding a node below ( 33, 23) because:
   Eb=0.25032E+06 V/cm  vb=0.12475E+09 cm/s  v/E=  498.4 cm^2/V-s  hb=0.26E-01 um
   Db=  18.1 cm^2/sec  2D/v=0.29065E-02 um  h*E=25.566 kT  dV=18.202 kT
 Nodes were added, so now xnodes= 34 and ynodes= 68.
 Iter. # 61 took  0:00+ 0:01 with norm =  1.61     at V( 34, 24).
 Iter. # 62 took  0:00+ 0:01 with norm =  2.46     at n( 34, 26).
 Iter. # 63 took  0:00+ 0:01 with norm =  1.16     at V( 34, 24).
 Iter. # 64 took  0:00+ 0:01 with norm =  1.87     at n( 33, 24).
 Iter. # 65 took  0:00+ 0:01 with norm =  14.0     at n( 34, 24).
 Iter. # 66 took  0:00+ 0:01 with norm =  3.54     at n( 34, 24).
 Iter. # 67 took  0:00+ 0:01 with norm = 0.183     at n( 34, 25).
 Iter. # 68 took  0:00+ 0:01 with norm = 0.134E-02 at n( 22,  7).
 Iter. # 69 took  0:00+ 0:01 with norm = 0.877E-01 at n( 24, 50).
 Iter. # 70 took  0:00+ 0:01 with norm = 0.333E-01 at n( 17, 48).
 Iter. # 71 took  0:00+ 0:01 with norm = 0.132E-01 at n( 33, 57).
 Iter. # 72 took  0:00+ 0:01 with norm = 0.270E-01 at n( 15, 46).
 Iter. # 73 took  0:00+ 0:01 with norm = 0.216E-01 at n( 14, 46).
 Iter. # 74 took  0:00+ 0:01 with norm = 0.396E-01 at n( 17, 47).
 Iter. # 75 took  0:00+ 0:01 with norm = 0.187E-03 at n( 16,  7).
 Iter. # 76 took  0:00+ 0:01 with norm = 0.331E-01 at n( 20, 47).
 Iter. # 77 took  0:00+ 0:01 with norm = 0.521E-01 at n( 18, 51).
 Iter. # 78 took  0:00+ 0:01 with norm = 0.612E-04 at n( 14,  7).
 Iter. # 79 took  0:00+ 0:01 with norm = 0.607E-01 at n( 25, 49).
 Iter. # 80 took  0:00+ 0:01 with norm = 0.407E-04 at V( 22,  7).
 
  +------> 1 x nodes/character         maximum voltage =   36.79 kT/q =   0.953 V
  !                    at  34, 23 going v the local dV =   16.00 kT/q =   0.415 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.25E-10 nio = 0.38     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNd
 sNNNNNNNNNNNNNNNNNNn----------NNNNNd
 sNNN---------------           -NNNNd
 sN---                          -NNNd
 sN-                            --NNd
 s-                              -NNd
 s-                               -Nd
 s                                 Nd
 s     8899988888888                d
       9ppppppp9999998888            
     8pppppppppppppppppppp98         
  ppppppppppppppppppppppppppp8       
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
 
 
 
 Adding a node below ( 34, 23) because:
   Eb=0.40006E+06 V/cm  vb=0.21937E+09 cm/s  v/E=  548.4 cm^2/V-s  hb=0.13E-01 um
   Db=  18.1 cm^2/sec  2D/v=0.16529E-02 um  h*E=20.430 kT  dV=16.004 kT
 Nodes were added, so now xnodes= 34 and ynodes= 69.
 Iter. # 81 took  0:00+ 0:01 with norm =  1.22     at V( 34, 24).
 Iter. # 82 took  0:00+ 0:01 with norm =  1.09     at n( 34, 23).
 Iter. # 83 took  0:00+ 0:01 with norm = 0.395     at V( 34, 24).
 Iter. # 84 took  0:00+ 0:01 with norm = 0.566     at n( 34, 24).
 Iter. # 85 took  0:00+ 0:01 with norm = 0.515E-01 at n( 34, 50).
 Iter. # 86 took  0:00+ 0:01 with norm = 0.683E-01 at n( 20, 49).
 Iter. # 87 took  0:00+ 0:01 with norm = 0.973E-03 at n( 23,  7).
 Iter. # 88 took  0:00+ 0:01 with norm = 0.453E-01 at n( 26, 51).
 Iter. # 89 took  0:00+ 0:01 with norm = 0.221E-01 at n( 15, 47).
 Iter. # 90 took  0:00+ 0:01 with norm = 0.177E-01 at n( 34, 67).
 Iter. # 91 took  0:00+ 0:01 with norm = 0.321E-01 at n( 15, 47).
 Iter. # 92 took  0:00+ 0:01 with norm = 0.249E-03 at n( 19,  7).
 Iter. # 93 took  0:00+ 0:01 with norm = 0.435E-01 at n( 16, 47).
 Iter. # 94 took  0:00+ 0:01 with norm = 0.354E-01 at n( 16, 48).
 Iter. # 95 took  0:00+ 0:01 with norm = 0.103E-03 at n( 17,  7).
 Iter. # 96 took  0:00+ 0:01 with norm = 0.535E-01 at n( 26, 51).
 Iter. # 97 took  0:00+ 0:01 with norm = 0.217E-01 at n( 34, 58).
 Iter. # 98 took  0:00+ 0:01 with norm = 0.728E-01 at n( 17, 48).
 Iter. # 99 took  0:00+ 0:01 with norm = 0.225E-04 at V( 22,  7).
 Iter. #100 took  0:00+ 0:01 with norm = 0.551E-01 at n( 26, 51).
 
  +------> 1 x nodes/character         maximum voltage =   36.79 kT/q =   0.953 V
  !                    at  31, 22 going < the local dV =   13.77 kT/q =   0.357 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.37E-10 nio = 0.56     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNd
 sNNNNNNNNNNNNNNNNNNn----------NNNNNd
 sNNN---------------           -NNNNd
 sN---                          -NNNd
 sN-                            --NNd
 s-                              -NNd
 s-                               -Nd
 s                                 Nd
 s     8899988888888                d
       8999999999988888              
      8pppppppppppppppp998           
  899ppppppppppppppppppppppp8        
  pppppppppppppppppppppppppppp998888 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
 
 
 
 Iter. #101 took  0:00+ 0:01 with norm = 0.203E-04 at V( 21,  7).
 Iter. #102 took  0:00+ 0:01 with norm = 0.328E-01 at n( 26, 51).
 Iter. #103 took  0:00+ 0:01 with norm = 0.471E-01 at n( 17, 49).
 Iter. #104 took  0:00+ 0:01 with norm = 0.924E-02 at n( 15, 47).
 Iter. #105 took  0:00+ 0:01 with norm = 0.204E-01 at n( 19, 49).
 Iter. #106 took  0:00+ 0:01 with norm = 0.217E-01 at n( 15, 47).
 Iter. #107 took  0:00+ 0:01 with norm = 0.498E-01 at n( 34, 57).
 Iter. #108 took  0:00+ 0:01 with norm = 0.496E-02 at n(  3, 48).
 Iter. #109 took  0:00+ 0:01 with norm = 0.455E-01 at n( 17, 48).
 Iter. #110 took  0:00+ 0:01 with norm = 0.322E-01 at n( 16, 47).
 Iter. #111 took  0:00+ 0:01 with norm = 0.147E-01 at n( 34, 64).
 Iter. #112 took  0:00+ 0:01 with norm = 0.468E-01 at n( 17, 48).
 Iter. #113 took  0:00+ 0:01 with norm = 0.114E-01 at n(  1, 51).
 Iter. #114 took  0:00+ 0:01 with norm = 0.222E-01 at n( 17, 49).
 Iter. #115 took  0:00+ 0:01 with norm = 0.226E-01 at n( 15, 47).
 Iter. #116 took  0:00+ 0:01 with norm = 0.201E-01 at n( 25, 50).
 Iter. #117 took  0:00+ 0:01 with norm = 0.579E-01 at n( 17, 48).
 Iter. #118 took  0:00+ 0:01 with norm = 0.641E-06 at n( 19,  7).
    The Newton solution converged in  118 iterations which took   2:58.
 
 
 ==== Shark   v2.0                                                                   REB 1988         charge map ====
 
  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region
 
 
 
  +------> 1 x nodes/character         maximum voltage =   36.79 kT/q =   0.953 V
  !                    at  31, 22 going < the local dV =   13.77 kT/q =   0.357 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.30E-10 nio = 0.44     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNd
 sNNNNNNNNNNNNNNNNNNn----------NNNNNd
 sNNN---------------           -NNNNd
 sN---                          -NNNd
 sN-                            --NNd
 s-                              -NNd
 s-                               -Nd
 s                                 Nd
 s     8899988888888                d
       8999999999988888              
      8pppppppppppppppp998           
  899ppppppppppppppppppppppp8        
  pppppppppppppppppppppppppppp998888 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
 
 
 
 
 
 ==== Shark   v2.0                                                                  REB 1988 starting transient ====
 
 Outputting the simulation results for t =   0.00 (psec).
 Iter. #  1 took  0:00+ 0:01 with norm = 0.378E-02 at n( 14, 46).
 Iter. #  2 took  0:00+ 0:01 with norm = 0.174E-01 at n( 14, 47).
 Iter. #  3 took  0:00+ 0:01 with norm = 0.905E-01 at n( 34, 55).
 Iter. #  4 took  0:00+ 0:01 with norm = 0.111E-01 at n( 34, 62).
 Iter. #  5 took  0:00+ 0:01 with norm = 0.542E-01 at n( 34, 56).
 Iter. #  6 took  0:00+ 0:01 with norm = 0.111E-06 at n( 15,  7).
    The Newton solution converged in    6 iterations which took   0:09.
 
 
 ==== Shark   v2.0                                                                   REB 1988         charge map ====
 
  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region
 
 
 
  +------> 1 x nodes/character         maximum voltage =   36.79 kT/q =   0.953 V
  !                    at  30, 22 going > the local dV =   13.77 kT/q =   0.357 V
  !                                        elec. conc. = 0.52E+08 nio = 0.78E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.28E-10 nio = 0.43     /cm**3
  V 
 
  gggggggggggggggggggggggggggggggggg 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii 
 sNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNd
 sNNNNNNNNNNNNNNNNNNn----------NNNNNd
 sNNN---------------           -NNNNd
 sN---                          -NNNd
 sN-                            --NNd
 s-                              -NNd
 s-                               -Nd
 s                                 Nd
 s     8899988888888                d
       8999999999988888              
      8pppppppppppppppp998           
  899ppppppppppppppppppppppp8        
  pppppppppppppppppppppppppppp998888 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
  pppppppppppppppppppppppppppppppppp 
 
 
 
 Outputting the simulation results for t = ****** (psec).
  Total CPU time used:  0 days, 00:03:10
