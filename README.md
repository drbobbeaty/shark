# Shark - A Thesis Project

<p align="center">
  <img src="doc/img/shark.png" width="256" height="256" border="0" />
</p>

In 1988, I received my Ph.D. from Purdue University, and as part of that
thesis, this code was written to simulate, in two dimensions, and time, the
internal state of the device as it was under bias. The details of the code
are really best described in the Thesis, which is still available from
_University Microfilms_, which is interesting to me today... but hey, it's not
a lot of storage, so why not?

## Input File Format

_Shark_ is based on an _input deck_ like so many codebases that were meant
to run on supercomputers without even so much as terminal access. One such
input file is:
```
*
* This will be a comparative test of Sutherland's work for a Silicon
* MOSFET with SiO2, biased into weak inversion, triode mode
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
bias vsource=0.00 vdrain=0.35 eavg=0.0 vg1=2.5
* do the steady state operation: delt -> infinity
time tout=10.0 tstop=100.0 delt=1.0e38
solution type=linpk
start with=DC maxit=300 debug=sortof guess=MOS
```
and simulates a simple MISFET transistor in weak-inversion, triode operation.
This was a similar test case given in Sutherland's work, and was used to check
that _Shark_ was matching his work for correctness.

When run with something like:
```bash
$ shark < wito > wito.out
```
you will get something like:
```
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
 * MOSFET with SiO2, biased into weak inversion, triode mode
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
 bias vsource=0.00 vdrain=0.35 eavg=0.0 vg1=2.5
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
         drain voltage =  0.350 (V).
         source end of the gate =  2.500 (V).
         drain end of the gate =  2.500 (V).
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
       Vgmax =  2.500 V
  Vg+Qss/Cox =  2.500 V
         Cox = 0.345E-07 F/cm**2
          Vo = 2.7840 V
          V1 = 0.0000 V
        phis =  0.629 V
          xd =  0.202 um

 Because Shark is selecting the mesh, Shark is setting the thickness of the
 substrate to  4.413 um

    automatically selecting the y mesh for a MODFET channel.

    placing y nodes uniformly in tox...
    placing y nodes quadratically in xd...
    placing y nodes uniformly past xd to tsub...

    placing x nodes non-uniformly in the channel...
    setting xnodes =  32
    setting ynodes =  78
    maximum hx ratio =  1.111 from:
             hx max. = 0.1962 (um)
             hx min. = 0.0404 (um)
    maximum hy ratio =  3.252 from:
             hy max. = 0.0181 (um)
             hy min. = 0.0051 (um)


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
 Iter. #  1 took   0.06 msec. with norm = 0.800     at V(  1,  6).
 Iter. #  2 took   0.07 msec. with norm = 0.203     at V(  1, 21).
 Iter. #  3 took   0.06 msec. with norm = 0.106E-01 at V(  1, 21).
 Iter. #  4 took   0.06 msec. with norm = 0.284E-04 at V(  1, 21).
 Iter. #  5 took   0.06 msec. with norm = 0.215E-09 at V(  1, 21).
    The Newton solution converged in    5 iterations which took  0.00 sec.

 Solving for the exact n and V for the left boundary...
 Iter. #  1 took    0.4 msec. with norm = 0.177E-05 at n(  1, 53).
 Iter. #  2 took    0.2 msec. with norm = 0.582E-09 at n(  1, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.

 Solving the electrostatics for the right boundary...
 Iter. #  1 took   0.05 msec. with norm = 0.810     at V( 32,  6).
 Iter. #  2 took   0.06 msec. with norm = 0.205     at V( 32, 21).
 Iter. #  3 took   0.07 msec. with norm = 0.108E-01 at V( 32, 21).
 Iter. #  4 took   0.07 msec. with norm = 0.295E-04 at V( 32, 21).
 Iter. #  5 took   0.06 msec. with norm = 0.233E-09 at V( 32, 21).
    The Newton solution converged in    5 iterations which took  0.00 sec.

 Solving for the exact n and V for the right boundary...
 Iter. #  1 took    0.2 msec. with norm = 0.141E-05 at n( 32, 19).
 Iter. #  2 took    0.2 msec. with norm = 0.118E-11 at n( 32, 21).
    The Newton solution converged in   2 iterations which took   0.0 sec.

 Solving the electrostatic problem for slice i=  2.
 Iter. #  1 took   0.05 msec. with norm = 0.594E-09 at V(  2,  8).
    The Newton solution converged in    1 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  2.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  2, 66).
 Iter. #  2 took    0.2 msec. with norm = 0.581E-09 at n(  2, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  3.
 Iter. #  1 took   0.04 msec. with norm = 0.594E-09 at V(  3,  8).
    The Newton solution converged in    1 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  3.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  3, 53).
 Iter. #  2 took    0.2 msec. with norm = 0.582E-09 at n(  3, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  4.
 Iter. #  1 took   0.05 msec. with norm = 0.238E-04 at V(  4,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.402E-11 at V(  4, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  4.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  4, 66).
 Iter. #  2 took    0.2 msec. with norm = 0.580E-09 at n(  4, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  5.
 Iter. #  1 took   0.05 msec. with norm = 0.305E-03 at V(  5,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.666E-09 at V(  5, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  5.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  5, 73).
 Iter. #  2 took    0.2 msec. with norm = 0.569E-09 at n(  5, 41).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  6.
 Iter. #  1 took   0.05 msec. with norm = 0.332E-03 at V(  6,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.795E-09 at V(  6, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  6.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  6, 63).
 Iter. #  2 took    0.2 msec. with norm = 0.558E-09 at n(  6, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  7.
 Iter. #  1 took   0.04 msec. with norm = 0.376E-03 at V(  7,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.103E-08 at V(  7, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  7.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  7, 66).
 Iter. #  2 took    0.2 msec. with norm = 0.543E-09 at n(  7, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  8.
 Iter. #  1 took   0.05 msec. with norm = 0.413E-03 at V(  8,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.125E-08 at V(  8, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  8.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  8, 54).
 Iter. #  2 took    0.2 msec. with norm = 0.528E-09 at n(  8, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i=  9.
 Iter. #  1 took   0.05 msec. with norm = 0.460E-03 at V(  9,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.156E-08 at V(  9, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i=  9.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n(  9, 55).
 Iter. #  2 took    0.2 msec. with norm = 0.512E-09 at n(  9, 40).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 10.
 Iter. #  1 took   0.04 msec. with norm = 0.514E-03 at V( 10,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.197E-08 at V( 10, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 10.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 10, 66).
 Iter. #  2 took    0.2 msec. with norm = 0.493E-09 at n( 10, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 11.
 Iter. #  1 took   0.05 msec. with norm = 0.569E-03 at V( 11,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.245E-08 at V( 11, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 11.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 11, 76).
 Iter. #  2 took    0.2 msec. with norm = 0.473E-09 at n( 11, 76).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 12.
 Iter. #  1 took   0.05 msec. with norm = 0.631E-03 at V( 12,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.305E-08 at V( 12, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 12.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 12, 62).
 Iter. #  2 took    0.2 msec. with norm = 0.449E-09 at n( 12, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 13.
 Iter. #  1 took   0.05 msec. with norm = 0.705E-03 at V( 13,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.385E-08 at V( 13, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 13.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 13, 60).
 Iter. #  2 took    0.2 msec. with norm = 0.423E-09 at n( 13, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 14.
 Iter. #  1 took   0.05 msec. with norm = 0.781E-03 at V( 14,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.480E-08 at V( 14, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 14.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 14, 51).
 Iter. #  2 took    0.2 msec. with norm = 0.394E-09 at n( 14, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 15.
 Iter. #  1 took   0.05 msec. with norm = 0.870E-03 at V( 15,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.605E-08 at V( 15, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 15.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 15, 57).
 Iter. #  2 took    0.2 msec. with norm = 0.362E-09 at n( 15, 36).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 16.
 Iter. #  1 took   0.05 msec. with norm = 0.966E-03 at V( 16,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.759E-08 at V( 16, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 16.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 16, 57).
 Iter. #  2 took    0.2 msec. with norm = 0.326E-09 at n( 16, 39).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 17.
 Iter. #  1 took   0.05 msec. with norm = 0.107E-02 at V( 17,  7).
 Iter. #  2 took   0.05 msec. with norm = 0.958E-08 at V( 17, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 17.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 17, 74).
 Iter. #  2 took    0.2 msec. with norm = 0.286E-09 at n( 17, 40).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 18.
 Iter. #  1 took   0.05 msec. with norm = 0.968E-03 at V( 18,  7).
 Iter. #  2 took   0.09 msec. with norm = 0.791E-08 at V( 18, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 18.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 18, 49).
 Iter. #  2 took    0.2 msec. with norm = 0.251E-09 at n( 18, 40).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 19.
 Iter. #  1 took   0.05 msec. with norm = 0.872E-03 at V( 19,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.653E-08 at V( 19, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 19.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 19, 51).
 Iter. #  2 took    0.2 msec. with norm = 0.218E-09 at n( 19, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 20.
 Iter. #  1 took   0.05 msec. with norm = 0.786E-03 at V( 20,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.537E-08 at V( 20, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 20.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 20, 59).
 Iter. #  2 took    0.2 msec. with norm = 0.189E-09 at n( 20, 35).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 21.
 Iter. #  1 took   0.05 msec. with norm = 0.707E-03 at V( 21,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.440E-08 at V( 21, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 21.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 21, 53).
 Iter. #  2 took    0.2 msec. with norm = 0.162E-09 at n( 21, 36).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 22.
 Iter. #  1 took   0.05 msec. with norm = 0.636E-03 at V( 22,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.360E-08 at V( 22, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 22.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 22, 67).
 Iter. #  2 took    0.2 msec. with norm = 0.139E-09 at n( 22, 38).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 23.
 Iter. #  1 took   0.05 msec. with norm = 0.573E-03 at V( 23,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.295E-08 at V( 23, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 23.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 23, 49).
 Iter. #  2 took    0.2 msec. with norm = 0.117E-09 at n( 23, 35).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 24.
 Iter. #  1 took   0.04 msec. with norm = 0.516E-03 at V( 24,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.241E-08 at V( 24, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 24.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 24, 50).
 Iter. #  2 took    0.2 msec. with norm = 0.979E-10 at n( 24, 36).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 25.
 Iter. #  1 took   0.05 msec. with norm = 0.465E-03 at V( 25,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.197E-08 at V( 25, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 25.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 25, 61).
 Iter. #  2 took    0.2 msec. with norm = 0.806E-10 at n( 25, 35).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 26.
 Iter. #  1 took   0.05 msec. with norm = 0.418E-03 at V( 26,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.161E-08 at V( 26, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 26.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 26, 53).
 Iter. #  2 took    0.2 msec. with norm = 0.651E-10 at n( 26, 76).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 27.
 Iter. #  1 took   0.04 msec. with norm = 0.377E-03 at V( 27,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.131E-08 at V( 27, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 27.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 27, 56).
 Iter. #  2 took    0.2 msec. with norm = 0.510E-10 at n( 27, 37).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 28.
 Iter. #  1 took   0.04 msec. with norm = 0.339E-03 at V( 28,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.107E-08 at V( 28, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 28.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 28, 63).
 Iter. #  2 took    0.2 msec. with norm = 0.377E-10 at n( 28, 34).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 29.
 Iter. #  1 took   0.05 msec. with norm = 0.305E-03 at V( 29,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.870E-09 at V( 29, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 29.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 29, 51).
 Iter. #  2 took    0.2 msec. with norm = 0.265E-10 at n( 29, 35).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 30.
 Iter. #  1 took   0.04 msec. with norm = 0.275E-03 at V( 30,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.709E-09 at V( 30, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 30.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 30, 51).
 Iter. #  2 took    0.2 msec. with norm = 0.158E-10 at n( 30, 33).
    The Newton solution converged in   2 iterations which took   0.0 sec.
 Solving the electrostatic problem for slice i= 31.
 Iter. #  1 took   0.05 msec. with norm = 0.247E-03 at V( 31,  7).
 Iter. #  2 took   0.06 msec. with norm = 0.577E-09 at V( 31, 20).
    The Newton solution converged in    2 iterations which took  0.00 sec.
 Solving the exact n and V problem for slice i= 31.
 Iter. #  1 took    0.2 msec. with norm = 0.177E-05 at n( 31, 75).
 Iter. #  2 took    0.2 msec. with norm = 0.691E-11 at n( 31, 32).
    The Newton solution converged in   2 iterations which took   0.0 sec.

  +------> 1 x nodes/character         maximum voltage =   -9.18 kT/q =  -0.238 V
  !                    at  32, 14 going v the local dV =    2.32 kT/q =   0.060 V
  !                                        elec. conc. = 0.20E+05 nio = 0.30E+15 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V

  gggggggggggggggggggggggggggggggg
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
 s                                d
 s                                d
 s                                d
 s                                d
 s                                d
 s                                d
 s                                d
 s                                d
 s                                d
 s88888888888888888888888888888888d
 s99999999999999999999999999999999d
 sppppppppppppppppppppppppppppppppd
 sppppppppppppppppppppppppppppppppd
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp




 Now attempting to use DC() to solve for the biases:
   Vs =   0.000   Vd =   0.350   Vg1 =   2.500   Eavg =     0.0
 Solving for t=  0.000     sec. with dt= 0.1000E+33 sec.
 Iter. #  1 took  0:00+ 0:00+ 0:00 with norm =  7.77     at n( 32, 19).
 Iter. #  2 took  0:00+ 0:00+ 0:00 with norm =  5.11     at n( 32, 23).
 Iter. #  3 took  0:00+ 0:00+ 0:00 with norm =  9.71     at n( 32, 24).
 Iter. #  4 took  0:00+ 0:00+ 0:00 with norm =  4.45     at n(  2, 23).
 Iter. #  5 took  0:00+ 0:00+ 0:00 with norm =  3.08     at V( 32, 23).
 Iter. #  6 took  0:00+ 0:00+ 0:00 with norm =  5.09     at n( 27,  7).
 Iter. #  7 took  0:00+ 0:00+ 0:00 with norm =  9.43     at n( 27, 16).
 Iter. #  8 took  0:00+ 0:00+ 0:00 with norm =  3.00     at n( 16, 67).
 Iter. #  9 took  0:00+ 0:00+ 0:00 with norm =  9.98     at n(  6, 19).
 Iter. # 10 took  0:00+ 0:00+ 0:00 with norm =  19.2     at n(  5, 18).
 Iter. # 11 took  0:00+ 0:00+ 0:00 with norm =  22.6     at n(  4, 18).
 Iter. # 12 took  0:00+ 0:00+ 0:00 with norm =  34.9     at n( 30,  7).
 Iter. # 13 took  0:00+ 0:00+ 0:00 with norm =  29.6     at n( 30, 12).
 Iter. # 14 took  0:00+ 0:00+ 0:00 with norm =  25.5     at n(  3, 18).
 Iter. # 15 took  0:00+ 0:00+ 0:00 with norm =  25.3     at n(  3, 19).
 Iter. # 16 took  0:00+ 0:00+ 0:00 with norm =  28.5     at n( 30, 13).
 Iter. # 17 took  0:00+ 0:00+ 0:00 with norm =  29.0     at n( 30, 14).
 Iter. # 18 took  0:00+ 0:00+ 0:00 with norm =  29.9     at n( 31, 14).
 Iter. # 19 took  0:00+ 0:00+ 0:00 with norm =  31.9     at n( 31, 15).
 Iter. # 20 took  0:00+ 0:00+ 0:00 with norm =  26.0     at n( 30, 16).

  +------> 1 x nodes/character         maximum voltage =    9.73 kT/q =   0.252 V
  !                    at  31, 28 going > the local dV =   11.26 kT/q =   0.292 V
  !                                        elec. conc. = 0.16E+08 nio = 0.25E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V

  gggggggggggggggggggggggggggggggg
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
 sNNNN--                        -Nd
 sNN--                          -Nd
 sN-                             nd
 sN                              -d
 s-                              -d
 s-                               d
 s                                d
 s                                d
 s                                d
 s          88888888888           d
 s      8999999999999999998       d
 s     8ppppppppppppppppppp9      d
 s     9pppppppppppppppppppp8     d
      8ppppppppppppppppppppp9
      9pppppppppppppppppppppp8
     9ppppppppppppppppppppppp9
  889ppppppppppppppppppppppppp9
  ppppppppppppppppppppppppppppp988
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp
  pppppppppppppppppppppppppppppppp



 Adding a node to the left of ( 32, 28) because:
   Er=0.10469E+06 V/cm  vr=-.71552E+08 cm/s  v/E= -683.4 cm^2/V-s  hr=0.40E-01 um
   Dr=  18.1 cm^2/sec  2D/v=0.50676E-02 um  h*E=  16.325 kT  dV=  15.939 kT
 A node was added at the right side, so now xnodes= 33 and ynodes= 78.
 Iter. # 21 took  0:00+ 0:00+ 0:00 with norm =  32.8     at n( 31, 16).
 Iter. # 22 took  0:00+ 0:00+ 0:00 with norm =  4.09     at n( 27, 22).
 Iter. # 23 took  0:00+ 0:00+ 0:00 with norm =  2.28     at n( 27, 22).
 Iter. # 24 took  0:00+ 0:00+ 0:00 with norm =  2.39     at n( 28, 26).
 Iter. # 25 took  0:00+ 0:00+ 0:00 with norm =  2.97     at n( 28, 26).
 Iter. # 26 took  0:00+ 0:00+ 0:00 with norm =  5.93     at n( 30, 30).
 Iter. # 27 took  0:00+ 0:00+ 0:00 with norm =  16.4     at n( 32, 30).
 Iter. # 28 took  0:00+ 0:00+ 0:00 with norm =  6.56     at n( 31, 33).
 Iter. # 29 took  0:00+ 0:00+ 0:00 with norm =  4.67     at n( 32, 34).
 Iter. # 30 took  0:00+ 0:00+ 0:00 with norm =  4.68     at n( 32, 33).
 Iter. # 31 took  0:00+ 0:00+ 0:00 with norm =  5.79     at n( 32, 32).
 Iter. # 32 took  0:00+ 0:00+ 0:00 with norm =  13.6     at n( 32, 31).
 Iter. # 33 took  0:00+ 0:00+ 0:00 with norm =  7.48     at n( 32, 31).
 Iter. # 34 took  0:00+ 0:00+ 0:00 with norm = 0.614E-01 at n( 32, 32).
 Iter. # 35 took  0:00+ 0:00+ 0:00 with norm = 0.437E-05 at V( 33,  7).
 Iter. # 36 took  0:00+ 0:00+ 0:00 with norm = 0.182E-05 at V( 33,  7).
 Iter. # 37 took  0:00+ 0:00+ 0:00 with norm = 0.753E-06 at V( 33,  7).
    The Newton solution converged in   37 iterations which took   0:02.


 ==== Shark   v2.0                                                                   REB 1988         charge map ====

  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region



  +------> 1 x nodes/character         maximum voltage =   10.41 kT/q =   0.270 V
  !                    at  33, 31 going v the local dV =   10.85 kT/q =   0.281 V
  !                                        elec. conc. = 0.16E+08 nio = 0.25E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V

  ggggggggggggggggggggggggggggggggg
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
 sNNNN--                       -NNNd
 sNN--                          -NNd
 sN-                            --Nd
 sN                              -Nd
 s-                              -Nd
 s-                               Nd
 s                                -d
 s                                -d
 s                                -d
 s          88888888888           -d
 s      8999999999999999998        d
 s     8ppppppppppppppppppp8       d
 s     9ppppppppppppppppppp9       d
      8ppppppppppppppppppppp8
      9ppppppppppppppppppppp9
     9ppppppppppppppppppppppp8
  889pppppppppppppppppppppppp98
  pppppppppppppppppppppppppppp988
  pppppppppppppppppppppppppppppp999
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp





 ==== Shark   v2.0                                                                  REB 1988 starting transient ====

 Outputting the simulation results for t =   0.00 (psec).
 Iter. #  1 took  0:00+ 0:00+ 0:00 with norm = 0.311E-06 at V( 33,  7).
    The Newton solution converged in    1 iterations which took   0:00.


 ==== Shark   v2.0                                                                   REB 1988         charge map ====

  The following is the charge map for equilibrium:
     g = resistive gate
     i = insulator, or oxide
     n = n-type
     p = p-type
     blank = depletion region



  +------> 1 x nodes/character         maximum voltage =   10.41 kT/q =   0.270 V
  !                    at  33, 31 going v the local dV =   10.85 kT/q =   0.281 V
  !                                        elec. conc. = 0.16E+08 nio = 0.25E+18 /cm**3
  ! 2 y nodes/character                minimum voltage =  -34.42 kT/q =  -0.892 V
  !                                        elec. conc. = 0.67E-11 nio = 0.10     /cm**3
  V

  ggggggggggggggggggggggggggggggggg
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
 sNNNN--                       -NNNd
 sNN--                          -NNd
 sN-                            --Nd
 sN                              -Nd
 s-                              -Nd
 s-                               Nd
 s                                -d
 s                                -d
 s                                -d
 s          88888888888           -d
 s      8999999999999999998        d
 s     8ppppppppppppppppppp8       d
 s     9ppppppppppppppppppp9       d
      8ppppppppppppppppppppp8
      9ppppppppppppppppppppp9
     9ppppppppppppppppppppppp8
  889pppppppppppppppppppppppp98
  pppppppppppppppppppppppppppp988
  pppppppppppppppppppppppppppppp999
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp
  ppppppppppppppppppppppppppppppppp



 Outputting the simulation results for t = ****** (psec).
  Total CPU time used:  0 days, 00:00:02
```

## Development

The `Makefile` is set up for a typical MacOS box, but there are versions of
the makefile for a lot of different platforms. I just kept making new ones
as I moved to new systems.

```bash
$ make
```

will build everything, and then you're ready to run the code.

### Necessary Tools

**[Homebrew](http://brew.sh/)** - all of the following tools can be installed with Homebrew. If it's not already installed on you laptop, it's easy enough to go to the website, run the command to install it, and then continue to the next step. If you are on linux, expect to be installing some scripts and RPMs, but everything is available there as well - just not as convenient.

A `gfortran` or `f77` compiler. I have installed with:
```bash
$ brew install gfortran
```
but it seems to be now part of the GCC collection, and is installed with:
```bash
$ brew install gcc
```
but you're going to get a lot more than just `gfortran`.

