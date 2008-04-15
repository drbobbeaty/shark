#
#   This is the makefile for Shark for Mac OS X
#
#	XMAX = maximum number of nodes in the x direction
#	YMAX = maximum number of nodes in the y direction
#
CC = gcc
F77 = gfortran
ifeq ($(shell uname -m),Power Macintosh)
	CPU_OPTS = -faltivec -ftree-vectorize
else
ifeq ($(shell uname -m),i386)
	CPU_OPTS = -mfpmath=sse -msse3
else
ifeq ($(shell uname -m),sun3m)
	F77 = f77
	CPU_OPTS = -f68881
endif
else
ifeq ($(shell uname -m),NP1)
	F77 = fort10.3
	CPU_OPTS = -DNP1
endif
endif
endif
CCOPTS = -O3 $(CPU_OPTS)
F77OPTS = -O3 -Wall $(CPU_OPTS)
LDOPTS = 
CTAGS = ectags  --excmd=number --tag-relative=no  --fields=+a+m+n+S -R
SIZE = -DXMAX=75 -DYMAX=100
#
#   These are all the non-optimized object files
#
F77_OBJS = shark.o depmap.o newton.o getdk.o mkdev.o mkmesh.o \
       parser.o guess.o bias.o sharkio.o equil.o DC.o \
       oneDes.o full1D.o orthomin.o minres.o gauss.o linpk.o \
       sgtsl.o sgbsl.o isamax.o

CC_OBJS = sscal.o sdot.o saxpy.o sgbfa.o

OBJS = $(F77_OBJS) $(CC_OBJS)

#   This is the general way to make a .o from a .f
#    -Wimplicit => make all undefined variables cause error messages
#    -c => do not link the output
.f.o :
	@ $(F77) $(F77OPTS) -c $*.f -o $*.o

#   This is the general way to make a .o from a .c
.c.o :
	@ $(CC) $(CCOPTS) -c $*.c -o $*.o

#
#   The -DDIRECT means we have enough memory for the
#   direct solvers. -DLAPACK means that we're using the
#   LAPACK solvers.
#
ifeq ($(shell uname -m),NP1)
#   Notice the -DNP1 in the compile line!
.F.o :
	@ $(F77) -DNP1 -F $(SIZE) $*.F
	@ mv $*.f $*.fs
	@ sed 's/^#/c/' < $*.fs > $*.f
	@ $(F77) -w77 -r8 -i4 -u -c $*.f
	@ rm $*.f $*.fs
else
.F.o :
	$(F77) $(F77OPTS) -DDIRECT -c $(SIZE) $*.F -o $*.o
endif

#
#   This is the final compilation of all the modules
#   The '-framework Accelerate' is needed for LAPACK
#
shark : $(OBJS)
	$(F77) $(CPU_OPTS) $(OBJS) -o shark $(LDOPTS)

clean :
	rm -rf shark $(OBJS) $(OPT_OBJS) core

tags :
	$(CTAGS) `pwd`

#
# These are the custom compiles that are needed for the gfortran 4.3.0
#
#newton.o : newton.F
#	$(F77) -O3 -Wall -ftree-vectorize -DDIRECT -c $(SIZE) $*.F -o $*.o

mkmesh.o : mkmesh.F
	$(F77) $(F77OPTS) -DDIRECT -c $(SIZE) $*.F -o $*.o

ifeq ($(shell uname -m),NP1)
#   These routines can be compiled with optimization
minres.o orthomin.o :
	@ $(F77) -DNP1 -F $(SIZE) $*.F
	@ mv $*.f $*.fs
	@ sed 's/^#/c/' < $*.fs > $*.f
	@ $(F77) -w77 -r8 -i4 -u -O -c $*.f
	@ rm $*.f $*.fs

#   These routines can be compiled with vectorization and optimization
linpk.o gauss.o sgtsl.o :
	@ $(F77) -DNP1 -F $(SIZE) $*.F
	@ mv $*.f $*.fs
	@ sed 's/^#/c/' < $*.fs > $*.f
	@ $(F77) -w77 -r8 -i4 -u -V -O -c $*.f
	@ rm $*.f $*.fs $*.fv
endif
