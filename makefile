################################################################################
# 
#   MAKEFILE for NARC propagator  
#
#   2011
################################################################################


# Standard fortran compiler 
FC = gfortran

# Flags for compiling
#OPT = -fbounds-check -ffixed-line-length-none -g -Wall
OPT = -fbounds-check -ffixed-line-length-none -Wall
OPT77 =
OPT90 =

narc: narc.o extrcart.o propagat.o inputfrq.o incidfrq.o eigenpbl.o fd_ana.o \
	ehomogen.o crystal.o elastic.o inhomg_dis.o inhomg_ana.o utils.o extrcart_sac.o \
	squirt_lib.o computestiff.o VTI_lib.o cextrcart.o cextrcart_sac.o cpropagat.o ceigenpbl.o
	$(FC) $(OPT) -o ../Bin/narc narc.o extrcart.o propagat.o inputfrq.o incidfrq.o eigenpbl.o fd_ana.o \
	ehomogen.o crystal.o elastic.o inhomg_dis.o inhomg_ana.o utils.o extrcart_sac.o \
	squirt_lib.o computestiff.o VTI_lib.o cextrcart.o cextrcart_sac.o cpropagat.o ceigenpbl.o

clean:
	rm -f *.o


# 
#	Compile Instuctions
#
%.o: %.f90
	$(FC) $(OPT) $(OPT90) -c $(MODULES) $*.f90 
%.o: %.f
	$(FC) $(OPT) $(OPT77) -c $(MODULES) $*.f
