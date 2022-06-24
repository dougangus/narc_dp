################################################################################
# 
#   MAKEFILE for NARC propagator  
#
#   2007
################################################################################

### CONFIGURATION  

FC = gfortran
FFLAGS = -O2

#
# Source objects
#
OBJS= narc.o extrcart.o propagat.o inputfrq.o incidfrq.o eigenpbl.o fd_ana.o ehomogen.o crystal.o elastic.o inhomg_dis.o inhomg_ana.o utils.o
#
# Executable name
#
EXEC=../narc
#
#narc: $(EXEC)
$(EXEC): ${OBJS}
	$(FC) $(INCS) $(OPT) -o ${EXEC} ${OBJS}
## rules for compiling .f files:

.f.o:
	$(FC) $(FFLAGS) -c $*.f


clean:
	rm -f *.o *.mod
