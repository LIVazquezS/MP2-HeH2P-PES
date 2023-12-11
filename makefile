F90 = gfortran #CHOOSE YOUR FORTRAN COMPILER HERE

FFLAGS = -O3   #OPTIONAL: OPTIMIZATION FLAG
##########################
# Object Files for build #
##########################


OBJS = \
evaluatekernel.o \
RKHS.o \


 all: evaluatekernel.x clean
 evaluatekernel.x : $(OBJS1)
 	 ${F90}  -o $@ $(OBJS1)

#######################################
# Object dependencies and compilation #
#######################################
evaluatekernel.o : evaluate_kernel.f90 \
 RKHS.o
 	$(F90) -c $(FFLAGS) -o $@ evaluate_kernel.f90

RKHS.o : RKHS.f90
	$(F90) -c $(FFLAGS) -o $@ RKHS.f90

.PHONY: clean
clean:
	rm *.mod *.o