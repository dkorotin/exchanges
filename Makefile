PROG = exchanges.x 

#gfortan Mac Os
# FC = gfortran 
# FFLAGS = -fopenmp -Og -O2 -g
# LIBS = -fopenmp -framework Accelerate #-L/usr/local/opt/lapack/lib/ -llapack -lblas

# Intel fortran linux
 FC = ifort 
 FFLAGS = -qopenmp -O2 -g -mavx
 LIBS = -qopenmp -lmkl_intel_lp64  -lmkl_sequential -lmkl_core
#uncomment for debug:
# FFLAGS = -qopenmp -O0 -g -mavx -traceback -check

# gfortran linux
# FC = gfortran 
# FFLAGS = -fopenmp -Og -O2 -g
# LIBS = -L/usr/local/opt/lapack/lib/ -llapack -lblas


LFLAGS =

OBJ = parameters.o general.o iomodule.o find_nnbrs.o green_function.o

all: $(PROG)
 
exchanges.x:  $(OBJ) exchanges.o
	$(FC) $(LFLAGS) -o $@ exchanges.o $(OBJ) $(LIBS)

clean:
	rm *.o *.mod
	
%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(*F).o $<
	
	
