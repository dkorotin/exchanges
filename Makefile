PROG = exchanges.x 

#gfortan Mac Os
FC = gfortran 
FFLAGS = -fopenmp -Og -O2 -g
LIBS = -fopenmp -framework Accelerate #-L/usr/local/opt/lapack/lib/ -llapack -lblas

# Intel fortran linux
# FC = ifort 
# FFLAGS = -openmp -Og -O2 -g
# LIBS = -openmp -L/usr/local/opt/lapack/lib/ -llapack -lblas

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
	
	
