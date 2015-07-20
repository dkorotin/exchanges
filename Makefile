FC = gfortran 
PROG = exchanges.x 
FFLAGS = -Og
LIBS = -llapack -lblas
LFLAGS = 

OBJ = parameters.o general.o iomodule.o find_nnbrs.o

all: $(PROG)
 
exchanges.x:  $(OBJ) exchanges.o
	$(FC) $(LFLAGS) -o $@ exchanges.o $(OBJ) $(LIBS)

clean:
	rm *.o
	
%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(*F).o $<
	
	
