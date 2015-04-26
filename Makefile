FC = gfortran 
PROG = aleip.x 
FFLAGS = -Og
LIBS = -llapack -lblas
LFLAGS = 

OBJ =  kinds.o read_hamilt.o compute_delta.o compute_g.o compute_gloc.o \
       read_pwscf.o 

all: $(PROG)
 
aleip.x:  $(OBJ) aleip.o
	$(FC) $(LFLAGS) -o $@ aleip.o $(OBJ) $(LIBS)

clean:
	rm *.o
	
%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(*F).o $<
	
	
