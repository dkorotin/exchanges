FC = mpif90  
#FC = gfortran  
PROG = aleip 
FFLAGS = -traceback
#FFLAGS = -traceback -check all
#LIBS = 
LFLAGS = -mkl=sequential
LFLAGS = -llapack
PREFLAGS = 
#
OBJ =  kinds.o read_hamilt.o compute_delta.o compute_g.o compute_gloc.o \
       read_pwscf.o 
#
all: $(PROG)
# 
aleip:  $(LIBS) aleip.o $(OBJ)
	$(FC) $(LFLAGS) -o $@ aleip.o $(OBJ) $(LIBS)
	strip $@
#
#
clean:
	rm *.o
	
%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(*F).o $<
	
	
