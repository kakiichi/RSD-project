# For GCC Suite
FC = gfortran
#FC = ifort
FFLAGS = -fdefault-real-8 -O3
#FFLAGS = -g -pg -Wall -fbounds-check -fdefault-real-8
#FFLAGS = -pg
#FFLAGS = -O3


SRC = param.f90 \
      main.f90 \
      absorber.f90 \
      corrfunc.f90 \
      interp.f90 \
      tau_eff.f90 \
      quadpack.f90

OBJS = $(SRC:.f90=.o)

IGMtransmission: $(OBJS) 
	  $(FC) $(FFLAGS) -o IGMtransmission.exe *.o

%.mod: %.o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	-rm $(OBJS)
	-rm *~                                                                       

