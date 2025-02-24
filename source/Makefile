FC=ifort
# FC=ifx
# FC=gfortran


FLAGS = -O3 -lfftw3 -llapack -lblas #-fbounds-check -DDEBUG -pg         
# FLAGS = -O3 -lfftw3 -llapack -lblas #-fbounds-check -DDEBUG -pg
# FFLAGS = -O3 -I$(FFTW_INCLUDE) -lfftw3 -llapack -lblas    # Para ifort
# FFLAGS = -O3 -I$(FFTW_INCLUDE) -lfftw3 -llapack -lblas    # Para gfortran

SRC = global_module.f90 \
      modelvar_module.f90 \
      parse_module.f90 \
      sysdef_module.f90 \
      initial_module.f90 \
      models_module.f90 \
      rpmd_module.f90 \
      print_module.f90 \
      propagation_module.f90 \
      nonadiabatic_module.f90 \
      zhunakamura_module.f90 \
      runtraj_module.f90

MAIN =  rpmdsh.f90

XOBJ = ../bin/sharp.x

OBJS =  ${SRC:.f90=.o}


%.o : %.f90
	$(FC) -c $(FLAGS) $< -o $@

%.mod : %.f90
	$(FC) -c $(FLAGS) $< -o $@


all:  $(OBJS)
	$(FC) $(OBJS) $(MAIN) -o $(XOBJ) $(FLAGS)
#	$(FC) $(OBJS) $(MAIN) -o $(XOBJ) $(FLAGS)  $(INCLUDE) $(LIB)

clean:
	@echo "==> Cleaning ALL intermediate and executable files"
	rm -f *.mod *.o $(XOBJ)
