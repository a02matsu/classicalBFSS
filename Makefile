#FC=ifort
FC=gfortran
#FLAGS_IFORT= -mkl -CB -traceback -g 
FLAGS_IFORT= -mkl -O2
FLAGS_GCC= -O2 -llapack -lblas 
# コンパイルのために順番が大事。下層ほど先に書く。 
SRCS=\
     global_parameters.f90 \
     matrix_functions.f90 \
     SUN_generators.f90 \
     subroutines.f90
OBJS=$(SRCS:.f90=.o)
MAIN_SRCS=classicalBFSS.f90
MAIN_OBJ=classicalBFSS.o
PROG=classicalBFSS.exe
#########################
CORR_SRC=calc_correlation.f90 
CORR_OBJ=$(CORR_SRC:.f90=.o)
CORR_PROG=calc_correlations.exe
########################
SFF_SRC=calc_spectral_form_factor.f90 
SFF_OBJ=$(SFF_SRC:.f90=.o)
SFF_PROG=calc_spectral_form_factor.exe
#########################
LIBS= ~/lib/liblapack.a


#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) 

$(PROG): $(OBJS) $(MAIN_OBJ)
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -o $@ $(OBJS) $(MAIN_OBJ) $(LIBS)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(MAIN_OBJ) 
endif


corr:$(CORR_PROG)

$(CORR_PROG): $(OBJS) $(CORR_OBJ)
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -o $@ $(OBJS) $(CORR_OBJ) $(LIBS)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(CORR_OBJ)
endif


sff:$(SFF_PROG)

$(SFF_PROG): $(OBJS) $(SFF_OBJ)
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -o $@ $(OBJS) $(SFF_OBJ) $(LIBS)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(SFF_OBJ)
endif

# moduleをコンパイルするときの依存性を解消
%.o: %.f90
ifeq ($(FC),gfortran)
	$(FC) $(FLAGS_GCC) -c $<
else
	$(FC) $(FLAGS_IFORT) -c $<
endif
%.mod: %.f90 %.o
	@true

# moduleの依存性
subroutines.o: \
  global_parameters.o \
  SUN_generators.o \
  matrix_functions.o
  calc_correlation.o: \
  matrix_functions.o
calc_correlation.o: \
	global_parameters.o \
	subroutines.o \
	matrix_functions.o 
calc_spectral_form_factor.o: \
	global_parameters.o \
	subroutines.o




.PHONY: clean
clean:
	mv $(PROG) $(PROG).bak; rm -f *.o *.mod core 
