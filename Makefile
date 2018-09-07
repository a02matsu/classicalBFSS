#FC=ifort
FC=gfortran
FLAGS_IFORT= -mkl -CB -traceback -g 
FLAGS_GCC= -I~/lib/ -O2 -llapack -lblas
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
OBS_SRC=calc_correlation.f90 
OBS_OBJ=$(OBS_SRC:.f90=.o)
OBS_PROG=calc_obs.exe
########################


#.SUFFIXES : .o .f90 # .oを作るときは必ず.f90から作るよ
.SUFFIXES : .f90 # .oを作るときは必ず.f90から作るよ
 
all:$(PROG) 

$(PROG): $(OBJS) $(MAIN_OBJ)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(MAIN_OBJ) 
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(MAIN_OBJ) 
endif


obs:$(OBS_PROG)

$(OBS_PROG): $(OBJS) $(OBS_OBJ)
ifeq ($(FC),gfortran)
	$(FC) -O2 $(FLAGS_GCC) -o $@ $(OBJS) $(OBS_OBJ) $(LIB)
else
	$(FC) $(FLAGS_IFORT) -o $@ $(OBJS) $(OBS_OBJ) $(LIB)
endif



# moduleをコンパイルするときの依存性を解消
%.o: %.f90
ifeq ($(FC),gfortran)
	$(FC) -c $<
else
	$(FC) $(FLAGS_IFORT) -c $<
endif
%.mod: %.f90 %.o
	@true

# moduleの依存性
subroutines.o: \
  global_parameters.o \
  matrix_functions.o


.PHONY: clean
clean:
	mv $(PROG) $(PROG).bak; rm -f *.o *.mod core 
