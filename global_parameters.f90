!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! global parameters
module global_parameters
implicit none

integer :: check_gauss
integer :: check_ham

integer :: NMAT
integer, parameter :: DIM=9
!double precision :: time

integer :: job_number 
integer :: new_config
double precision :: deltaT
double precision :: totalT
integer :: write_output

double precision :: MASS2
integer :: Ntau  ! nint(totalT / deltaT)
integer :: dimG  ! NMAT*NMAT-1
integer :: matrix_size ! dimG * DIM


character(10), parameter :: FMT_time="E15.8"
character(10), parameter :: FMT_vals="E12.5"


!character(128) :: INPUT_FILE_NAME="input.dat"
character(128) :: INPUT_FILE_NAME="input.dat"
character(128) :: PAR_FILE_NAME="theory_parameters.dat"
character(128) :: Inconf_FILE_NAME
character(128) :: Outconf_FILE_NAME
character(128) :: Xmat_FILE_NAME
character(128) :: Vmat_FILE_NAME
character(128) :: Fmat_FILE_NAME
character(128) :: Xmode_FILE_NAME
character(128) :: Vmode_FILE_NAME
character(128) :: Fmode_FILE_NAME

!integer, parameter :: INPUT_FILE=10
integer, parameter :: INPUT_FILE=10
integer, parameter :: PAR_FILE=11
integer, parameter :: Inconf_FILE=12
integer, parameter :: Outconf_FILE=13
integer, parameter :: Xmat_FILE=14
integer, parameter :: Vmat_FILE=15
integer, parameter :: Fmat_FILE=16
integer, parameter :: Xmode_FILE=17
integer, parameter :: Vmode_FILE=18
integer, parameter :: Fmode_FILE=19



end module global_parameters


