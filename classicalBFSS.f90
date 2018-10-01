!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAIN
program main
use global_parameters
use subroutines
use matrix_functions, only : check_hermitian
implicit none

!!! X(s), V(s), F(s)
complex(kind(0d0)), allocatable :: Xmat1(:,:,:)
complex(kind(0d0)), allocatable :: Vmat1(:,:,:)
complex(kind(0d0)), allocatable :: Fmat1(:,:,:)
!!! X(s+t), V(s+t), F(s+t)
complex(kind(0d0)), allocatable :: Xmat2(:,:,:)
complex(kind(0d0)), allocatable :: Vmat2(:,:,:)
complex(kind(0d0)), allocatable :: Fmat2(:,:,:)
!!! \int ds X(s), \int ds V(s), \int ds F(s), 
complex(kind(0d0)), allocatable :: Xmat1_int(:)
complex(kind(0d0)), allocatable :: Vmat1_int(:)
complex(kind(0d0)), allocatable :: Fmat1_int(:)
!!! \int ds X(s+t), \int ds V(s+t), \int ds F(s+t), 
complex(kind(0d0)), allocatable :: Xmat2_int(:)
complex(kind(0d0)), allocatable :: Vmat2_int(:)
complex(kind(0d0)), allocatable :: Fmat2_int(:)
!!! \int ds X(s)X(s+t),  \int ds V(s)V(s+t),  \int ds F(s)F(s+t), 
complex(kind(0d0)), allocatable :: XX_int(:,:)
complex(kind(0d0)), allocatable :: VV_int(:,:)
complex(kind(0d0)), allocatable :: FF_int(:,:)
double precision :: time, Ham0, Ham1

integer :: i,j,m,n,k,ii,jj,nn,ele_I,ele_J, counter
character(128) OUT_FMT

!! read parameter and input files
call set_parameters

! set matrices
allocate( Xmat1(1:NMAT,1:NMAT,1:DIM) )
allocate( Vmat1(1:NMAT,1:NMAT,1:DIM) )
allocate( Fmat1(1:NMAT,1:NMAT,1:DIM) )
allocate( Xmat2(1:NMAT,1:NMAT,1:DIM) )
allocate( Vmat2(1:NMAT,1:NMAT,1:DIM) )
allocate( Fmat2(1:NMAT,1:NMAT,1:DIM) )
!!!
allocate( Xmat1_int(1:matrix_size) )
allocate( Vmat1_int(1:matrix_size) )
allocate( Fmat1_int(1:matrix_size) )
allocate( Xmat2_int(1:matrix_size) )
allocate( Vmat2_int(1:matrix_size) )
allocate( Fmat2_int(1:matrix_size) )
!!!
allocate( XX_int(matrix_size, matrix_size) )
allocate( VV_int(matrix_size, matrix_size) )
allocate( FF_int(matrix_size, matrix_size) )

! set initial values of Xmat1 
! and 
! make Xmat2 by evoluting Xmat1 by dulationT 
  !open(unit=23, file="test", status='replace', action='write')
  !write(23,*) "###"
call set_initial(time,Xmat1,Vmat1,Fmat1,Xmat2,Vmat2,Fmat2)

if( check_ham == 1 ) then
  call check_hamiltonian(Xmat1,Vmat1)
endif

OUT_FMT=trim('(a,I3,2X,a,f6.2,2X,a,E12.2,2X,a,E12.2,2X,a,E12.2)')
if( write_output == 0 ) then
  open(unit=MXX_FILE, file=MXX_FILE_NAME, status='replace', action='write')
  open(unit=MVV_FILE, file=MVV_FILE_NAME, status='replace', action='write')
  open(unit=MFF_FILE, file=MFF_FILE_NAME, status='replace', action='write')
  write(MXX_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"totalT=",totalT, "integ_time=",integrationT, "delay_time=", dulationT
  write(MVV_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"totalT=",totalT, "integ_time=",integrationT, "delay_time=", dulationT
  write(MFF_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"totalT=",totalT, "integ_time=",integrationT, "delay_time=", dulationT
!  open(unit=Xmat_FILE, file=Xmat_FILE_NAME, status='replace', action='write')
!  open(unit=Vmat_FILE, file=Vmat_FILE_NAME, status='replace', action='write')
!  open(unit=Fmat_FILE, file=Fmat_FILE_NAME, status='replace', action='write')
!  !!!
!  open(unit=XXint_FILE, file=XXint_FILE_NAME, status='replace', action='write')
!  open(unit=VVint_FILE, file=VVint_FILE_NAME, status='replace', action='write')
!  open(unit=FFint_FILE, file=FFint_FILE_NAME, status='replace', action='write')
!  !!!
!  write(Xmat_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"TIME=",totalT,"deltaT=", deltaT, "integration_time=",integrationT
!  write(Vmat_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"TIME=",totalT,"deltaT=", deltaT, "integration_time=",integrationT
!  write(Fmat_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"TIME=",totalT,"deltaT=", deltaT, "integration_time=",integrationT
!  !!!
!  write(XXint_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"TIME=",totalT,"deltaT=", deltaT, "integration_time=",integrationT
!  write(VVint_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"TIME=",totalT,"deltaT=", deltaT, "integration_time=",integrationT
!  write(FFint_FILE,OUT_FMT) "# NMAT=",NMAT,"M=",Mass,"TIME=",totalT,"deltaT=", deltaT, "integration_time=",integrationT
endif

call calc_hamiltonian(Ham0,Xmat1,Vmat1)
Xmat1_int=(0d0,0d0)
Vmat1_int=(0d0,0d0)
Fmat1_int=(0d0,0d0)
Xmat2_int=(0d0,0d0)
Vmat2_int=(0d0,0d0)
Fmat2_int=(0d0,0d0)
XX_int=(0d0,0d0)
VV_int=(0d0,0d0)
FF_int=(0d0,0d0)
counter=0
do k=0,Ntau-1
  !!! integration 
  call integration(Xmat1_int, Xmat2_int, XX_int, Xmat1, Xmat2, counter)
  call integration(Vmat1_int, Vmat2_int, VV_int, Vmat1, Vmat2, counter)
  call integration(Fmat1_int, Fmat2_int, FF_int, Fmat1, Fmat2, counter)
  !!!! time evolution 
  time=time+deltaT 
  call time_evolution_LeapFrog(Xmat1,Vmat1,Fmat1)
  call time_evolution_LeapFrog(Xmat2,Vmat2,Fmat2)
  !! check gauss's law constraint
  if( check_gauss == 1 ) then
    call check_gauss_law(Xmat1,Vmat1,time)
  endif
  counter = counter + 1
  if( counter == num_integration ) then
    !! write out result
    !call check_hermitian( XX_int )
    !call check_hermitian( VV_int )
    !call check_hermitian( FF_int )
    if( write_output == 0 ) then 
      !call write_matrix(time,Xmat_int,Xmat_FILE)
      !call write_matrix(time,Vmat_int,Vmat_FILE)
      !call write_matrix(time,Fmat_int,Fmat_FILE)
      !!!!
      call write_matrix2(time, Xmat1_int, Xmat2_int, XX_int, MXX_FILE)
      call write_matrix2(time, Vmat1_int, Vmat2_int, VV_int, MVV_FILE)
      call write_matrix2(time, Fmat1_int, Fmat2_int, FF_int, MFF_FILE)
    endif
    !! reset integration data
    Xmat1_int=(0d0,0d0)
    Vmat1_int=(0d0,0d0)
    Fmat1_int=(0d0,0d0)
    Xmat2_int=(0d0,0d0)
    Vmat2_int=(0d0,0d0)
    Fmat2_int=(0d0,0d0)
    XX_int=(0d0,0d0)
    VV_int=(0d0,0d0)
    FF_int=(0d0,0d0)
    counter = 0
  endif
  !!!!!!!!!!!!!!!!!!!!
enddo
call calc_hamiltonian(Ham1,Xmat1,Vmat1)

write(*,*) "# H(final)-H(initial)=",Ham1-Ham0

if( write_output == 0 ) then
  close(MXX_FILE)
  close(MVV_FILE)
  close(MFF_FILE)
endif

open(unit=Outconf_FILE, file=Outconf_FILE_NAME, status='replace', form='unformatted')
write(Outconf_File) job_number
write(Outconf_FILE) time 
write(Outconf_FILE) Xmat1
write(Outconf_FILE) Vmat1
close(Outconf_FILE)

call system('&
  cd CONFIG; &
  FILE=$(ls lastconfig_* | tail -1); &
  LINK="lastconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
  fi; ln -s $FILE $LINK')


end program main


