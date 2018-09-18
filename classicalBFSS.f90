!!! module for global parameters
!include "global_parameters.f90"
!include "subroutines.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MAIN
program main
use global_parameters
use subroutines
implicit none

complex(kind(0d0)), allocatable :: Xmat(:,:,:)
complex(kind(0d0)), allocatable :: Vmat(:,:,:)
complex(kind(0d0)), allocatable :: Fmat(:,:,:)
double precision :: time

integer :: i,j,m,n,k

call set_parameters
call set_initial(time,Xmat,Vmat,Fmat)
call calc_force(Fmat,Xmat)

if( check_ham == 1 ) then
  call check_hamiltonian(Xmat,Vmat)
endif


if( write_output == 0 ) then
  open(unit=Xmat_FILE, file=Xmat_FILE_NAME, status='replace', action='write')
  open(unit=Vmat_FILE, file=Vmat_FILE_NAME, status='replace', action='write')
  open(unit=Fmat_FILE, file=Fmat_FILE_NAME, status='replace', action='write')
endif
do k=0,Ntau-1
  time=time+deltaT
  call time_evolution_LeapFrog(Xmat,Vmat,Fmat)
  if( check_gauss == 1 ) then
    call check_gauss_law(Xmat,Vmat,time)
  endif
  if( write_output == 0 ) then 
    call write_matrix(time,Xmat,Xmat_FILE)
    call write_matrix(time,Vmat,Vmat_FILE)
    call write_matrix(time,Fmat,Fmat_FILE)
  endif
enddo
if( write_output == 0 ) then
  close(Xmat_FILE)
  close(Vmat_FILE)
  close(Fmat_FILE)
endif

open(unit=Outconf_FILE, file=Outconf_FILE_NAME, status='replace', form='unformatted')
write(Outconf_File) job_number
write(Outconf_FILE) time 
write(Outconf_FILE) Xmat
write(Outconf_FILE) Vmat
close(Outconf_FILE)

call system('&
  cd CONFIG; &
  FILE=$(ls lastconfig_* | tail -1); &
  LINK="lastconf"; if [ -e "$LINK" ]; then /bin/rm $LINK; &
  fi; ln -s $FILE $LINK')


end program main


