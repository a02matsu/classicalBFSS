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
complex(kind(0d0)), allocatable :: Amat(:,:,:)

integer :: i,j,m,n,k

call set_parameters
call set_initial(Xmat,Vmat,Amat)
call calc_force(Amat,Xmat)

open(unit=Xmat_FILE, file=Xmat_FILE_NAME, status='replace', action='write')
open(unit=Vmat_FILE, file=Vmat_FILE_NAME, status='replace', action='write')
do k=0,Ntau-1
  time=time+deltaT
  call time_evolution_LeapFrog(Xmat,Vmat,Amat)

  call write_matrix(Xmat,Xmat_FILE)
  call write_matrix(Vmat,Vmat_FILE)
enddo
close(Xmat_FILE)
close(Vmat_FILE)

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


