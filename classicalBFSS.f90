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

integer :: i,j,m,n,k

call set_parameters
call set_initial(Xmat,Vmat,Fmat)
call calc_force(Fmat,Xmat)

if( write_output == 0 ) then
  open(unit=Xmat_FILE, file=Xmat_FILE_NAME, status='replace', action='write')
  open(unit=Vmat_FILE, file=Vmat_FILE_NAME, status='replace', action='write')
  open(unit=Fmat_FILE, file=Fmat_FILE_NAME, status='replace', action='write')
  open(unit=Xmode_FILE, file=Xmode_FILE_NAME, status='replace', action='write')
  open(unit=Vmode_FILE, file=Vmode_FILE_NAME, status='replace', action='write')
  open(unit=Fmode_FILE, file=Fmode_FILE_NAME, status='replace', action='write')
endif
do k=0,Ntau-1
  time=time+deltaT
  call time_evolution_LeapFrog(Xmat,Vmat,Fmat)
  if( write_output == 0 ) then 
    call write_matrix(Xmat,Xmat_FILE)
    call write_matrix(Vmat,Vmat_FILE)
    call write_matrix(Fmat,Fmat_FILE)
    call write_modes(Xmat,Xmode_FILE)
    call write_modes(Fmat,Fmode_FILE)
  endif
enddo
if( write_output == 0 ) then
  close(Xmat_FILE)
  close(Vmat_FILE)
  close(Fmat_FILE)
  close(Xmode_FILE)
  close(Vmode_FILE)
  close(Fmode_FILE)
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


