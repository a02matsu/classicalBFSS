program calc_correlation
use global_parameters
use subroutines
implicit none

character(128) :: arg
complex(kind(0d0)), allocatable :: Xmat(:,:,:)
complex(kind(0d0)), allocatable :: Vmat(:,:,:)


call getarg(1,arg)

Xmat_FILE_NAME="OUTPUT/Xmat_" // arg
Vmat_FILE_NAME="OUTPUT/Vmat_" // arg



open(PAR_FILE, file=PAR_FILE_NAME, status='old', action='READ')
read(PAR_FILE,*) NMAT
close(PAR_FILE)
allocate(XMAT(1:NMAT,1:NMAT,1:DIM))


open(Xmat_FILE, file=Xmat_FILE_NAME, status='old', action='READ')
open(Vmat_FILE, file=Vmat_FILE_NAME, status='old', action='READ')
call read_matrix(time,Xmat,Xmat_FILE)
call read_matrix(time,Vmat,Vmat_FILE)

write(*,*) Xmat

close(Xmat_FILE)

end program calc_correlation
