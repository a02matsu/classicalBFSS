module subroutines
use global_parameters
use matrix_functions
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine set_parameters
use global_parameters
implicit none

!!!!! read parameter file !!!!
open(PAR_FILE, file=PAR_FILE_NAME, status='old', action='READ')

read(PAR_FILE,*) NMAT
read(PAR_FILE,*) job_number
read(PAR_FILE,*) new_config
read(PAR_FILE,*) totalT
read(PAR_FILE,*) deltaT

close(PAR_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NTAU=totalT/deltaT

if( job_number == 0 ) then
  Inconf_FILE_NAME="CONFIG/lastconf.dat"
else
  write(Inconf_FILE_NAME, '("CONFIG/config_",i4.4)') job_number-1
endif

end subroutine set_parameters


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initial values
subroutine set_initial(Xmat,Vmat,Amat)
use global_parameters
implicit none

complex(kind(0d0)), allocatable :: Xmat(:,:,:)
complex(kind(0d0)), allocatable :: Vmat(:,:,:)
complex(kind(0d0)), allocatable :: Amat(:,:,:)
integer :: n, i,j
double precision :: r

integer :: seedsize
integer, allocatable :: seed(:)

double precision :: rmat(1:NMAT,1:NMAT,1:DIM)
double precision :: trace


! ディレクトリ生成
call system('./make_directory.sh')

! set matrices
allocate( Xmat(1:NMAT,1:NMAT,1:DIM) )
allocate( Vmat(1:NMAT,1:NMAT,1:DIM) )
allocate( Amat(1:NMAT,1:NMAT,1:DIM) )

if( new_config == 1 ) then 
  job_number=0
  time=0d0

  ! random numberの初期設定
  call random_seed(size=seedsize)
  allocate( seed(seedsize) )
  do i=1,seedsize
    call system_clock(count=seed(i))
  enddo
  call random_seed(put=seed(:))

  call random_number(rmat)
  rmat = 2d0*rmat - 1d0

  ! make rmat traceless
  do n=1,DIM
    trace=0d0
    do i=1,NMAT
      trace=trace+rmat(i,i,n)
    enddo
    do i=1,NMAT
      rmat(i,i,n)=rmat(i,i,n)-trace/dble(NMAT)
    enddo
  enddo

  ! from rmat to Xmat
  do n=1,DIM
    do i=1,NMAT
      Xmat(i,i,n)=dcmplx(rmat(i,i,n))
      do j=i+1,NMAT
        XMAT(i,j,n)=dcmplx(rmat(i,j,n)) + dcmplx(rmat(j,i,n))*(0d0,1d0)
        XMAT(j,i,n)=dcmplx(rmat(i,j,n)) - dcmplx(rmat(j,i,n))*(0d0,1d0)
      enddo
    enddo
  enddo

  !!!!!!!!!!!!!!!!!
  !! TEMPORARY
  !! Vmat must satisry [X_M, V_M]=0 
  Vmat=Xmat

else
  if( job_number== 0 ) then 
    Inconf_FILE_NAME="CONFIG/lastconf"
  else
    write(Inconf_FILE_NAME,'("CONFIG/lastconfig_",i4.4)') job_number-1
  endif

  open(unit=Inconf_FILE, file=Inconf_FILE_NAME, status='old', action='read',form='unformatted')

  read(Inconf_File) job_number
  job_number=job_number+1

  read(Inconf_File) time
  read(Inconf_File) Xmat
  read(Inconf_File) Vmat

  close(Inconf_FILE)
endif


write(Outconf_FILE_NAME,'("CONFIG/lastconfig_",i4.4)') job_number
write(Xmat_FILE_NAME,'("OUTPUT/Xmat_",i4.4)') job_number
write(Vmat_FILE_NAME,'("OUTPUT/Vmat_",i4.4)') job_number
end subroutine set_initial



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate forces
subroutine calc_force(Amat,Xmat)
use global_parameters
use matrix_functions, only : Matrix_Commutator
implicit none

complex(kind(0d0)), intent(out) :: Amat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)), intent(in) :: Xmat(1:NMAT,1:NMAT,1:DIM)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) tmpmat2(1:NMAT,1:NMAT)

integer :: n,m,i,j

Amat=(0d0,0d0)
do m=1,DIM
  do n=1,DIM
    if ( n .ne. m ) then
      call Matrix_commutator(tmpmat,Xmat(:,:,m),Xmat(:,:,n))
      call Matrix_commutator(tmpmat2,Xmat(:,:,n),tmpmat)
      Amat(:,:,m)=Amat(:,:,m)+tmpmat2
    endif
  enddo
enddo

end subroutine calc_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2nd order time evolution
subroutine time_evolution_LeapFrog(Xmat,Vmat,Amat)
use global_parameters
implicit none

complex(kind(0d0)), intent(inout) :: Xmat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)), intent(inout) :: Vmat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)), intent(inout) :: Amat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)) :: Amat2(1:NMAT,1:NMAT,1:DIM)
integer n

Xmat = Xmat + deltaT*VMat + dcmplx(0.5d0*deltaT*deltaT)*Amat

call calc_force(Amat2,Xmat)
Vmat = Vmat + dcmplx(0.5d0*deltaT)*(Amat+Amat2)
Amat=Amat2

end subroutine time_evolution_LeapFrog

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write matrix to file
subroutine write_matrix(MAT,F_NUM)
use global_parameters
implicit none

complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT,1:DIM)
integer, intent(in) :: F_NUM

integer m,j,i

write(F_NUM,'(E15.8,2X)',advance='no') time 
do m=1,dim
  do j=1,NMAT
    do i=1,NMAT
      write(F_NUM,'(E12.5,2X,E12.5,2x)',advance='no') MAT(i,j,m)
    enddo
  enddo
enddo
write(F_NUM,*)

end subroutine write_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read matrix from file
subroutine read_matrix(t, MAT,F_NUM)
use global_parameters
implicit none

double precision, intent(out) :: t
complex(kind(0d0)), intent(out) :: MAT(1:NMAT,1:NMAT,1:DIM)
integer, intent(in) :: F_NUM
integer m,j,i

double precision :: VAL(2*NMAT*NMAT*DIM+1)

read(F_NUM,*) VAL
t=VAL(1)
do m=1,dim
  do j=1,NMAT
    do i=1,NMAT
      MAT(i,j,m)= &
        dcmplx(VAL(2*NMAT*NMAT*(m-1)+2*NMAT*(j-1)+2*i)) &
        + dcmplx(VAL(2*NMAT*NMAT*(m-1)+2*NMAT*(j-1)+2*i+1)) * (0d0,1d0)
    enddo
  enddo
enddo



end subroutine read_matrix


end module subroutines
