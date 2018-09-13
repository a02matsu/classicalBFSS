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
read(PAR_FILE,*) MASS2

read(PAR_FILE,*) job_number
read(PAR_FILE,*) new_config
read(PAR_FILE,*) write_output
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
subroutine set_initial(Xmat,Vmat,Fmat)
use global_parameters
implicit none

complex(kind(0d0)), allocatable :: Xmat(:,:,:)
complex(kind(0d0)), allocatable :: Vmat(:,:,:)
complex(kind(0d0)), allocatable :: Fmat(:,:,:)
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
allocate( Fmat(1:NMAT,1:NMAT,1:DIM) )

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

  !call random_number(rmat)
  !rmat = 2d0*rmat - 1d0

  call BoxMuller(rmat)

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
  !Vmat=Xmat
  Vmat=(0d0,0d0)

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
write(Fmat_FILE_NAME,'("OUTPUT/Fmat_",i4.4)') job_number
write(Xmode_FILE_NAME,'("OUTPUT/Xmode_",i4.4)') job_number
write(Vmode_FILE_NAME,'("OUTPUT/Vmode_",i4.4)') job_number
write(Fmode_FILE_NAME,'("OUTPUT/Fmode_",i4.4)') job_number
end subroutine set_initial



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate forces
subroutine calc_force(Fmat,Xmat)
use global_parameters
use matrix_functions, only : Matrix_Commutator
implicit none

complex(kind(0d0)), intent(out) :: Fmat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)), intent(in) :: Xmat(1:NMAT,1:NMAT,1:DIM)

complex(kind(0d0)) tmpmat(1:NMAT,1:NMAT)
complex(kind(0d0)) tmpmat2(1:NMAT,1:NMAT)

integer :: n,m,i,j

Fmat=(0d0,0d0)
do m=1,DIM
  do n=1,DIM
    if ( n .ne. m ) then
      call Matrix_commutator(tmpmat,Xmat(:,:,m),Xmat(:,:,n))
      call Matrix_commutator(tmpmat2,Xmat(:,:,n),tmpmat)
      Fmat(:,:,m)=Fmat(:,:,m)+tmpmat2
    endif
  enddo
  Fmat(:,:,m)=Fmat(:,:,m) - dcmplx( MASS2 )*Xmat(:,:,m)
enddo

end subroutine calc_force


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2nd order time evolution
subroutine time_evolution_LeapFrog(Xmat,Vmat,Fmat)
use global_parameters
implicit none

complex(kind(0d0)), intent(inout) :: Xmat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)), intent(inout) :: Vmat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)), intent(inout) :: Fmat(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)) :: Fmat2(1:NMAT,1:NMAT,1:DIM)
integer n

Xmat = Xmat + deltaT*VMat + dcmplx(0.5d0*deltaT*deltaT)*Fmat

call calc_force(Fmat2,Xmat)
Vmat = Vmat + dcmplx(0.5d0*deltaT)*(Fmat+Fmat2)
Fmat=Fmat2

end subroutine time_evolution_LeapFrog

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write matrix to file
subroutine write_matrix(MAT,FILE_NUM)
use global_parameters
implicit none

complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT,1:DIM)
integer, intent(in) :: FILE_NUM

integer m,j,i

write(FILE_NUM,'(E15.8,2X)',advance='no') time 
do m=1,dim
  do j=1,NMAT
    do i=1,NMAT
      write(FILE_NUM,'(E12.5,2X,E12.5,2x)',advance='no') MAT(i,j,m)
    enddo
  enddo
enddo
write(FILE_NUM,*)

end subroutine write_matrix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read matrix from file
subroutine read_matrix(t, MAT,FILE_NUM)
use global_parameters
implicit none

double precision, intent(out) :: t
complex(kind(0d0)), intent(out) :: MAT(1:NMAT,1:NMAT,1:DIM)
integer, intent(in) :: FILE_NUM
integer m,j,i

double precision :: VAL(2*NMAT*NMAT*DIM+1)

read(FILE_NUM,*) VAL
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! write modes to file
!!  The argument is MATRIX. The subroutine computes modes and 
!!  write it to the file. 
subroutine write_modes(MAT,FILE_NUM)
use global_parameters
use SUN_generators, only : trace_MTa
implicit none

complex(kind(0d0)), intent(in) :: MAT(1:NMAT,1:NMAT,1:DIM)
integer, intent(in) :: FILE_NUM
complex(kind(0d0)) :: trace

integer m,a

write(FILE_NUM,'(E15.8,2X)',advance='no') time 
do m=1,dim
  do a=1,NMAT*NMAT-1
    call trace_MTa(trace,MAT(:,:,m),a)
    write(FILE_NUM,'(E12.5,2X)',advance='no') dble(trace)
  enddo
enddo
write(FILE_NUM,*)

end subroutine write_modes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Gaussian random number
!!  BoxMuller(gauss)
!!
!! generage ensemble exp(-1/2(x^2) and exp(-1/2(y^2))
!! 
!! output 2N gaussian randum numbers
!! gauss is an allocatable array.
!! It will be reallocated after calling this routine.
subroutine BoxMuller(gauss)
implicit none

double precision, intent(out) :: gauss(:,:,:)
double precision, parameter :: PI=dacos(-1d0)
double precision, allocatable :: rand(:)
double precision, allocatable :: tmp_gauss(:)

integer :: K,L,M,N
integer i,a,b,c

K=size(gauss,1)
L=size(gauss,2)
M=size(gauss,3)

N=K*L*M
if( mod(N,2) == 0 ) then
  allocate( rand(1:N) )
  allocate( tmp_gauss(1:N) )
else
  allocate( rand(1:N+1) )
  allocate( tmp_gauss(1:N+1) )
endif

call random_number( rand )
!call genrand_real3(rand)
!write(*,*) rand

do i=1,N/2
  tmp_gauss(2*i-1) = dsqrt(-2d0*dlog(rand(2*i-1)))*dsin(2d0*Pi*rand(2*i))
  tmp_gauss(2*i) = dsqrt(-2d0*dlog(rand(2*i-1)))*dcos(2d0*Pi*rand(2*i))
enddo
if( mod(N,2) /= 0 ) then
  tmp_gauss(N) = dsqrt(-2d0*dlog(rand(N)))*dcos(2d0*Pi*rand(N+1))
endif

do c=1,M
  do b=1,L
    do a=1,K
      i = K*L*(c-1) + K*(b-1) + a
      gauss(a,b,c) = tmp_gauss(i)
    enddo
  enddo
enddo

end subroutine BoxMuller



end module subroutines
