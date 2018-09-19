program calc_correlation
use global_parameters
use subroutines
use matrix_functions 
implicit none

character(128) :: arg
character :: op_XVF ! X,V,F
character :: op_GM ! gauge or matrix
character :: op_DN ! <X_i^dagger X_j> or <X_i X_j>
double precision :: time
integer :: Num_Lines

character(128) :: DELAY_FILE_NAME
integer, parameter :: DELAY_FILE=100

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision, parameter :: Delta = 1.0d0 ! \Delta
double precision, parameter :: DULATION = 0.0d0 ! t
integer :: NUM_Delta
integer :: NUM_DULATION 
integer :: NUM_SAMPLES
integer :: INI1, FIN1, INI2, FIN2
integer :: counter, traj, k, i

integer, parameter :: EIGEN=102
character(128) :: EIGEN_NAME
character(20) :: FMT
complex(kind(0d0)), allocatable :: eigenvalues(:)
double precision, allocatable :: singularvalues(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(kind(0d0)), allocatable :: mode1(:,:)
complex(kind(0d0)), allocatable :: mode2(:,:)
complex(kind(0d0)), allocatable :: matrixM(:,:) !(1:matrix_size, 1:matrix_size)
complex(kind(0d0)), allocatable :: tmp_mat(:,:) !(1:matrix_size, 1:matrix_size)
complex(kind(0d0)), allocatable :: tmp_vec1(:) ! (1:matrix_size)
complex(kind(0d0)), allocatable :: tmp_vec2(:) !(1:matrix_size)
double precision :: rate
!!!!!!!!!!

integer t1,t2,t_rate, t_max


!! specify the output files
call getarg(1,arg)
call getarg(2,op_XVF)
call getarg(3,op_GM)
call getarg(4,op_DN)

if( op_XVF == "X" ) then 
  Xmat_FILE_NAME="OUTPUT/Xmat_" // arg
elseif( op_XVF == "V" ) then 
  Xmat_FILE_NAME="OUTPUT/Vmat_" // arg
elseif( op_XVF == "F" ) then 
  Xmat_FILE_NAME="OUTPUT/Fmat_" // arg
else
  write(*,*) "give me X or V or F"
  stop
endif
!!
if( op_GM /= "G" .and. op_GM /= "M" ) then
  write(*,*) "selsect gauge mode (G) or matrix mode (M)"
  stop
endif

if( op_DN /= "D" .and. op_DN /= "N" ) then
  write(*,*) "Warnning: settig to <X_i^\dagger X_j> mode"
  op_DN = "D" 
endif

EIGEN_NAME = "EIGENS/" // op_XVF // op_GM // op_DN // "_SV_" // arg
DELAY_FILE_NAME = "EIGENS/" // op_XVF // op_GM // op_DN // "tmp3.dat"
FMT='(' // FMT_vals // ',2X)'

!! read theory data from theory_parameters.dat
open(PAR_FILE, file=PAR_FILE_NAME, status='old', action='READ')
read(PAR_FILE,*) NMAT
read(PAR_FILE,*) MASS2
read(PAR_FILE,*) temperature
read(PAR_FILE,*) deltaT
close(PAR_FILE)

!! set variables
if( op_GM == "G" ) then 
  dimG = NMAT*NMAT-1
elseif( op_GM == "M" ) then
  dimG = NMAT*NMAT
endif
matrix_size = dimG * DIM


allocate( mode1(1:dimG,1:DIM) )
allocate( mode2(1:dimG,1:DIM) )
allocate( eigenvalues(1:matrix_size) )
allocate( singularvalues(1:matrix_size) )
allocate( matrixM(1:matrix_size, 1:matrix_size) )
allocate( tmp_mat(1:matrix_size, 1:matrix_size) )
allocate( tmp_vec1(1:matrix_size) )
allocate( tmp_vec2(1:matrix_size) )


!! count number of lines
!!  and
!! prepare shfted data
NUM_DULATION = nint( DULATION/deltaT )
NUM_Delta = nint( Delta/deltaT )
!call system_clock(t1)
!write(*,*) NUM_DULATION, NUM_Delta
call make_delay_file(Num_Lines, Delay_FILE, Delay_FILE_NAME, &
  Xmat_FILE, Xmat_FILE_NAME, NUM_DULATION )
!call make_delay_file_sys(Num_Lines, Delay_FILE_NAME, &
  !Xmat_FILE_NAME, NUM_DULATION )
NUM_SAMPLES = (Num_lines - NUM_DULATION) / NUM_Delta
!call system_clock(t2,t_rate,t_max)
!if( t2<t1) then
  !write(*,*) dble((t_max-t1) + t2 + 1)/dble(t_rate)
!else
  !write(*,*) dble(t2-t1) / dble(t_rate)
!endif




!! open input files
open(Xmat_FILE,file=Xmat_FILE_NAME,status='old')
open(Delay_FILE,file=Delay_FILE_NAME,status='old')
!! open output file
open(unit=EIGEN, file=EIGEN_NAME, status='replace', action='write')

write(EIGEN,'(a,I6,2X,a,E10.3,2X,a,E10.3)') "# #(samples) = ", NUM_SAMPLES, "Delta = ", Delta, "t=", DULATION
do traj=1, NUM_SAMPLES
  tmp_mat=(0d0,0d0)
  tmp_vec1=(0d0,0d0)
  tmp_vec2=(0d0,0d0)
  do k=1, NUM_Delta
    if( OP_GM == "G" ) then 
      call read_modes(mode1,Xmat_FILE)
      call read_modes(mode2,Delay_FILE) 
    elseif( op_GM =="M" ) then
      call read_mat(mode1,Xmat_FILE)
      call read_mat(mode2,Delay_FILE) 
    endif
    !write(*,*) mode1
    !! 
    if(k==1 .or. k==NUM_Delta) then
      rate=0.5d0
    else
      rate=1d0
    endif
    call integration_step(tmp_mat,tmp_vec1,tmp_vec2,mode1,mode2,rate,op_DN)
    !call calc_eigenvalues_of_correlations(eigenvalues,mode1,mode2)
    !call check_hermitian(tmp_mat)
  enddo
  tmp_mat = tmp_mat / Delta
  tmp_vec1 = tmp_vec1 / Delta
  tmp_vec2 = tmp_vec2 / Delta
  !call calc_eigenvalues(eigenvalues, tmp_mat, tmp_vec1, tmp_vec2)
  !write(EIGEN,*) dble(eigenvalues)
  call calc_singularvalues(singularvalues, tmp_mat, tmp_vec1, tmp_vec2)
  do i=1,matrix_size
    write(EIGEN,FMT,advance='no') singularvalues(i)
  enddo
  write(EIGEN,*) 
enddo

close(Xmat_FILE)
close(Delay_FILE)
close(EIGEN)

call system( '/bin/rm ' // trim(DELAY_FILE_NAME) )


end program calc_correlation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Count number of lines of a given file
subroutine count_lines(num_lines, NFILE,NFILE_NAME)
implicit none

integer, intent(out) :: num_lines
integer, intent(in) :: NFILE
character(128), intent(in) :: NFILE_NAME

open(NFILE,file=NFILE_NAME,status='old')
num_lines=0

do 
  read(NFILE,'()',end=100)
  num_lines = num_lines + 1
enddo

100 close(NFILE)

end subroutine count_lines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Count number of lines of a given file
subroutine count_lines_sys(num_lines, NFILE_NAME)
use global_parameters
implicit none

integer, intent(out) :: num_lines
character(128), intent(in) :: NFILE_NAME

call system( 'cat ' // trim(NFILE_NAME) // ' | grep "" -c > ' // trim(TMP_FILE_NAME) )
open(TMP_FILE,file=TMP_FILE_NAME,status='old')
read(TMP_FILE,*) num_lines
close(TMP_FILE)
call system( 'rm ' // trim(TMP_FILE_NAME) )

end subroutine count_lines_sys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_delay_file(NLINE, OUT_FILE, OUT_FILE_NAME, IN_FILE, IN_FILE_NAME, num_delay)
implicit none

integer, intent(out) :: NLINE
integer, intent(in) :: IN_FILE, OUT_FILE, num_delay
character(128), intent(in) :: IN_FILE_NAME, OUT_FILE_NAME
integer,parameter :: max_line_len = 10000
character(max_line_len) linebuf
integer :: i

NLINE=0
open(IN_FILE,file=IN_FILE_NAME,status='old')
open(unit=OUT_FILE, file=OUT_FILE_NAME, status='replace', action='write')

do i=1,num_delay
  read(IN_FILE,'()',end=100)
  NLINE=NLINE+1
enddo

do 
  read ( IN_FILE, '(a)', end=100) linebuf
  NLINE=NLINE+1
  write( OUT_FILE, '(a)' ) trim(linebuf)
enddo

100 close(IN_FILE)
close(OUT_FILE)


end subroutine make_delay_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_delay_file_sys(num_lines, OUT_FILE_NAME, IN_FILE_NAME, num_delay)
use global_parameters
implicit none

integer, intent(out) :: NUM_LINES
integer, intent(in) :: num_delay
character(128), intent(in) :: IN_FILE_NAME, OUT_FILE_NAME
character(50) :: c_num
integer :: i


call system( 'cat ' // trim(IN_FILE_NAME) // ' | grep "" -c > ' // trim(TMP_FILE_NAME) )
open(TMP_FILE,file=TMP_FILE_NAME,status='old')
read(TMP_FILE,*) num_lines
close(TMP_FILE)
call system( 'rm ' // trim(TMP_FILE_NAME) )

write(c_num,*) num_lines-num_delay
call system( 'tail -n  ' // trim(c_num) // "  " // &
  trim(IN_FILE_NAME) // ' > ' // trim(OUT_FILE_NAME) )


end subroutine make_delay_file_sys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_partial_data(OUT_FILE, OUT_FILE_NAME, INI, FIN, IN_FILE, IN_FILE_NAME)
implicit none

integer, intent(in) :: INI, FIN
integer, intent(in) :: IN_FILE, OUT_FILE
character(128), intent(in) :: IN_FILE_NAME, OUT_FILE_NAME

integer,parameter :: max_line_len = 10000
character(max_line_len) linebuf
integer :: i

open(IN_FILE,file=IN_FILE_NAME,status='old')
open(unit=OUT_FILE, file=OUT_FILE_NAME, status='replace', action='write')


do i=1, INI-1
  read( IN_FILE, '()' )
enddo
do i=INI,FIN
  read ( IN_FILE, '(a)' ) linebuf
  write( OUT_FILE, '(a)' ) trim(linebuf)
enddo

close( OUT_FILE )
close( IN_FILE )


end subroutine make_partial_data

!!!!!!!!!!!!!!!!
subroutine read_modes(mode,NFILE)
use global_parameters
use subroutines, only : matrix_to_modes
use matrix_functions, only :  check_hermitian
implicit none

complex(kind(0d0)), intent(out) :: mode(1:dimG,1:DIM)
integer, intent(in) :: NFILE

complex(kind(0d0)) :: mat(1:NMAT,1:NMAT,1:DIM)
integer :: m,j,i,pos
double precision :: tm, re(1:2*NMAT*NMAT*DIM)

read(NFILE,*) tm,re
pos=1
do m=1,DIM
  do j=1,NMAT
    do i=1,NMAT
      MAT(i,j,m) = dcmplx(re(pos)) + (0d0,1d0)*dcmplx(re(pos+1))
      pos=pos+2
    enddo
  enddo
enddo
!do m=1,DIM
  !call check_hermitian(mat(:,:,m))
!enddo
call matrix_to_modes(mode,mat)
!mode=dble(c_mode)
end subroutine read_modes

!!!!!!!!!!!!!!!!
subroutine read_mat(mat,NFILE)
use global_parameters
use subroutines, only : matrix_to_modes
use matrix_functions, only :  check_hermitian
implicit none

complex(kind(0d0)) :: mat(1:dimG,1:DIM)
integer, intent(in) :: NFILE

integer :: m,j,i,pos
double precision :: tm, re(1:2*NMAT*NMAT*DIM)

read(NFILE,*) tm,re
pos=1
do m=1,DIM
  do j=1,NMAT
    do i=1,NMAT
      MAT(NMAT*(j-1)+i,m) = dcmplx(re(pos)) + (0d0,1d0)*dcmplx(re(pos+1))
      pos=pos+2
    enddo
  enddo
enddo

end subroutine read_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integration_step(tmp_mat,tmp_vec1,tmp_vec2,mode1,mode2,rate,op_DN)
use global_parameters
implicit none
complex(kind(0d0)), intent(out) :: tmp_mat(1:matrix_size, 1:matrix_size)
complex(kind(0d0)), intent(out) :: tmp_vec1(1:matrix_size)
complex(kind(0d0)), intent(out) :: tmp_vec2(1:matrix_size)
complex(kind(0d0)), intent(in) :: mode1(1:dimG,1:DIM)
complex(kind(0d0)), intent(in) :: mode2(1:dimG,1:DIM)
double precision, intent(in) :: rate
character, intent(in) :: op_DN

integer :: a,b,m,n

do n=1,DIM
  do b=1,dimG
    do m=1,DIM
      do a=1,dimG
        if( op_DN == "D" ) then 
          tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) &
          = tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) &
            + dconjg(mode1(a,m))*mode2(b,n)*dcmplx(deltaT*rate)
        else
          tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) &
          = tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) &
            + mode1(a,m)*mode2(b,n)*dcmplx(deltaT*rate)
        endif
      enddo
    enddo
  enddo
enddo
!!!
do m=1,DIM
  do a=1,dimG
    if( op_DN == "D" ) then 
      tmp_vec1(dimG*(m-1)+a) = tmp_vec1(dimG*(m-1)+a) + dconjg(mode1(a,m))*deltaT*rate
    else
      tmp_vec1(dimG*(m-1)+a) = tmp_vec1(dimG*(m-1)+a) + mode1(a,m)*deltaT*rate
    endif
    tmp_vec2(dimG*(m-1)+a) = tmp_vec2(dimG*(m-1)+a) + mode2(a,m)*deltaT*rate
  enddo
enddo


end subroutine integration_step

!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_eigenvalues(eigenvalues, tmp_mat, tmp_vec1, tmp_vec2)
use global_parameters
use matrix_functions, only : matrix_eigenvalues
implicit none

complex(kind(0d0)), intent(out) :: Eigenvalues(1:matrix_size)
complex(kind(0d0)), intent(in) :: tmp_mat(1:matrix_size, 1:matrix_size)
complex(kind(0d0)), intent(in) :: tmp_vec1(1:matrix_size)
complex(kind(0d0)), intent(in) :: tmp_vec2(1:matrix_size)

complex(kind(0d0)) :: matrixM(1:matrix_size, 1:matrix_size)

integer :: a,b,m,n

do n=1,DIM
  do b=1,dimG
    do m=1,DIM
      do a=1,dimG
        matrixM(dimG*(m-1)+a, dimG*(n-1)+b) = &
          tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) &
          - tmp_vec1(dimG*(m-1)+a) * tmp_vec2(dimG*(n-1)+b)
      enddo
    enddo
  enddo
enddo
call matrix_eigenvalues(eigenvalues, matrixM)
!write(*,*) eigenvalues


end subroutine calc_eigenvalues


!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_singularvalues(singularvalues, tmp_mat, tmp_vec1, tmp_vec2)
use global_parameters
use matrix_functions, only : matrix_singularvalues
implicit none

double precision, intent(out) :: singularvalues(1:matrix_size)
complex(kind(0d0)), intent(in) :: tmp_mat(1:matrix_size, 1:matrix_size)
complex(kind(0d0)), intent(in) :: tmp_vec1(1:matrix_size)
complex(kind(0d0)), intent(in) :: tmp_vec2(1:matrix_size)

complex(kind(0d0)) :: matrixM(1:matrix_size, 1:matrix_size)

integer :: a,b,m,n

do n=1,DIM
  do b=1,dimG
    do m=1,DIM
      do a=1,dimG
        matrixM(dimG*(m-1)+a, dimG*(n-1)+b) = &
          tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) &
          - tmp_vec1(dimG*(m-1)+a) * tmp_vec2(dimG*(n-1)+b)
      enddo
    enddo
  enddo
enddo
call matrix_singularvalues(singularvalues, matrixM)
!write(*,*) eigenvalues


end subroutine calc_singularvalues

!include "calc_eigenvalues_of_correlations.f90"



