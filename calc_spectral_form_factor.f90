!!!!!!!!!!!!!!!!!!
!! call as 
!! % calc_spectral_form_factor.exe [SV file name(input)] [SFF file name(output)] Dtau, ini_tau, fin_tau
program main
use global_parameters
use subroutines
implicit none

double precision :: Dtau !=1d-1
integer :: Ndiv
double precision :: ini_tau !=0d0
double precision :: fin_tau !=10d0
character(20) :: c_Ndiv
character(20) :: c_ini_tau
character(20) :: c_fin_tau

double precision, allocatable :: SV(:)
double precision :: sum_exp
character(128) :: IN_FILE
character(128) :: OUT_FILE

double precision :: tau
integer :: ite, i
complex(kind(0d0)) :: tmp


!! determine the file name including singular values
call getarg(1,IN_FILE) ! input singular values
call getarg(2,OUT_FILE) ! output spectral form factors
call getarg(3,c_Ndiv) ! output spectral form factors
call getarg(4,c_ini_tau) ! output spectral form factors
call getarg(5,c_fin_tau) ! output spectral form factors

read(c_Ndiv,*) Ndiv
read(c_ini_tau,*) ini_tau
read(c_fin_tau,*) fin_tau

Dtau = (fin_tau - ini_tau) / dble(Ndiv) 

!! set parameters
call set_parameters
matrix_size = NMAT*NMAT*DIM
allocate( SV(1:matrix_size) )

!! set the singular values
open(10,file=IN_FILE,status='old')
read(10,*) SV
close(10)
!sum_exp=0d0
!do i=1,matrix_size
  !sum_exp = sum_exp + dlog( sv(i) )
!enddo


!! g
open(unit=11, file=OUT_FILE, status='replace', action='write')
tau=ini_tau
do ite = 0, Ndiv
  tmp=(0d0,0d0)
  do i=1, matrix_size
    tmp = tmp + exp( (0d0,1d0) * dcmplx((SV(i))*tau) )
    !tmp = tmp + exp( (0d0,1d0) * dcmplx(dlog(SV(i))*tau) )
  enddo
  write(11, *) tau, dble(tmp*dconjg(tmp))/dble(matrix_size**2)
  !tmp=cdexp( (0d0,1d0) * sum_exp * tau )
  !write(11, *) tau, dble(tmp*dconjg(tmp))/dble(matrix_size**2)
  tau=tau+Dtau
enddo

end program main
