!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate interval average of 
!! M_{IJ} = 1/Delta( \int X_I(s) X_J(t+s) ds ) - 1/Delta( \int X_I(s)ds) 1/Delta(\int X_J(s) ds )
!! 
!! output is the eigenvalues of M_{IJ}
subroutine calc_eigenvalues_of_correlations(Eigenvalues, FILE1, FILE1_NAME, FILE2, FILE2_NAME, interval)
use global_parameters
use subroutines
use matrix_functions, only : matrix_eigenvalues
implicit none

complex(kind(0d0)), intent(out) :: Eigenvalues(1:matrix_size)
integer, intent(in) :: FILE1, FILE2
character(128), intent(in) :: FILE1_NAME, FILE2_NAME
double precision, intent(in) :: interval

complex(kind(0d0)) :: MAT1(1:NMAT,1:NMAT,1:DIM)
complex(kind(0d0)) :: MAT2(1:NMAT,1:NMAT,1:DIM)
double precision :: mode1(1:dimG, 1:DIM)
double precision :: mode2(1:dimG, 1:DIM)
double precision :: tm, re(1:2*NMAT*NMAT*DIM)
integer :: pos
integer :: num_data

complex(kind(0d0)) :: matrixM(1:matrix_size, 1:matrix_size)
double precision :: tmp_mat(1:matrix_size, 1:matrix_size)
double precision :: tmp_vec1(1:matrix_size)
double precision :: tmp_vec2(1:matrix_size)

character(50) :: FMT1, FMT2

integer :: a,b,i,j,m,n,d

!FMT1="'(" // trim(FMT_time) // ")'"
!FMT2="'(" // trim(FMT_vals) // ")'"

open(FILE1,file=FILE1_NAME,status='old')
open(FILE2,file=FILE2_NAME,status='old')

num_data = nint(interval/deltaT)
tmp_mat=0d0
tmp_vec1=0d0
tmp_vec2=0d0
do d=1,num_data
  read(FILE1,*) tm,re
  pos=1
  do m=1,DIM
    do j=1,NMAT
      do i=1,NMAT
        MAT1(i,j,m) = dcmplx(re(pos)) + (0d0,1d0)*dcmplx(re(pos+1))
        pos=pos+2
      enddo
    enddo
  enddo
  call matrix_to_modes(mode1,mat1)
  !!
  read(FILE2,*) tm,re
  pos=1
  do m=1,DIM
    do j=1,NMAT
      do i=1,NMAT
        MAT2(i,j,m) = dcmplx(re(pos)) + (0d0,1d0)*dcmplx(re(pos+1))
        pos=pos+2
      enddo
    enddo
  enddo
  call matrix_to_modes(mode2,mat2)
  !!!!!!!!!!!!!!
  do n=1,DIM
    do b=1,dimG
      do m=1,DIM
        do a=1,dimG
          tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) &
            = tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) + mode1(a,m)*mode2(b,n)*deltaT
        enddo
      enddo
    enddo
  enddo
  !!!
  do m=1,DIM
    do a=1,dimG
      tmp_vec1(dimG*(m-1)+a) = tmp_vec1(dimG*(m-1)+a) + mode1(a,m)*deltaT
      tmp_vec2(dimG*(m-1)+a) = tmp_vec2(dimG*(m-1)+a) + mode2(a,m)*deltaT
    enddo
  enddo
enddo

do n=1,DIM
  do b=1,dimG
    do m=1,DIM
      do a=1,dimG
        matrixM(dimG*(m-1)+a, dimG*(n-1)+b) = &
          dcmplx(tmp_mat(dimG*(m-1)+a, dimG*(n-1)+b) / interval) &
          - dcmplx(tmp_vec1(dimG*(m-1)+a)/interval * tmp_vec2(dimG*(n-1)+b)/interval)
      enddo
    enddo
  enddo
enddo
close( FILE1 )
close( FILE2 )

call matrix_eigenvalues(eigenvalues, matrixM)


end subroutine calc_eigenvalues_of_correlations



