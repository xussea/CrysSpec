subroutine CD_abs(N,U,Um,C,E,Intensity,R,ECDvel)
implicit none
integer :: i,j,k,N
double precision, dimension(N,N) :: C,H
double precision, dimension(N) :: DiabE,E,Intensity
double precision, dimension(N,3) :: U,R
double precision, dimension(3) :: tmp
double precision, parameter :: cc=137.0359999
double precision, parameter :: au2ev=27.2116
double precision :: f
double precision, dimension(N) :: II
double precision, dimension(N,3) :: Um

logical :: ECDvel

if ( ECDvel .eqv. .false. ) then
  II=0.0
  do i=1,N
     do j=1,N
        do k=1,N
          f=0.5*E(i)/(au2ev*cc)
          call cross(R(k,:)-R(j,:),U(k,1:3),tmp)
          II(i)=II(i)+C(j,i)*C(k,i)*dot_product(U(i,:),Um(i,:))
          II(i)=II(i)+f*C(j,i)*C(k,i)*dot_product(U(j,1:3),tmp)
        end do
     end do
     Intensity(i)=II(i)/maxval(abs(II))
  end do
else

  H=0.0
  do i=1,N
    H(i,i)=E(i)
  end do  
  H=matmul(matmul(C,H),transpose(C))

  II=0.0
  do i=1,N
     do j=1,N
        do k=1,N
          f=H(j,j)/H(k,k)
          II(i)=II(i)+C(j,i)*C(k,i)*dot_product(U(i,:),Um(i,:))
        end do
     end do
     Intensity(i)=II(i)
  end do

end if

end subroutine CD_abs

subroutine create_R(R,R0,N,Ne,NB,NB1,NB2,NB3,B1,B2,B3,maxv,maxvct)
implicit none
include "include/types.inc"
integer :: i
integer :: maxv,maxvct,N,Ne,NB1,NB2,NB3,NB
integer, dimension(1:N,0:maxv) :: IB1
integer, dimension(1:N,0:maxv,1:N,1:maxv) :: IB2
integer, dimension(1:N,0:maxvct,1:N,0:maxvct) :: IB3

double precision, dimension(Ne,3) :: R0
double precision, dimension(NB,3) :: R

type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2
type(P2), dimension(1:N*(maxvct+1)*(N-1)*(maxvct+1)) :: B3

do i=1,NB1
        R(i,:)=R0(B1(i)%N,:)
end do

do i=1,NB2
        R(i+NB1,:)=R0(B2(i)%N,:)
end do
end subroutine create_R
