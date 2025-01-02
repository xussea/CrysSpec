subroutine init_TDM(TDM,TDMe,N,Ne,maxv,maxvct,NB1,NB2,NB3,NB,B1,B2,B3,BCT,IB1,IB2,IB3,IBCT,SF,SC,SA)
!use types
implicit none
include "include/types.inc"
integer :: i,j,N1,M1,V1,U1
integer :: N,Ne,maxv,maxvct
integer :: NB1,NB2,NB3,NB
integer, dimension(1:N,1:N) :: IBCT
integer, dimension(1:N,0:maxv) :: IB1
integer, dimension(1:N,0:maxv,1:N,1:maxv) :: IB2
integer, dimension(1:N,0:maxvct,1:N,0:maxvct) :: IB3

double precision ::  SF,SC,SA
double precision, dimension(NB,3) :: TDM
double precision, dimension(Ne,3) :: TDMe

double precision, external :: FC

type(P1), dimension(1:N*(N-1)) :: BCT
type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2
type(P2), dimension(1:N*(maxvct+1)*(N-1)*(maxvct+1)) :: B3

TDM(:,:)=0.00
do i=1,NB1
  N1=B1(i)%N
  V1=B1(i)%V
  TDM(IB1(N1,V1),:)=TDMe(N1,:)*FC(0,V1,SF)
end do

!do i=1,NB3
!  N1=B3(i)%N
!  V1=B3(i)%V
!  M1=B3(i)%M
!  U1=B3(i)%U
!  TDM(IB3(N1,V1,M1,U1),:)=TDMe(N+IBCT(N1,M1),:)*FC(0,V1,SC)*FC(0,U1,SA)
!end do


end subroutine init_TDM

