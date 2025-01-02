subroutine basisGS(N,maxv,maxvct,B1,B2,IB1,IB2,NB,NB1,NB2,twop)
!use types
implicit none
include "include/types.inc"
integer :: i,j,k,l
integer :: N,NB1,NB2,NB
integer :: maxv,maxvct

integer, dimension(1:N,0:maxv) :: IB1
integer, dimension(1:N,0:maxv,1:N,1:maxv) :: IB2

type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2

character(len=2), external :: lowcase
logical :: twop
!ONE-PARTICLE BASIS SET 
NB1=0
do i=1,N
  do j=1,maxv
    NB1=NB1+1
    B1(NB1)%N=i
    B1(NB1)%V=j
    IB1(i,j)=NB1
  end do
end do

!TWO-PARTICLE BASIS SET
NB2=0
if (twop .eqv. .true. ) then
do i=1,N
  do j=1,maxv
    do k=i+1,N
      if (i .ne. k) then
        do l=1,maxv-j
          NB2=NB2+1
          B2(NB2)%N=i
          B2(NB2)%V=j
          B2(NB2)%M=k
          B2(NB2)%U=l
          IB2(i,j,k,l)=NB2+NB1
        end do
      end if
    end do
  end do
end do
end if

NB=NB1+NB2
end subroutine basisGS

subroutine create_emisionTDM(N,Ne,maxv,maxvct,NB1,NB2,NB3,NB,NGS,NGS1,NGS2,B1,B2,B3,GS1,GS2,IBCT,BCT,ETDMD,TDMe,SF,SC,SA)
implicit none
include "include/types.inc"

integer :: i,j,k
integer :: N,Ne,NB1,NB2,NB3,NB,NGS,NGS1,NGS2,NGS3
integer :: maxv,maxvct
integer :: N1,N2,V1,V2,M1,M2,W1,W2

integer, dimension(1:N,1:N) :: IBCT

type(P1), dimension(1:N*(N-1)) :: BCT
type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2
type(P2), dimension(1:N*(maxvct+1)*(N-1)*(maxvct+1)) :: B3

type(P1), dimension(1:N*(maxv+1)) :: GS1
type(P2), dimension(1:N*(maxv+1*(N-1)*maxv)) :: GS2

double precision :: SF,SA,SC
double precision, dimension(Ne,3) :: TDME
double precision, dimension(0:NGS,1:NB,3) :: ETDMD
double precision, external :: FC

ETDMD(:,:,:) = 0.00d00
do i=1,NB
  do j=0,NGS

    !print*, NB,NGS,i,j
    if ( j .eq. 0 ) then
      if ( i .LE. NB1 ) then
        N1=B1(i)%N
        V1=B1(i)%V
        ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(0,V1,SF)
      else if ( i .GT. NB1+NB2) then
        N1=B3(i-NB1-NB2)%N
        V1=B3(i-NB1-NB2)%V
        M1=B3(i-NB1-NB2)%M
        W1=B3(i-NB1-NB2)%U
        ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N+IBCT(N1,M1),1:3)*FC(0,V1,SC)*FC(0,W1,SA)
      end if
    else if ( j .GT. 0 .and. j .LE. NGS1 ) then
      if ( i .LE. NB1 ) then
        N1=B1(i)%N
        V1=B1(i)%V
        N2=GS1(j)%N
        V2=GS1(j)%V
        IF (N1 .EQ. N2) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(V2,V1,SF)
      else if ( i .GT. NB1 .and. i .LE. NB1+NB2) then
        N1=B2(i-NB1)%N
        V1=B2(i-NB1)%V
        M1=B2(i-NB1)%M
        W1=B2(i-NB1)%U
        N2=GS1(j)%N
        V2=GS1(j)%V
        IF (N1 .EQ. N2) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(V2,V1,SF)*FC(0,W1,0.0d00)
        IF (M1 .EQ. N2) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(0,V1,SF)*FC(V2,W1,0.0d00)
      else if ( i .GT. NB1+NB2) then
        N1=B3(i-NB1-NB2)%N
        V1=B3(i-NB1-NB2)%V
        M1=B3(i-NB1-NB2)%M
        W1=B3(i-NB1-NB2)%U
        N2=GS1(j)%N
        V2=GS1(j)%V
        IF (N1 .EQ. N2 ) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N+IBCT(N1,M1),1:3)*FC(V2,V1,SC)*FC(0,W1,SA)
      end if
    else if ( j .GT. NGS1 ) then 
      if ( i .LE. NB1 ) then
        N1=B1(i)%N
        V1=B1(i)%V
        N2=GS2(j-NGS1)%N
        V2=GS2(j-NGS1)%V
        M2=GS2(j-NGS1)%M
        W2=GS2(j-NGS1)%U
        IF (N1 .EQ. N2 ) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(V2,V1,SF)*FC(W2,0,0.0d0)
        IF (N1 .EQ. M2 ) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(W2,V1,SF)*FC(V2,0,0.0d0)
      else if ( i .GT. NB1 .and. i .LE. NB1+NB2) then
        N1=B2(i-NB1)%N
        V1=B2(i-NB1)%V
        M1=B2(i-NB1)%M
        W1=B2(i-NB1)%U
        N2=GS2(j-NGS1)%N
        V2=GS2(j-NGS1)%V
        M2=GS2(j-NGS1)%M
        W2=GS2(j-NGS1)%U
        IF (N1 .EQ. N2 .AND. M1 .EQ. M2) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(V2,V1,SF)*FC(W2,W1,0.0d0)
        IF (N1 .EQ. M2 .AND. M1 .EQ. N2) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N1,1:3)*FC(W2,V1,SF)*FC(V2,W1,0.0d0)
      else if ( i .GT. NB1+NB2 ) then
        N1=B3(i-NB1-NB2)%N
        V1=B3(i-NB1-NB2)%V
        M1=B3(i-NB1-NB2)%M
        W1=B3(i-NB1-NB2)%U
        N2=GS2(j-NGS1)%N
        V2=GS2(j-NGS1)%V
        M2=GS2(j-NGS1)%M
        W2=GS2(j-NGS1)%U
        IF (N1 .EQ. N2 .AND. M1 .EQ. M2) ETDMD(j,i,:)=ETDMD(j,i,:)+TDME(N+IBCT(N1,M1),1:3)*FC(V2,V1,SC)*FC(W2,W1,SA)
      end if
    end if
  end do
end do
end subroutine create_emisionTDM

subroutine emi_AE(N,maxv,NB,NGS,NGS1,NGS2,GS1,GS2,E,AE,OM0)
implicit none
include "include/types.inc"
integer :: i,j,k
integer :: N,maxv,NB,NGS,NGS1,NGS2
integer :: N1,V1,N2,V2

real*8, dimension(NB) ::E

double precision :: OM0 
double precision, dimension(0:NGS,NB) :: AE

type(P1), dimension(1:N*(maxv+1)) :: GS1
type(P2), dimension(1:N*(maxv+1*(N-1)*maxv)) :: GS2

do i=1,NB
  do j=0,NGS
    if (j .eq. 0 ) then
        AE(j,i)=E(i)
    else if (j .GT. 0 .AND. j .LE. NGS1 ) then
        N1=GS1(j)%N
        V1=GS1(j)%V
       AE(j,i)=E(i)-V1*OM0
    else if (j .GT. NGS1 ) then
        N1=GS2(j-NGS1)%N
        V1=GS2(j-NGS1)%V
        N2=GS2(j-NGS1)%M
        V2=GS2(j-NGS1)%U
       AE(j,i)=E(i)-(V1+V2)*OM0
    end if
  end do
end do

end subroutine emi_AE

subroutine calc_W(NB,E,W,T)
implicit none
integer :: i,j,k
integer :: NB
double precision :: T
real*8, dimension(NB) :: E,W

real*8 :: Q
real*8, parameter :: kb=8.617333E-05

W=0.00
Q=0.00
if (T .GT. 0.00) then
  do i=1,NB
        Q=Q+dexp(-(E(i)-E(1))/(kb*T))
  end do
  do i=1,NB
        W(i)=dexp(-(E(i)-E(1))/(kb*T))/Q
  end do
else
  W(1)=1.00
end if

end subroutine calc_W
