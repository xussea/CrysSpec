subroutine spc(N,Ne,NB1,NB2,NB3,NB,B1,B2,B3,IBCT,BCT,maxv,maxvct,TDMe,OM0,&
&          SF,SC,SA,E,C,T,NP,Emin,Emax,sigma,func,maxS,outp,shift,emi,redem)
implicit none
include "include/types.inc"

integer :: i,j,k
integer :: V1,V2
integer :: N,Ne,NB1,NB2,NB3,NB,NGS,NGS1,NGS2,NGS3
integer :: maxv,maxvct
integer :: NP

type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2
type(P2), dimension(1:N*(maxvct+1)*(N-1)*(maxvct+1)) :: B3

type(P1), dimension(:), allocatable :: GS1
type(P2), dimension(:), allocatable :: GS2
integer, dimension(:,:), allocatable :: IGS1
integer, dimension(:,:,:,:), allocatable :: IGS2

integer, dimension(1:N,1:N) :: IBCT
type(P1), dimension(1:N*(N-1)) :: BCT

real*8, dimension(NB) :: E
real*8, dimension(:), allocatable :: W,EGS

double precision :: sigma, Emin, Emax
double precision :: SF,SA,SC,OM0,T,shift,ss
double precision, dimension(NB,NB) :: C
double precision, dimension(Ne,3) :: TDMe
double precision, dimension(:,:), allocatable :: AE
double precision, dimension(:,:,:), allocatable :: ETDMD,ETDMA,WETDM
double precision, dimension(:), allocatable :: EE,F

character(len=10) :: func,maxS
character(len=50) :: outp

logical :: emi,redem

allocate(GS1(1:N*(maxv+1)))
allocate(GS2(1:N*(maxv+1)*(N-1)*maxv))
allocate(IGS1(1:N,0:maxv))
allocate(IGS2(1:N,0:maxv,1:N,1:maxv))

call basisGS(N,maxv,maxvct,GS1,GS2,IGS1,IGS2,NGS,NGS1,NGS2)

allocate(ETDMD(0:NGS,1:NB,3))
allocate(ETDMA(0:NGS,1:NB,3))
allocate(WETDM(0:NGS,1:NB,3))
allocate(AE(0:NGS,1:NB))

call create_emisionTDM(N,Ne,maxv,maxvct,NB1,NB2,NB3,NB,NGS,NGS1,NGS2,B1,B2,B3,GS1,GS2,IBCT,BCT,ETDMD,TDMe,SF,SC,SA)
ETDMA=0.00d00

do i=1,3
        call DGEMM('N','N',NGS+1,NB,NB,1.0d00,ETDMD(:,:,i),NGS+1,C,NB,0.0d00,ETDMA(:,:,i),NGS+1)
        !ETDMA(:,:,i)=matmul(ETDMD(:,:,i),C)
end do

call emi_AE(N,maxv,NB,NGS,NGS1,NGS2,GS1,GS2,E,AE,OM0,emi)

if (emi .eqv. .true. ) then
   
   if (redem .eqv. .true. ) then
   !!!!Cubic term of reduced emission
   do i=1,NGS1
     do j=1,NB
        V1=GS1(i)%V
        ETDMA(i,j,:)=ETDMA(i,j,:)*sqrt(abs(((AE(0,j)+shift)-V1*OM0)**3))
     end do
   end do

   do i=1,NGS2
     do j=1,NB
        V1=GS2(i)%V
        V2=GS2(i)%U
        ETDMA(NGS1+i,j,:)=ETDMA(NGS1+i,j,:)*sqrt(abs(((AE(0,j)+shift)-(V1+V2)*OM0)**3))
     end do
   end do
  end if
  !!!!Thermal average
  allocate(W(NB))
  call calc_W(NB,E,W,T)
  do i=1,NB
          WETDM(:,i,:)=ETDMA(:,i,:)*sqrt(W(i))
  end do
  
else 

   allocate(W(NGS+1))
   allocate(EGS(NGS+1))
   do i=0,NGS
        EGS(i+1)=E(1)-AE(i,1)
   end do
   call calc_W(NGS+1,EGS,W,T)
   do i=0,NGS
        WETDM(i,:,:)=ETDMA(i,:,:)*sqrt(W(i+1))
   end do
   deallocate(EGS)

end if

allocate(EE(NB*(NGS+1)))
allocate(F(NB*(NGS+1)))

k=0
do i=1,NB
  do j=0,NGS
    k=k+1
    EE(k)=AE(j,i)
    F(k)=dot_product(WETDM(j,i,:),WETDM(j,i,:))
    !print*, k, EE(k),F(k),W(i)
  end do
end do

call conv(EE,F,Emin,Emax,sigma,NP,NB*(NGS+1),func,maxS,outp)

deallocate(GS1)
deallocate(GS2)
deallocate(IGS1)
deallocate(IGS2)
deallocate(ETDMA)
deallocate(ETDMD)
deallocate(WETDM)
deallocate(AE)
deallocate(W)

end subroutine spc

