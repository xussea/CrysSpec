subroutine spc(N,Ne,NB1,NB2,NB3,NB,B1,B2,B3,IBCT,BCT,maxv,maxvct,TDMe,OM0,&
&          SF,SC,SA,E,C,T,NP,Emin,Emax,sigma,func,maxS,outp,shift,emi,redem,Iosc,stick)
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
double precision, dimension(:,:), allocatable :: AE,WETDM
double precision, dimension(:,:,:), allocatable :: ETDMD,ETDMA
double precision, dimension(:), allocatable :: EE,F

double precision, parameter :: ev2au=1.0d0/27.2116d0

double precision, external :: osc
character(len=10) :: func,maxS
character(len=50) :: outp

logical :: twopemi, emi,redem,Iosc,stick

allocate(GS1(1:N*(maxv+1)))
allocate(GS2(1:N*(maxv+1)*(N-1)*maxv))
allocate(IGS1(1:N,0:maxv))
allocate(IGS2(1:N,0:maxv,1:N,1:maxv))

twopemi=.false.
if (NB2 .gt. 0 ) twopemi=.true.

call basisGS(N,maxv,maxvct,GS1,GS2,IGS1,IGS2,NGS,NGS1,NGS2,twopemi)
print*, "Ground state basis... OK"
call flush(6)

allocate(ETDMD(0:NGS,1:NB,3))
allocate(ETDMA(0:NGS,1:NB,3))
allocate(WETDM(0:NGS,1:NB))
allocate(AE(0:NGS,1:NB))
allocate(EE(NB*(NGS+1)))
allocate(F(NB*(NGS+1)))

call create_emisionTDM(N,Ne,maxv,maxvct,NB1,NB2,NB3,NB,NGS,NGS1,NGS2,B1,B2,B3,GS1,GS2,IBCT,BCT,ETDMD,TDMe,SF,SC,SA)
print*, "Extended ground state intensities... OK"
call flush(6)

ETDMA=0.00d00
WETDM=0.00d00
EE=0.0d00
F=0.0d00

do i=1,3
        call DGEMM('N','N',NGS+1,NB,NB,1.0d00,ETDMD(:,:,i),NGS+1,C,NB,0.0d00,ETDMA(:,:,i),NGS+1)
        !ETDMA(:,:,i)=matmul(ETDMD(:,:,i),C)
end do

call emi_AE(N,maxv,NB,NGS,NGS1,NGS2,GS1,GS2,E,AE,OM0,emi)
print*, "Calculate extended transition energies... OK"
call flush(6)


if (emi .eqv. .true. ) then
 
   do i=1,NB
     do j=0,NGS
       WETDM(j,i)=dot_product(ETDMA(j,i,:),ETDMA(j,i,:)) 
     end do
   end do
   print*, "Emission intensity... OK"
   call flush(6)

   if (redem .eqv. .true. ) then
   !!!!Cubic term of reduced emission
   do i=1,NGS1
     do j=1,NB
        V1=GS1(i)%V
        WETDM(i,j)=WETDM(i,j)*abs((1-V1*OM0/(AE(0,j)+shift))**3)
     end do
   end do

   do i=1,NGS2
     do j=1,NB
        V1=GS2(i)%V
        V2=GS2(i)%U
        WETDM(NGS1+i,j)=WETDM(NGS1+i,j)*abs((1-(V1+V2)*OM0/(AE(0,j)+shift))**3)
     end do
   end do

    print*, "Reduced emission intensity... OK"
    call flush(6)
  end if

  if (Iosc .eqv. .true. ) then
    do i=1,NB
      WETDM(:,i)=WETDM(:,i)*((AE(0,j)+shift)*ev2au)**3
    end do
  end if

  !!!!Thermal average
  allocate(W(NB))
  call calc_W(NB,E,W,T)
  print*, "Calculated Thermal weights... OK"
  call flush(6)

  k=0
  do i=1,NB
    do j=0,NGS
      k=k+1
      EE(k)=AE(j,i)
      F(k)=WETDM(j,i)*W(i)
      !print*, k, EE(k),F(k),W(i)
    end do
  end do
  print*, "Emission thermalized... OK"
  call flush(6)
 
  if (stick .eqv. .true. ) then
    open(12,file='emi_stick_temp.dat')
    k=0
    do i=1,NB
      do j=0,NGS
        write(12,*) EE(k), F(k)
      end do
    end do
    close(12)
  end if  

else 

   allocate(W(NGS+1))
   allocate(EGS(NGS+1))
   do i=0,NGS
        EGS(i+1)=E(1)-AE(i,1)
   end do
   call calc_W(NGS+1,EGS,W,T)
   deallocate(EGS)
  print*, "Calculated Thermal weights... OK"
  call flush(6)


   k=0
   do i=1,NB
     do j=0,NGS
       k=k+1
       EE(k)=AE(j,i)
       if (Iosc .eqv. .true. ) then
         F(k)=osc(E(i),ETDMA(j,i,:))*W(j+1)
       else
         F(k)=dot_product(ETDMA(j,i,:),ETDMA(j,i,:))*W(j+1)
       end if
     end do
   end do
   Print*, "Thermalized Absorption... OK"
   Call flush(6)

  if (stick .eqv. .true. ) then
    open(12,file='abs_stick_temp.dat')
    k=0
    do i=1,NB
      do j=0,NGS
        write(12,*) EE(k), F(k)
      end do
    end do
    close(12)
  end if  


end if


call conv(EE,F,Emin,Emax,sigma,NP,NB*(NGS+1),func,maxS,outp)

return
end subroutine spc

