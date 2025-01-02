program Spano
!use types
implicit none
include "include/types.inc"

integer :: i,j,k
integer :: N,Ne,NB1,NB2,NB3,NB4,NB
integer :: maxv,maxvct,Npoints
integer, dimension(:,:), allocatable :: IB1,IBCT
integer, dimension(:,:,:,:), allocatable :: IB2,IB3

type(P1), dimension(:), allocatable :: B1,BCT
type(P2), dimension(:), allocatable :: B2,B3

double precision :: T,sigma,Emin,Emax,Eshift
double precision :: om0,QF,QC,QA,SF,SC,SA,SCA,SFC,SFA
double precision, dimension(:), allocatable :: F,CD_F
double precision, dimension(:,:), allocatable :: TDMe,TDMd,TDMa,R,R0
double precision, dimension(:,:), allocatable :: He, H

double precision, dimension(:,:), allocatable :: TDMm,TDMmd

character(len=200) :: HINP,cutfile
character(len=50) :: outp
character(len=10) :: func,maxS
character(len=2), dimension(3) :: block
character(len=1) :: ani
character(len=1), dimension(3) :: aniso
character(len=10), external :: lowcase

logical :: Absosc,Emiosc,absorp, ECD, wfa,emi,redem,Abstick,Emistick,AbsTemp,ECDvel,CPL
logical :: H_CD

double precision, external :: osc

real*8, dimension(:), allocatable :: work,eigen
integer :: lwork, info

namelist /INP/ block,N,maxv,maxvct,HINP,om0,QF,QC,QA,aniso,cutfile
namelist /PLOT/ emin,emax,sigma,npoints,func,maxS
namelist /SPEC/ absorp, wfa, ECD,emi,T,Eshift,redem,Abstick,Emistick, &
        & Absosc,Emiosc,AbsTemp,ECDvel,CPL

block(:)=''
block(1)='1p'
Emin=-9.0d99
Emax=-8.0d99
sigma=0.05
NPoints=100000
T=0.00
func='gauss'
maxS='H'
Eshift=0.0
absorp=.true.
Absosc=.false.
Emiosc=.false.
AbsTemp=.false.
ECD=.false.
wfa=.false.
emi=.false.
redem=.false.
Abstick=.false.
Emistick=.false.
ECDvel=.false.
CPL=.false.
H_CD=.false.
cutfile=''
aniso=''
aniso(1)='i'

read(*,NML=INP)
read(*,NML=PLOT)
read(*,NML=SPEC)

if ( ECDvel .eqv. .true. ) then
   print*, "WARNING: activated ECD spectra in velocity. ECD=.t."
   print*, "This is incompatible with linear absorption"
   print*, "electric velocity TDM and normal magnetic iTDM (origin dependent) must be used"
   print*, "for ECD in lenght, compatible with linear absorption use electric lenght TDM (usual)"
   print*, "and intrinsic magnetic TDM, estimated as in SI of 10.1021/acs.jctc.0c01100"
end if

allocate(B1(1:N*(maxv+1)))
allocate(IB1(1:N,0:maxv))
allocate(B2(1:N*(maxv+1)*(N-1)*maxv))
allocate(IB2(1:N,0:maxv,1:N,1:maxv))
allocate(B3(1:N*(maxvct+1)*(N-1)*(maxvct+1)))
allocate(IB3(1:N,0:maxvct,1:N,0:maxvct))
allocate(BCT(1:N*(N-1)))
allocatE(IBCT(1:N,1:N))

if (trim(cutfile) == '' ) then
        call create_basis(block,N,Ne,maxv,maxvct,B1,B2,B3,IB1,IB2,IB3,BCT,IBCT,NB,NB1,NB2,NB3)
else
        call cutbasis(block,N,Ne,maxv,maxvct,B1,B2,B3,IB1,IB2,IB3,BCT,IBCT,NB,NB1,NB2,NB3,cutfile)
end if
print*, "Create basis... OK"
call flush(6)

allocate(He(Ne,Ne))
allocate(R0(Ne,3))
allocate(H(1:NB,1:NB))
allocate(TDMe(NE,3))
allocate(TDMm(NE,3))
allocate(TDMd(1:NB,3))
allocate(TDMmd(1:NB,3))
allocate(R(1:NB,3))
allocate(TDMa(1:NB,3))
allocate(F(1:NB))
allocate(CD_F(1:NB))
allocate(eigen(NB))
He(:,:)=0.00d00
H(:,:)=0.00d00

!!!!
TDMe(:,:)=0.00d00
TDMd(:,:)=0.00d00
TDMa(:,:)=0.00d00
R(:,:)= 0.0d00
R0(:,:)= 0.0d00
!TDMe(1:N,1)=1.00d00
!!!!

if ( ECD .eqv. .true. ) H_CD=.true.
if ( CPL .eqv. .true. ) H_CD=.true.

call read_He(Ne,He,HINP,TDMe,TDMm,block,R0,H_CD,cutfile)
call calculate_HR(QF,QC,QA,SF,SC,SA,SCA,SFC,SFA)

call init_TDM(TDMd,TDMe,N,Ne,maxv,maxvct,NB1,NB2,NB3,NB,B1,B2,B3,BCT,IB1,IB2,IB3,IBCT,SF,SC,SA)
call init_H(N,Ne,maxv,maxvct,NB,NB1,NB2,NB3,IB1,IB2,IB3,IBCT,B1,B2,B3,om0,SF,SC,SA,SCA,SFC,SFA,He,H)
print*, "Generated Hamiltonian... OK"
call flush(6)

!do i=1,NB
!do j=1,i
!        print*, i,j,H(i,j)
!end do
!end do
!stop
if (trim(aniso(1)) /= 'i') then 
        call polarized(NB,TDMd,aniso)
        print*, "Polarized intensity... OK"
        call flush(6)
end if

allocate(work(1))
info=0
lwork=-1
call DSYEV('V','L',NB,H,NB,eigen,work,lwork,info)
if (info /= 0) then
   stop 'dsyev 1 failed'
endif
lwork = int(work(1))
deallocate(work)

allocate(work(lwork))
call DSYEV('V','U',NB,H,NB,eigen,work,lwork,info)
if (info /= 0) then
   stop 'dsyev 2 failed'
endif
deallocate(work)

print*, "Diagonalize Hamiltonian... OK"
call flush(6)
!call DGEMV('T',size(H,1),size(H,2),1.0d00,H,size(H,1),TDMd,1,0.0d00,TDMa,1)

if (absorp .eqv. .true. ) then
outp='abs.dat'
if (AbsTemp .eqv. .true. ) then
  call spc(N,NB1,NB2,NB3,NB,B1,B2,B3,IBCT,BCT,maxv,maxvct,TDMe,OM0,SF,SC,SA,&
  & eigen,H,T,Npoints,Emin,Emax,sigma,func,maxS,outp,Eshift,.false.,.false.,Absosc,Abstick)

else
  TDMa(:,:)=0.0d00
  TDMa=matmul(transpose(H),TDMd)
  F=0.0
  if (AbsOsc .eqv. .true. ) then
    do i=1,NB
      F(i)=osc(Eshift+eigen(i),TDMa(i,:))
    end do
  else
    do i=1,NB
      F(i)=dot_product(TDMa(i,:),TDMa(i,:))
    end do
  end if
  
  print*, "Get Absorption intensities... OK"
  call flush(6)

  call conv(eigen,F,Emin,Emax,sigma,Npoints,NB,func,maxS,outp)
  print*, "Absorption Spectra... OK"
  call flush(6)

  if ( abstick .eqv. .true. ) then
    open(12,file='abs_stick.dat')
    do i=1,NB
      write(12,*) eigen(i),F(i)
    end do
    close(12)
  end if
end if
end if

if ( wfa .eqv. .true. ) then
        call wvf_analisis(NB,NB1,NB2,NB3,H,eigen)
        print*, "Wfa Analisys... OK"
        call flush(6)
end if

if ( ECD .eqv. .true. ) then
        outp='CD.dat'
        call init_TDM(TDMmd,TDMm,N,Ne,maxv,maxvct,NB1,NB2,NB3,NB,B1,B2,B3,BCT,IB1,IB2,IB3,IBCT,SF,SC,SA)
        call create_R(R,R0,N,Ne,NB,NB1,NB2,NB3,B1,B2,B3,maxv,maxvct)
        call CD_abs(NB,TDMd,TDMmd,H,eigen,CD_F,R,ECDvel) 
        print*, "ECD Intensity... OK"
        call flush(6)

        print*, CD_F
        call conv(eigen,CD_F,Emin,Emax,sigma,Npoints,NB,func,maxS,outp)
        print*, "ECD Spectra... OK"
        call flush(6)
end if

if ( emi .eqv. .true. ) then
        outp='emi.dat'
        call spc(N,Ne,NB1,NB2,NB3,NB,B1,B2,B3,IBCT,BCT,maxv,maxvct,TDMe,TDMm,R,OM0,SF,SC,SA,&
                & eigen,H,T,Npoints,Emin,Emax,sigma,func,maxS,outp,Eshift,.true.,&
                & ECDvel,.false.,redem,Emiosc,Emistick)
end if

if ( CPL .eqv. .true. ) then
        outp='CPL.dat'
        call create_R(R,R0,N,Ne,NB,NB1,NB2,NB3,B1,B2,B3,maxv,maxvct)
        call spc(N,Ne,NB1,NB2,NB3,NB,B1,B2,B3,IBCT,BCT,maxv,maxvct,TDMe,TDMm,R,OM0,SF,SC,SA,&
                & eigen,H,T,Npoints,Emin,Emax,sigma,func,maxS,outp,Eshift,.true.,&
                & ECDvel,.true.,redem,Emiosc,Emistick)
end if

end program 
