subroutine create_basis(bb,N,Ne,maxv,maxvct,B1,B2,B3,IB1,IB2,IB3,BCT,IBCT,NB,NB1,NB2,NB3)
implicit none
include "include/types.inc"
integer :: i,j,k,l
integer :: N,Ne,NB1,NB2,NB3,NB
integer :: maxv,maxvct
integer :: p,pp,ct,tb

integer, dimension(1:N,1:N) :: IBCT
integer, dimension(1:N,0:maxv) :: IB1
integer, dimension(1:N,0:maxv,1:N,1:maxv) :: IB2
integer, dimension(1:N,0:maxvct,1:N,0:maxvct) :: IB3

type(P1), dimension(1:N*(N-1)) :: BCT
type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2
type(P2), dimension(1:N*(maxvct+1)*(N-1)*(maxvct+1)) :: B3

character(len=2), dimension(3) :: bb
character(len=2), external :: lowcase

p=0
pp=0
ct=0
tb=0

do i=1,3
        if ( lowcase(bb(i)) == '1p' ) then
                p=1
        else if ( lowcase(bb(i)) == '2p' ) then
                pp=1
        else if ( lowcase(bb(i)) == 'ct' ) then
                ct=1
        end if
end do

if ( ct == 1 ) then
   Ne=N*N
else
   Ne=N
end if

!Electronic indexes for CT states
if ( ct .eq. 1 ) then
k=0
do i=1,N
  do j=1,N
    if (i .ne. j) then
      k=k+1
      BCT(k)%N=i  
      BCT(k)%V=j 
      IBCT(i,j)=k
    end if
  end do
end do
end if


!ONE-PARTICLE BASIS SET 
NB1=0
if ( p .eq. 1 ) then
do i=1,N
  do j=0,maxv
    NB1=NB1+1
    B1(NB1)%N=i
    B1(NB1)%V=j
    IB1(i,j)=NB1
  end do
end do
end if

!TWO-PARTICLE BASIS SET
NB2=0
if ( pp .eq. 1 ) then
do i=1,N
  do j=0,maxv
    do k=1,N
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

!CT BASIS
NB3=0
if ( ct .eq. 1 ) then
do i=1,N
  do j=0,maxvct
    do k=1,N
      if (i.ne.k) then
        do l=0,maxvct-j
          NB3=NB3+1
          B3(NB3)%N=i
          B3(NB3)%V=j
          B3(NB3)%M=k
          B3(NB3)%U=l
          IB3(i,j,k,l)=NB2+NB1+NB3
        end do
     end if 
    end do
  end do
end do
else if ( tb .eq. 1 ) then
do i=1,N
  do j=0,maxvct
    do k=1,N
      if (i.ne.k) then
      if (abs(i-k) .eq. 1 .or. abs(i-k) .eq. N-1 ) then
        do l=0,maxvct-j 
          NB3=NB3+1
          B3(NB3)%N=i
          B3(NB3)%V=j
          B3(NB3)%M=k
          B3(NB3)%U=l
          IB3(i,j,k,l)=NB2+NB1+NB3
        end do
     end if 
     end if 
    end do
  end do
end do
end if

NB=NB1+NB2+NB3
end subroutine create_basis 

subroutine cutbasis(bb,N,Ne,maxv,maxvct,B1,B2,B3,IB1,IB2,IB3,BCT,IBCT,NB,NB1,NB2,NB3,cutfile)
!use types
implicit none
include "include/types.inc"
integer :: i,j,k,l
integer :: N,Ne,NB1,NB2,NB3,NB
integer :: maxv,maxvct
integer :: p,pp,ct,tb
integer :: NF,NCT

integer, dimension(1:N,1:N) :: IBCT
integer, dimension(1:N,0:maxv) :: IB1
integer, dimension(1:N,0:maxv,1:N,1:maxv) :: IB2
integer, dimension(1:N,0:maxvct,1:N,0:maxvct) :: IB3

type(P1), dimension(1:N*(N-1)) :: BCT
type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2
type(P2), dimension(1:N*(maxvct+1)*(N-1)*(maxvct+1)) :: B3

character(len=200) :: cutfile
character(len=2), dimension(3) :: bb
character(len=2), external :: lowcase

p=0
pp=0
ct=0
tb=0

do i=1,3
        if ( lowcase(bb(i)) == '1p' ) then
                p=1
        else if ( lowcase(bb(i)) == '2p' ) then
                pp=1
        else if ( lowcase(bb(i)) == 'ct' ) then
                ct=1
        else if ( lowcase(bb(i)) == 'tb' ) then
                tb=1
                ct=0
        end if
end do

open(10,file=trim(cutfile))
read(10,*) NF,NCT
do i=1,NF
        read(10,*)
end do
do i=1,NCT
      read(10,*) BCT(i)%N ,BCT(i)%V
      IBCT(BCT(i)%N,BCT(i)%V)=i
      !print*, i,BCT(i)%N ,BCT(i)%V
end do
Ne=NF+NCT

!ONE-PARTICLE BASIS SET 
NB1=0
if ( p .eq. 1 ) then
do i=1,N
  do j=0,maxv
    NB1=NB1+1
    B1(NB1)%N=i
    B1(NB1)%V=j
    IB1(i,j)=NB1
  end do
end do
end if

!TWO-PARTICLE BASIS SET
NB2=0
if ( pp .eq. 1 ) then
do i=1,N
  do j=0,maxv
    do k=1,N
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

!CT BASIS
NB3=0
if ( ct .eq. 1 ) then
do i=1,NCT
  do j=0,maxvct
        do l=0,maxvct-j
          NB3=NB3+1
          B3(NB3)%N=BCT(i)%N
          B3(NB3)%V=j
          B3(NB3)%M=BCT(i)%V
          B3(NB3)%U=l
          !print*, i,B3(NB3)
          IB3(BCT(i)%N,j,BCT(i)%V,l)=NB2+NB1+NB3
        end do
  end do
end do
end if

NB=NB1+NB2+NB3
end subroutine cutbasis 


