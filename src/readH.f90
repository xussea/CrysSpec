subroutine read_He(N,He,input,TDM,TDMm,block,R,ECD,cutfile)
implicit none
integer :: i,j,k
integer :: N,NN,NF,NCT

double precision, dimension(N,N) :: He
double precision, dimension(N,3) :: TDM,R

double precision, dimension(N,3) :: TDMm

character(len=200) :: input,cutfile
character(len=2), dimension(3) :: block

character(len=2), external :: lowcase

logical :: ECD

open(10,file=trim(input))

He(:,:)=0.00d00
do i=1,N
  read(10,*) He(i,1:N)
end do

!do j=1,NN
!  do i=1,i
!    H(i,j)=H(j,i)
!  end do  
!end do

TDM(:,:)=0.00d00
do i=1,N
  read(10,*) TDM(i,1:3)
end do

if (ECD .eqv. .true. ) then
TDMm=0.0d0
read(10,*)
do i=1,N
        read(10,*) TDMm(i,:)
end do
read(10,*) 

R(:,:)=0.00d00
do i=1,N
  read(10,*) R(i,1:3)
end do
end if

close(10)

end subroutine read_He


