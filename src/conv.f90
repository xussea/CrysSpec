subroutine conv(E,F,Emin,Emax,a,NP,NS,ftype,Intinfo,outp)
implicit none
integer :: i,j,k
integer :: NP,NS

double precision :: a,Emin,Emax,E0,dx
double precision, dimension(NS) :: E
double precision, dimension(:), allocatable :: x,y
double precision, dimension(NS) :: F

double precision, external :: gauss, lorentz

character(len=50) :: outp
character(len=10) :: ftype,Intinfo
character(len=10), external :: lowcase

open(10,file=trim(outp))

if ( Emin == -9.0d99 .and. Emax == -8.0d99 ) then
  E0=minval(E)-2.5*a
  dx=(maxval(E)+2.5*a-E0)/(NP-1)
else if ( Emin /= 9.0d99 .and. Emax == 8.0d99 ) then
  E0=Emin
  dx=(maxval(E)+2.5*a-E0)/(NP-1)
else if ( Emin == 9.0d99 .and. Emax /= 8.0d99 ) then
  E0=minval(E)-2.5*a
  dx=(Emax-E0)/(NP-1)
else if ( Emin /= 9.0d99 .and. Emax /= 8.0d99 ) then
  if (Emax <= Emin ) then
    print*, "ERROR limits of x axis are not correct"
    print*, "from", Emin, "to",Emax
    stop
  else
    E0=Emin
    dx=(Emax-E0)/(NP-1)
  end if
end if

allocate(x(NP))
allocate(y(NP))
x(:)=0.0
y(:)=0.0


!CONVOLUTION
!$OMP PARALLEL DEFAULT(SHARED) private(i,j)
if (trim(lowcase(ftype)) == 'gauss') then
!$OMP DO
  do i=1,NP
    x(i)=E0+i*dx
    do j=1,NS
      y(i)=y(i)+gauss(x(i),E(j),F(j),a,Intinfo)
    end do
  end do
!$OMP END DO
else if (trim(lowcase(ftype)) == 'lorentz') then
!$OMP DO
  do i=1,NP
    x(i)=E0+i*dx
    do j=1,NS
      y(i)=y(i)+lorentz(x(i),E(j),F(j),a,Intinfo)
    end do
  end do
!$OMP END DO
end if
!$OMP BARRIER
!$OMP END PARALLEL
!y=y/maxval(y)

!Print output
do i=1,NP
  write(10,'(2(F15.5,3X))') x(i),y(i)
end do
close(10)
end subroutine conv


