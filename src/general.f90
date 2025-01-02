
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CALCULATE THE FACTORIAL OF AN INTEGER %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

function FAC(N) 
implicit none
integer :: N,II

double precision :: FAC

FAC=1.0D0
do II=1,N
  FAC=FAC*II
end do
end function FAC

pure Function Lowcase(str) Result (string)

!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(len(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(cap, str(i:i))
        if (ic > 0) string(i:i) = low(ic:ic)
    end do

End Function Lowcase  

function gauss(x,x0,I,a,info) result(y)
implicit none
real*8 :: x,x0,y,I,a
real*8, parameter :: pi=3.14159265359
character(len=1) :: info 
if ( a < 0.0 ) then
	write(*,*) "ERROR: Sigma of gaussian fucntion must be positive"
	stop
end if

if ( info == 'I' ) then
y=I*(1.0/(a*dsqrt(5.0d00*pi)))*dexp(-((x-x0)**2)/(2*a*a))
else
y=I*dexp(-((x-x0)**2)/(2*a*a))
end if
return
end function gauss

function lorentz(x,x0,I,a,info) result(y)
implicit none
real*8 :: x,x0,y,I,a
real*8, parameter :: pi=3.14159265359
character(len=1) :: info
if ( a < 0.0 ) then
	write(*,*) "ERROR: Sigma of gaussian fucntion must be positive"
	stop
end if

if ( info == 'I' ) then
y=I*(1.0/pi)*(0.5*a)/(((x-x0)**2)+((0.5*a)**2))
else
print*, "NOT IMPLEMENTED INT=H for Lorentz"
stop
end if
return
end function lorentz

function osc(E,U)
implicit none
double precision, parameter :: ev2au=1.0/27.2116
double precision :: E,osc
double precision, dimension(3) :: U

osc=2.0/3.0*dot_product(U,U)*abs(E)*ev2au

end function osc

subroutine cross(a, b,c) 
double precision, DIMENSION(3), INTENT(OUT) :: c
double precision, DIMENSION(3), INTENT(IN) :: a, b

c(1) = a(2) * b(3) - a(3) * b(2)
c(2) = a(3) * b(1) - a(1) * b(3)
c(3) = a(1) * b(2) - a(2) * b(1)
END subroutine cross
