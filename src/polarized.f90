subroutine polarized(N,U,aniso)
implicit none
integer :: i,j,k
integer :: N
double precision, dimension(N,3) :: U
double precision, dimension(:,:), allocatable :: D

character(len=1), dimension(3) :: aniso
logical :: x,y,z
character(len=1), external :: lowcase

allocate(D(size(U,1),size(U,2)))
D=U
U(:,:)=0.00
x=.false.
y=.false.
z=.false.
do i=1,3
        if (lowcase(trim(aniso(i))) == 'x' ) then
                x=.true.
        else if (lowcase(trim(aniso(i))) == 'y' ) then
                y=.true.
        else if (lowcase(trim(aniso(i))) == 'z' ) then
                z=.true.
        end if
end do

if (x .eqv. .true. ) then
        U(:,1)=D(:,1)
end if
if (y .eqv. .true. ) then
        U(:,2)=D(:,2)
end if
if (z .eqv. .true. ) then
        U(:,3)=D(:,3)
end if
deallocate(D)

end subroutine polarized


