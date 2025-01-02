subroutine wvf_analisis(NB,NB1,NB2,NB3,U,E)
implicit none
integer :: i,j,k
integer :: NB,NB1,NB2,NB3

double precision, dimension(NB) :: E
double precision, dimension(NB,NB) :: U
double precision, dimension(:,:), allocatable :: Nat

open(10,file='wvf.dat')
allocate(Nat(Nb,3))
Nat(:,:)=0.00d00
do i=1,NB
        do j=1,NB1
                Nat(i,1)=Nat(i,1)+U(j,i)**2
        end do
        do j=1,NB2
                Nat(i,2)=Nat(i,2)+U(NB1+j,i)**2
        end do
        do j=1,NB3
                Nat(i,3)=Nat(i,3)+U(NB1+NB2+j,i)**2
        end do
        write(10,*) i,E(i),Nat(i,1:3)
end do
close(10)
deallocate(Nat)
end subroutine


