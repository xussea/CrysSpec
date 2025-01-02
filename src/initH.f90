subroutine init_H(N,Ne,maxv,maxvct,NB,NB1,NB2,NB3,IB1,IB2,IB3,IBCT,B1,B2,B3,om0,SF,SC,SA,SCA,SFC,SFA,JD,H)
!use types
implicit none
include "include/types.inc"
integer :: II,JJ
integer :: i,j,k,l
integer :: N,Ne,NB1,NB2,NB3,NB
integer :: N1,N2,M1,M2,V1,V2,U1,U2 
integer :: maxv,maxvct

double precision :: flagh, om0, SF,SC,SA,SCA,SFC,SFA
double precision, dimension(NB,NB) :: H
double precision, dimension(Ne,Ne) :: JD

integer, dimension(1:N,1:N) :: IBCT
integer, dimension(1:N,0:maxv) :: IB1
integer, dimension(1:N,0:maxv,1:N,1:maxv) :: IB2
integer, dimension(1:N,0:maxvct,1:N,0:maxvct) :: IB3

type(P1), dimension(1:N*(N-1)) :: BCT
type(P1), dimension(1:N*(maxv+1)) :: B1
type(P2), dimension(1:N*(maxv+1)*(N-1)*maxv) :: B2
type(P2), dimension(1:N*(maxvct+1)*(N-1)*(maxvct+1)) :: B3

double precision, external :: FC

do II=1,NB
  do JJ=1,II
    FLAGH=0.0D0

    !CALCULATE <N1*,V1|H|N2*,V2>
    IF (II .LE. NB1) THEN
    IF (JJ .LE. NB1) THEN
      N1=B1(II)%N
      V1=B1(II)%V
      N2=B1(JJ)%N
      V2=B1(JJ)%V
      IF(N1.EQ.N2 .AND. V1.EQ.V2) FLAGH=FLAGH+OM0*V2+JD(N2,N2)
      IF(N1.NE.N2) FLAGH=FLAGH+JD(N1,N2)*FC(0,V1,SF)*FC(0,V2,SF)
    END IF
    END IF

    !CALCULATE <N1*,V1;M1,U1|H|N2*,V2>
    IF (II .GT. NB1 .AND. II .LE. NB1+NB2) THEN
    IF (JJ .LE. NB1) THEN
      N1=B2(II-NB1)%N
      V1=B2(II-NB1)%V
      M1=B2(II-NB1)%M
      U1=B2(II-NB1)%U
      N2=B1(JJ)%N
      V2=B1(JJ)%V
      IF(M1.EQ.N2) FLAGH=FLAGH+JD(N1,N2)*FC(U1,V2,SF)*FC(0,V1,SF)
    END IF
    END IF

    !CALCULATE <N1*,V1;M1,U1|H|N2*,V2;M2,U2>
    IF (II .GT. NB1 .AND. II .LE. NB1+NB2) THEN
    IF (JJ .GT. NB1 .AND. JJ .LE. NB1+NB2) THEN
      N1=B2(II-NB1)%N
      V1=B2(II-NB1)%V
      M1=B2(II-NB1)%M
      U1=B2(II-NB1)%U
      N2=B2(JJ-NB1)%N
      V2=B2(JJ-NB1)%V
      M2=B2(JJ-NB1)%M
      U2=B2(JJ-NB1)%U
      IF(N1.EQ.N2 .AND. V1.EQ.V2 .AND. M1.EQ.M2 .AND. U1.EQ.U2) FLAGH=FLAGH+OM0*(V2+U2)+JD(N2,N2)
      IF(M1.EQ.N2 .AND. N1.EQ.M2) FLAGH=FLAGH+JD(N1,N2)*FC(U1,V2,SF)*FC(U2,V1,SF)
      IF(M1.EQ.M2 .AND. N1.NE.N2 .AND. U1.EQ.U2) FLAGH=FLAGH+JD(N1,N2)*FC(0,V1,SF)*FC(0,V2,SF)
    END IF
    END IF
    
    !CALCULATE <N1+,V1;M1-,U1|H|N2+,V2;M2-,U2>
    IF (II .GT. NB1+NB2 .AND. II .LE. NB1+NB2+NB3) THEN
    IF (JJ .GT. NB1+NB2 .AND. JJ .LE. NB1+NB2+NB3) THEN
      N1=B3(II-NB1-NB2)%N
      V1=B3(II-NB1-NB2)%V
      M1=B3(II-NB1-NB2)%M
      U1=B3(II-NB1-NB2)%U
      N2=B3(JJ-NB1-NB2)%N
      V2=B3(JJ-NB1-NB2)%V
      M2=B3(JJ-NB1-NB2)%M
      U2=B3(JJ-NB1-NB2)%U
      IF(N1.EQ.N2 .AND. V1.EQ.V2 .AND. M1.EQ.M2 .AND. U1.EQ.U2) FLAGH=FLAGH+OM0*(V2+U2)+JD(N+IBCT(N1,M1),N+IBCT(N2,M2)) 
      IF(M1.EQ.N2 .AND. N1.EQ.M2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N+IBCT(N2,M2))*FC(U1,V2,SCA)*FC(U2,V1,SCA)
      IF(M1.EQ.N2 .AND. N1.NE.M2 .AND. U1.EQ.V2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N+IBCT(N2,M2))*FC(U1,V2,SCA)*FC(0,V1,SC)*FC(0,U2,SA)
      IF(M2.EQ.N1 .AND. N2.NE.M1 .AND. V1.EQ.U2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N+IBCT(N2,M2))*FC(V1,U2,SCA)*FC(0,U1,SC)*FC(0,V2,SA)
      IF(M1.EQ.M2 .AND. N1.NE.N2 .AND. U1.EQ.U2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N+IBCT(N2,M2))*FC(0,V1,SC)*FC(0,V2,SC)
      IF(N1.EQ.N2 .AND. M1.NE.M2 .AND. V1.EQ.V2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N+IBCT(N2,M2))*FC(0,U1,SA)*FC(0,U2,SA)
      IF (N1 .NE. N2 .AND. N1 .NE. M1 .AND. N1 .NE. M2 ) then
      IF (N2 .NE. M1 .AND. N2 .NE. M2 ) then
              IF (M1 .NE. M2 )FLAGH=FLAGH+JD(N+IBCT(N1,M1),N+IBCT(N2,M2))*FC(0,U1,SA)*FC(0,U2,SA)*FC(0,V1,SC)*FC(0,V2,SC)
      END IF
      END IF
    END IF
    END IF

    !CALCULATE <N1+,V1;M1-,U1|H|N2*,V2>
    IF (II .GT. NB1+NB2 .AND. II .LE. NB1+NB2+NB3) THEN
    IF (JJ .LE. NB1) THEN
      N1=B3(II-NB1-NB2)%N
      V1=B3(II-NB1-NB2)%V
      M1=B3(II-NB1-NB2)%M
      U1=B3(II-NB1-NB2)%U
      N2=B1(JJ)%N
      V2=B1(JJ)%V
      IF(M1.EQ.N2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N2)*FC(U1,V2,SFA)*FC(0,V1,SC)
      IF(N1.EQ.N2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N2)*FC(V1,V2,SFC)*FC(0,U1,SA)
    END IF
    END IF

    !CALCULATE <N1+,V1;M1-,U1|H|N2*,V2;M2,U2>
    IF (II .GT. NB1+NB2 .AND. II .LE. NB1+NB2+NB3) THEN
    IF (JJ .GT. NB1 .AND. JJ .LE. NB1+NB2) THEN
      N1=B3(II-NB1-NB2)%N
      V1=B3(II-NB1-NB2)%V
      M1=B3(II-NB1-NB2)%M
      U1=B3(II-NB1-NB2)%U
      N2=B2(JJ-NB1)%N
      V2=B2(JJ-NB1)%V
      M2=B2(JJ-NB1)%M
      U2=B2(JJ-NB1)%U
      IF(N1.EQ.N2 .AND. M1.EQ.M2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N2)*FC(V1,V2,SFC)*FC(U2,U1,SA)
      IF(N1.EQ.M2 .AND. N2.EQ.M1) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N2)*FC(U2,V1,SC)*FC(U1,V2,SFA)
      IF(N1.NE.N2 .AND. M1.EQ.M2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N2)*FC(0,V1,SC)*FC(0,V2,SF)*FC(U2,U1,SA)
      IF(M1.NE.N2 .AND. N1.EQ.M2) FLAGH=FLAGH+JD(N+IBCT(N1,M1),N2)*FC(0,U1,SA)*FC(0,V2,SF)*FC(U2,V1,SC)
      IF(N1.NE.N2 .AND. N1.NE.N2 .AND. N2.NE.M1 .AND. M1.NE.M2) THEN 
        FLAGH=FLAGH+JD(N+IBCT(N1,M1),N2)*FC(0,U1,SA)*FC(0,V2,SF)*FC(0,V1,SF)*FC(0,U2,0.00d00)
      ENDIF
    END IF
    END IF
 
    !SET HAMILTONIAN
    H(II,JJ)=FLAGH
    H(JJ,II)=H(II,JJ)
  end do
end do

end subroutine init_H


