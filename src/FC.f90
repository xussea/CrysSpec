!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!% CALCULATE THE GENERAL FRANCK-CONDON FACTOR <U|V> %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

FUNCTION FC(W,L,S)
IMPLICIT NONE

INTEGER :: W,L,U,V,II
DOUBLE PRECISION ::FC,S,LAMTA
DOUBLE PRECISION, EXTERNAL :: FAC

IF ( S < 0.0 ) then
  U=L
  V=W
ELSE
  U=W
  V=L
END IF

LAMTA=dsqrt(abs(S))
FC=0.0D0
!IF ( V .GT. U ) THEN
 DO II=0,min(U,V)
   FC=FC+(-1.0)**(V-II)*LAMTA**(V+U-2*II)/(FAC(II)*FAC(U-II)*FAC(V-II))
 ENDDO
!ELSE
!  DO II=0,V
!    FC=FC+(-1.0)**(V-II)*LAMTA**(V+U-2*II)/(FAC(II)*FAC(U-II)*FAC(V-II))
!  ENDDO
!END IF
FC=DEXP(-LAMTA**2/2.0)*DSQRT(FAC(V)*FAC(U))*FC
ENDFUNCTION


