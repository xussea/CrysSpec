subroutine calculate_HR(QF,QC,QA,SF,SC,SA,SCA,SFC,SFA)
implicit none
double precision :: z,QF,QC,QA,SF,SC,SA,SCA,SFC,SFA

if ( abs(QF) > 1.0d-06 ) then
        z=QF/abs(QF)
        SF=z*QF**2/2
else 
        SF=0.00d00
end if

if ( abs(QA) > 1.0d-06 ) then
        z=QA/abs(QA)
        SA=z*QA**2/2
else 
        SA=0.00d00
end if

if ( abs(QC) > 1.0d-06 ) then
        z=QC/abs(QC)
        SC=z*QC**2/2
else
        SC=0.00d00
end if

if ( abs(QC-QA) > 1.0d-06 ) then
        z=(QC-QA)/abs(QC-QA)
        SCA=z*(QC-QA)**2/2
else 
        SCA=0.00d00
end if

if ( abs(QF-QA) > 1.0d-06 ) then
        z=(QF-QA)/abs(QF-QA)
        SFA=z*(QF-QA)**2/2
else
        SFA=0.00d00
end if

if ( abs(QF-QC) > 1.0d-06 ) then
        z=(QF-QC)/abs(QF-QC)
        SFC=z*(QF-QC)**2/2
else
        SFC=0.00d00  
end if
end subroutine calculate_HR


