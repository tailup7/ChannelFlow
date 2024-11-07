program main 
    PARAMETER (NX=64,NY=64,NZ=64,NY1=NY+1)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    CHARACTER*3  CT
    CHARACTER*4  FILE1
    COMMON /FILE1/FILE1 /STEP/ISTEP /TSTEP/TSTEP
    COMMON /U/U(NZ,NX,NY) /V/V(NZ,NX,0:NY) /W/W(NZ,NX,NY)
    COMMON /P/P(NZ,NX,NY) 
    ct='000'
    FILe1='dns2'
    ISEP=0
    TSTEP=0.0000
    DO 10 J=1,NY
        DO 10 I=1,NX
        DO 10 K=1,NZ
        U(K,I,J)=1
        W(K,I,J)=0.01*(65-I) 
        P(K,I,J)=100-I
 10 CONTINUE

    DO 20 J=1,NY-1
        DO 20 I=1,NX
        DO 20 K=1,NZ
        V(K,I,J)=0              
 20 CONTINUE

    OPEN(10,file=FILE1//CT//'u.d',status='new',form='formatted')            
      WRITE(10,1000) ISTEP,TSTEP
      WRITE(10,2000) U
      CLOSE(10)
    OPEN(10,file=FILE1//CT//'v.d',status='new',form='formatted')
      WRITE(10,1000) ISTEP,TSTEP
      WRITE(10,2000) V
      CLOSE(10)
    OPEN(10,file=FILE1//CT//'w.d',status='new',form='formatted')
      WRITE(10,1000) ISTEP,TSTEP
      WRITE(10,2000) W
      CLOSE(10)
    OPEN(10,file=FILE1//CT//'p.d',status='new',form='formatted')
      WRITE(10,1000) ISTEP,TSTEP
      WRITE(10,3000) P
      CLOSE(10)
1000 FORMAT(I10,F12.6)
2000 FORMAT(8F10.6)
3000 FORMAT(8F10.5)
    WRITE(6,1200) ISTEP,TSTEP
1200 FORMAT(/'Data SAVE',10X,'STEP=',I7,10X,'T=',F10.3/)
    RETURN
    END