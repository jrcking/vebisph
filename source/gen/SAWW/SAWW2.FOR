C*****SAWW
      PROGRAM SAWW
C                     ************
C
C        A PROGRAM TO GENERATE ARBITRARY WATER WAVES
C        USING THE METHOD OF RIENECKER AND FENTON
C
C                     ************
C
      REAL PI,H,T,ETABA,ETAM,H0,CES,MAX0,HSTR
      REAL D,C,K,PRT1,PRT2,FC,FDC,DC,B0,Q,R
      REAL ETA(101),B(100),SUM1,A(226,226)
      REAL F(226),Z(226),dx
      REAL MAX,RHO,CT,X,Y1,UB1,DH1,H1,SIGMA,STEEP,X1
      REAL DH,CSTR,QSTR,RSTR,KSTR,B0STR,ETASTR(101),BSTR(100)
      REAL H01,HMX1,LMDA,CTK,UA,VA,UB,VB,PR,TEMP1,TEMP2,Y(400)
      real temp,ncount,aver_u, pr_old ! calculate averaged horizontal velocity
C
      INTEGER I,J,N1,N2,N3,ND2,ND3,ND4,ND5,ND6,IAT
      INTEGER FLAG,ITERA1,ITERA2,N,NMAX
      INTEGER NT,TSTEP,YMX
      integer iT,itmax, icount, ii
C
      COMMON/BLK1/PI,K,D,ETA,N1
      COMMON/BLK2/B0,B,N
      COMMON/BLK3/Q
      OPEN(5,FILE='SAWW.D')
      OPEN(6,FILE='SAWW.RES')
C
      PI=3.1415927
      G=9.81
      D=1.0
      RHO=1000.0
      DH=0.0
      NMAX=226
      IAT=1
      MAX=0.0
C
C INPUT IN M.K.S. UNITS:-
C THE WAVE HEIGHT,H (PEAK TO TROUGH)
C THE PERIOD,T
C THE WATER DEPTH,ETABA
C THE CURRENT,CES (SET FLAG=0), OR THE MASS TRANSPORT
C     VELOCITY,CES (SET FLAG=1)
C THE NUMBER OF DEPTHS AT WHICH THE VELS. AND ACCELS.
C    ETC. ARE TO BE GIVEN,YMX
C THE NUMBER OF POINTS THROUGHOUT THE CYCLE AT WHICH
C     THE VELS. ACCELS. ETC ARE TO BE GIVEN,NT
C THE NUMBER OF FOURIER COEFF MINUS 1. N LE 30
C
      itmax=1
!      print*,'input period T'
 
!    Time period  
      T=1.6461056608!(h=0.1)
      print*,T
      open(200,file='T_L.dat',access='append')
      do iT=1,itmax

!       T=T+0.001
!   Wave height
       H=0.125*2.0 
!       T=1.56d0
!   Water depth  
       ETABA=1.0d0
       CES=0.0
       FLAG=1

!    Order of Fourier expansion
       N=30

!    Particle spacing
       dx=0.05

!   Number of points in y-direction
       YMX=(ETABA+H+5*dx)/dx
       NT=1200


C
C INPUT THE DEPTHS,Y(YMX)
C
!      READ(5,*)(Y(J),J=1,YMX)
!      rewind(5)

       WRITE(6,218)
       WRITE(6,219)
       WRITE(6,221)DH,H
C
C NONDIMENSIONLISE THE DEPTHS
C
       DO 50 J=1,YMX
   50 Y(J)=Y(J)/ETABA
C
      N1=N+1
      N2=N+2
      N3=N+3
      ND2=2*N+2
      ND3=2*N+3
      ND4=2*N+4
      ND5=2*N+5
      ND6=2*N+6
C
C
C NONDIMENSIONLISE THE INITIAL DATA
C
      H=H/ETABA
      T=T*(G/ETABA)**0.5
      CES=CES/(G*ETABA)**0.5
      H0=H
C
C ITERATION FOR THE INITIAL APPROXIMATION FOR
C C  AND HENCE K. NEWTONS METHOD IS USED.
C
   30 ITERA1=0
      C=1.0
    1 K=2.0*PI/(T*C)
      PRT1=(TANH(K)/K)**0.5
      PRT2=SINH(K)*COSH(K)
      FC=C-PRT1
      FDC=1.0-0.5*PRT1/C*(1.0-K/PRT2)
      DC=-FC/FDC
        IF(ABS(DC).LT.1.0E-8)THEN
      C=C+DC
      K=2.0*PI/(T*C)
       ELSE
      C=C+DC
      ITERA1=ITERA1+1
      GOTO 1
        ENDIF
C
C FURTHER INITIAL VALUES
C
      B0=-C
      Q=C
      R=1.0+0.5*C**2
      B(1)=-0.25*H/(C*K)
      ETA(1)=1.0+0.5*H
      DO 15 J=2,N
      B(J)=0.0
   15 ETA(J)=1.0+0.5*H*COS((J-1.0)*PI/N)
      ETA(N1)=1.0-0.5*H
C
C CALCULATE THE MEAN OF ETA.
C
   33 ITERA2=0
   14 SUM1=0.0
      DO 10 J=2,N
   10 SUM1=SUM1+ETA(J)
      ETAM=0.5/N*(ETA(1)+ETA(N1)+2.0*SUM1)
C
C SET THE ARRAYS A AND F TO 0.
C
      DO 20 J=1,ND6
      F(J)=0.0
      DO 20 I=1,ND6
   20 A(J,I)=0.0
C
C CREATE THE ARRAY A. I-ROWS:J-COLS.
C
      DO 16 I=1,N1
      A(I,I)=U(I)
      DO 2 J=N2,ND2
    2 A(I,J)=DFB1N1(I,J-N2)
      A(I,ND4)=DFK1N1(I)
   16 A(I,ND5)=1.0
C
      DO 3 I=N2,ND2
      A(I,I-N1)=DFNN2(I-N1)
      DO 4 J=N2,ND2
    4 A(I,J)=DFBN2(I-N1,J-N2)
      A(I,ND4)=DFKN2(I-N1)
    3 A(I,ND6)=-1.0
C
      A(ND3,1)=0.5/N
      A(ND3,N1)=0.5/N
      DO 5 J=2,N
    5 A(ND3,J)=1.0/N
C
      A(ND4,1)=1.0
      A(ND4,N1)=-1.0
C
      A(ND5,ND3)=K*T
      A(ND5,ND4)=C*T
C
        IF(FLAG.EQ.0)THEN
      A(ND6,N2)=1.0
       ELSE
      A(ND6,ND5)=-1.0
        ENDIF
      A(ND6,ND3)=1.0
C
C CREATE THE ARRAY F.
C
      DO 7 I=1,N1
    7 F(I)=F1N1(I)
      DO 8 I=N2,ND2
    8 F(I)=0.5*(U(I-N1)**2+V(I-N1)**2)+ETA(I-N1)-R
      F(ND3)=ETAM-1.0
      F(ND4)=ETA(1)-ETA(N1)-H
      F(ND5)=K*C*T-2.0*PI
        IF(FLAG.EQ.0)THEN
      F(ND6)=C-CES+B0
       ELSE
      F(ND6)=C-CES-Q
        ENDIF
      DO 21 J=1,ND6
   21 F(J)=-1.0*F(J)
C
C SOLVE THE SIMULTANEOUS EQATIONS F=AZ.
C
      CALL SOLVE(A,NMAX,Z,F,ND6)
C
C
C SEARCH FOR LARGEST VALUE OF Z
C
      MAX0=MAX
      MAX=Z(1)
      DO 11 J=2,ND6
   11 MAX=AMAX1(ABS(MAX),ABS(Z(J)))
      WRITE(6,220)MAX
      C=C+Z(ND3)
      K=K+Z(ND4)
      Q=Q+Z(ND5)
      R=R+Z(ND6)
      B0=B0+Z(N2)
      DO 12 J=1,N1
   12 ETA(J)=ETA(J)+Z(J)
      DO 13 J=N3,ND2
   13 B(J-N2)=B(J-N2)+Z(J)
      ITERA2=ITERA2+1
C
C ALLOW THE ITERATION TO SETTLE
C
      IF(ITERA2.LT.3)GOTO 14
C
      IF(MAX.GT.MAX0)THEN
C
C THE ITERATION IS NOT CONVERGING SO
C
       IF(IAT.EQ.1)THEN
C
C NO PREVIOUSLY CONVERGED SOLUTION EXISTS. HALVE THE
C WAVE HEIGHT AND TRY AGAIN.
C
       DH=H/2
      IF(DH.LT.(0.01*H0))THEN
      WRITE(6,227)
      WRITE(6,228)
      STOP
      ENDIF
       H=H-DH
      DH1=DH*ETABA
      H1=H*ETABA
        WRITE(6,221)DH1,H1
        GOTO 30
       ELSE
C
C A PREVIOUSLY CONVERGED SOLUTION DOES EXIST SO USE
C THIS AS THE STARTING POINT FOR A NEW ATTEMPT TO
C REACH THE REQUESTED H0. HALVE THE WAVE HEIGHT STEP
C INCREMENT,DH, ADD TO PREVIOUS H AND START AGAIN
C
       C=CSTR
       K=KSTR
       Q=QSTR
       R=RSTR
       B0=B0STR
       DO 31 J=1,N1
   31 ETA(J)=ETASTR(J)
       DO 32 J=1,N
   32 B(J)=BSTR(J)
       DH=DH/2
      IF(DH.LT.(0.01*H0))THEN
      WRITE(6,227)
      WRITE(6,229)
      WRITE(6,230)
      GOTO 500
      ENDIF
       H=HSTR+DH
       DH1=DH*ETABA
       H1=H*ETABA
        WRITE(6,221)DH1,H1
        GOTO 33
       ENDIF
      ELSE
C
       IF(MAX.LT.1.0E-6)THEN
C
C THE SOLUTION HAS CONVERGED FOR SOME H
C
        IF((H0-H).LT.1.0E-5)THEN
C
C THE SOLUTION HAS CONVERGED FOR THE REQUESTED H
C
      GOTO 500
        ELSE
C
C ESTIMATE THE MAXIMUM WAVE HEIGHT
C
         HMX=0.284*PI/K*TANH(K)
         IF(H.GT.0.9*HMX)THEN
         IF(H0.GT.HMX)THEN
C REQUESTED H EXCEEDS THE MAXIMUM
         H01=H0*ETABA
         HMX1=HMX*ETABA
         WRITE(6,200)H01,HMX1
         STOP
         ELSE
C REQUESTED H WITHIN 10% OF MAXIMUM
         HMX1=HMX*ETABA
         WRITE(6,201)HMX1
         GOTO 500
         ENDIF
         ENDIF
C
        IAT=2
        CSTR=C
        KSTR=K
        QSTR=Q
        RSTR=R
        B0STR=B0
        HSTR=H
        DO 34 J=1,N1
   34 ETASTR(J)=ETA(J)
        DO 35 J=1,N
   35 BSTR(J)=B(J)
        H=H+DH
      DH1=DH*ETABA
      H1=H*ETABA
         WRITE(6,221)DH1,H1
         GOTO 33
        ENDIF
C
       ELSE
       GOTO 14
       ENDIF
      ENDIF
C
C RESTORE THE DIMENSIONS
C
  500 CSTR=C*(G*ETABA)**0.5
      HSTR=H*ETABA
      KSTR=K/ETABA
      LMDA=2.0*PI/KSTR
      T=T*(ETABA/G)**0.5
      SIGMA=ETABA/LMDA
      STEEP=HSTR/LMDA
      IF((SIGMA.LE.0.05).AND.(STEEP.GE.0.01))THEN
      WRITE(6,210)
      WRITE(6,226)
      WRITE(6,210)
      ENDIF
      DO 36 J=1,N1
   36 ETASTR(J)=ETA(J)*ETABA
      CES=CES*(G*ETABA)**0.5
      WRITE(6,202)
      WRITE(6,203)
      WRITE(6,204)
      WRITE(6,205)HSTR,CSTR,LMDA,KSTR,T,CES,ETABA
      WRITE(6,206)
      WRITE(6,207)
      WRITE(6,208)
      DO 37 I=1,N1
      X=LMDA*(I-1.0)/(N*2.0)
   37 WRITE(6,209)X,ETASTR(I)
      WRITE(6,210)
      X=0.0
C
C CALCULATE THE VELS AND ACCELS AT THE CREST
C
      UA=0.0
      UB=0.0
      VA=0.0
      VB=0.0
      DO 41 J=1,N
      UB=UB+J*B(J)*COSH(J*K*ETA(1))/COSH(J*K*D)
   41 VA=VA-J*B(J)*SINH(J*K*ETA(1))/COSH(J*K*D)*J*K*C
      UB=(B0+C+K*UB)*(G*ETABA)**0.5
      VA=K*VA*G
      WRITE(6,222)
      WRITE(6,223)
      WRITE(6,224)
      WRITE(6,225)UB,VB,UA,VA
C
C CALCULATE ALL OTHER VELS AND ACCELS ETC
C
!       open(20,file='aver_U.dat')
 
    !Write positions and velocities to this file
      open(111,file='tecplot.dat')

1101  format(1x,1es18.10e3,1x,1es18.10e3,1x,
     &       1es18.10e3,1x,1es18.10e3,1x,1es18.10e3) 

       CT=(dx/2d0)/ETABA
      ! DO 38 TSTEP=1,NT
     
      DO while ( CT .lt. 3.0*LMDA/ETABA )

      !CT=4.0*PI*(TSTEP-1.0)/(K*NT)
        CT=CT+dx/ETABA

       CTK=CT*K
      WRITE(6,214)         ! comment cout to reduce the file size for aver_u output
      WRITE(6,215)CTK
      WRITE(6,211)
      WRITE(6,212)
      WRITE(6,213)

!first time round to find precise free surface location
      y=0d0
      do ii=1,YMX
        y(ii)=(ii-1d0)*dx/etaba
      enddo


      temp=0.0
      ncount=0
      DO I=1,YMX
      UA=0.0
      UB=0.0
      VA=0.0
      VB=0.0
      DO 40 J=1,N
      TEMP1=J*B(J)*COSH(J*K*Y(I))/COSH(J*K*D)
      TEMP2=J*B(J)*SINH(J*K*Y(I))/COSH(J*K*D)
      UB=UB+TEMP1*COS(J*K*(X+CT))
      VB=VB+TEMP2*SIN(J*K*(X+CT))
      UA=UA+TEMP1*SIN(J*K*(X+CT))*J*C*K
   40 VA=VA-TEMP2*COS(J*K*(X+CT))*J*K*C
      UB1=(B0+K*UB)*(G*ETABA)**0.5
      UB=(B0+C+K*UB)*(G*ETABA)**0.5
      VB=K*VB*(G*ETABA)**0.5
      UA=K*UA*G
      VA=K*VA*G
      PR=((R-Y(I))*G*ETABA-0.5*(UB1**2+VB**2))*RHO
      Y1=Y(I)*ETABA
      X1=(X+CT)*ETABA
      IF(PR.le.0.0)THEN
!        WRITE(6,217)Y1
 
        if( PR_OLD .ge. 0.0 ) then
            y(i) = y(i-1) + 
     &        (-pr_old*y(i)+pr_old*y(i-1))/ 
     &        (pr-pr_old)
           do ii=i-1, 2, -1
             y(ii) = y(ii+1)-dx/etaba
           end do
           ii = i
           exit
        endif
      ENDIF
        pr_old=pr
      end do

!second time round to write data to file

      temp=0.0
      ncount=0
      DO I=1,ii
      UA=0.0
      UB=0.0
      VA=0.0
      VB=0.0
      DO 43 J=1,N
      TEMP1=J*B(J)*COSH(J*K*Y(I))/COSH(J*K*D)
      TEMP2=J*B(J)*SINH(J*K*Y(I))/COSH(J*K*D)
      UB=UB+TEMP1*COS(J*K*(X+CT))
      VB=VB+TEMP2*SIN(J*K*(X+CT))
      UA=UA+TEMP1*SIN(J*K*(X+CT))*J*C*K
   43 VA=VA-TEMP2*COS(J*K*(X+CT))*J*K*C
      UB1=(B0+K*UB)*(G*ETABA)**0.5
      UB=(B0+C+K*UB)*(G*ETABA)**0.5
      VB=K*VB*(G*ETABA)**0.5
      UA=K*UA*G
      VA=K*VA*G
      PR=((R-Y(I))*G*ETABA-0.5*(UB1**2+VB**2))*RHO
      Y1=Y(I)*ETABA
      X1=(X+CT)*ETABA

!      write(6,*) y1, ub, vb, pr
      write(111,*) x1, y1, ub, vb, pr
    
      end do




!       aver_u=temp/real(ncount)
!       write(20,*)real(TSTEP)*T/real(NT),aver_u,ncount,
!      &           real(ncount)*0.6/real(YMX)
      WRITE(6,210)
      enddo

      close(111)

!       close(20)
  200 FORMAT(/10X,'REQUESTED WAVE HEIGHT,',F5.2,1X,
     +',EXCEEDS THE MAXIMUM,',F5.2)
  201 FORMAT(/10X,'REQUESTED WAVE IS WITHIN 10% OF MAX,',F5.2)
  202 FORMAT(//4X,4('WAVE',4X),'PERIOD',4X,'CURRENT OR',4X,'DEPTH')
  203 FORMAT(1H ,2X,'HEIGHT',3X,'SPEED',2X,'LENGTH',2X
     +,'NUMBER',11X,'MASS TRANSPORT')
  204 FORMAT(1H ,3X,'(M)',5X,'(M/S)',3X,'(M)',4X,'(RAD/M)',
     +4X,'(S)',9X,'(M/S)',7X,'(M)')
  205 FORMAT(/4X,F5.2,3X,F5.2,2X,F7.3,2X,F6.3,2X,F5.2,5X,F7.3,5X,F7.3)
  206 FORMAT(//4X,'HORIZONTAL',2X,'SURFACE')
  207 FORMAT(1H ,4X,'DISTANCE',2X,'ELEVATION')
  208 FORMAT(1H ,5X,'(M)',9X,'(M)'/)
  209 FORMAT(1H ,5X,F6.3,6X,F7.4)
  210 FORMAT(//)
  211 FORMAT(1H ,2('VERTICAL',3X,'HORIZONTAL',3X),'VERTICAL',3X,
     +'PRESSURE')
  212 FORMAT(1H ,'DISTANCE',2(4X,'VELOCITY'),6X,2('ACCEL',7X))
  213 FORMAT(1H ,3X,'(M)',8X,'(M/S)',7X,'(M/S)',7X,'(M/SS)',6X,'(M/SS)',
     +5X
     +,'(N/MM)')
  214 FORMAT(1H ,30X,'PHASE,(RAD)')
  215 FORMAT(1H ,30X,F5.3/)
  216 FORMAT(1H ,F12.7,2(3X,E14.7),2(2X,E14.7),1X,E14.7)
  217 FORMAT(1H ,F8.3,5X,'SELECTED HEIGHT LIES ABOVE THE WAVE')
  218 FORMAT(////7X,'ERROR',7X,'STEP     WAVE')
  219 FORMAT(1H ,17X,'LENGTH    HEIGHT'/)
  220 FORMAT(1H ,3X,E10.3)
  221 FORMAT(1H ,16X,2(F7.3))
  222 FORMAT(1H ,10X,'VELS AND ACCELS AT THE CREST'/)
  223 FORMAT(1H ,6X,'VELOCITY',18X,'ACCELERATIONS'/9X,'(M/S)',22X,
     +'(M/SS)')
  224 FORMAT(1H ,2X,'HORI',8X,'VERT',11X,'HORI',8X,'VERT'/)
  225 FORMAT(1H ,E10.3,2X,E10.3,5X,E10.3,2X,E10.3//)
  226 FORMAT(1H ,10X,'WARNING! THE WAVE APPROACHES ',
     +'THE SOLITARY'/21X,'WAVE LIMIT. TEST ACCURACY BY '
     +/21X,'INCREASING N.')
  227 FORMAT(1H ,'THE STEP LENGTH IS LESS THAN 1% OF ',
     +'THE'/2X,'REQUESTED WAVE HEIGHT')
  228 FORMAT(1H ,'NO CONVERGED SOLUTION HAS BEEN OBTAINED ',
     +/2X,'AT ANY WAVE HEIGHT.TRY INCREASING N, IF POSSIBLE.')
  229 FORMAT(1H ,'RESULTS FOR THE LAST CONVERGED SOLUTION FOLLOW')
  230 FORMAT(1H ,'TRY INCREASING N, IF POSSIBLE.')
      write(200,*)T,LMDA
      enddo ! end of iT=1,itmax
      close(200)
      STOP
      END
      REAL FUNCTION DFBN2(M,J)
      REAL K,PI,D,ETA(101)
      INTEGER J,M,N1
      COMMON/BLK1/PI,K,D,ETA,N1
        IF(J.EQ.0)THEN
      DFBN2=U(M)
       ELSE
      DFBN2=U(M)*J*K*C2(M,J)+V(M)*J*K*S2(M,J)
        ENDIF
      RETURN
      END
      REAL FUNCTION S1(M,J)
      REAL ETA(101),K,D,PI
      INTEGER J,M,N1
      COMMON/BLK1/PI,K,D,ETA,N1
      S1=SINH(J*K*ETA(M))/COSH(J*K*D)*COS(J*(M-1.0)*PI/(N1-1.0))
      RETURN
      END
      REAL FUNCTION S2(M,J)
      REAL ETA(101),PI,K,D
      INTEGER J,M,N1
      COMMON/BLK1/PI,K,D,ETA,N1
      S2=SINH(J*K*ETA(M))/COSH(J*K*D)*SIN(J*(M-1.0)*PI/(N1-1.0))
      RETURN
      END
      REAL FUNCTION C1(M,J)
      REAL ETA(101),PI,D,K
      INTEGER J,M,N1
      COMMON/BLK1/PI,K,D,ETA,N1
      C1=COSH(J*K*ETA(M))/COSH(J*K*D)*SIN(J*(M-1.0)*PI/(N1-1.0))
      RETURN
      END
      REAL FUNCTION C2(M,J)
      REAL ETA(101),D,PI,K
      INTEGER J,M,N1
      COMMON/BLK1/PI,K,D,ETA,N1
      C2=COSH(J*K*ETA(M))/COSH(J*K*D)*COS(J*(M-1.0)*PI/(N1-1.0))
      RETURN
      END
      REAL FUNCTION U(M)
      REAL B0,B(100),PI,K,D,ETA(101),SUM
      INTEGER J,N,M,N1
      COMMON/BLK1/PI,K,D,ETA,N1
      COMMON/BLK2/B0,B,N
      SUM=0.0
      DO 1 J=1,N
    1 SUM=J*B(J)*C2(M,J)+SUM
      U=B0+K*SUM
      RETURN
      END
      REAL FUNCTION V(M)
      REAL B0,B(100),PI,K,D,ETA(101),SUM
      INTEGER J,M,N1,N
      COMMON/BLK1/PI,K,D,ETA,N1
      COMMON/BLK2/B0,B,N
      SUM=0.0
      DO 1 J=1,N
    1 SUM=SUM+J*B(J)*S2(M,J)
      V=K*SUM
      RETURN
      END
      REAL FUNCTION DFB1N1(M,J)
      REAL ETA(101),K,D,PI
      INTEGER J,M,N1
      COMMON/BLK1/PI,K,D,ETA,N1
        IF(J.EQ.0)THEN
      DFB1N1=ETA(M)
       ELSE
      DFB1N1=S1(M,J)
        ENDIF
      RETURN
      END
      REAL FUNCTION DFK1N1(M)
      REAL SUM,ETA(101),B(100),K,D,PI,B0
      INTEGER N,J,N1,M
      COMMON/BLK1/PI,K,D,ETA,N1
      COMMON/BLK2/B0,B,N
      SUM=0.0
      DO 1 J=1,N
    1 SUM=SUM+B(J)*TANH(J*K*D)*S1(M,J)*J*D
      DFK1N1=(U(M)-B0)/K*ETA(M)-SUM
      RETURN
      END
      REAL FUNCTION DFKN2(M)
      REAL B(100),B0,K,D,ETA(101),SUM1,SUM2,SUM3,SUM4,PI
      INTEGER J,M,N,N1
      COMMON/BLK1/PI,K,D,ETA,N1
      COMMON/BLK2/B0,B,N
      SUM1=0.0
      SUM2=0.0
      SUM3=0.0
      SUM4=0.0
      DO 1 J=1,N
      SUM1=J**2*B(J)*S1(M,J)+SUM1
      SUM2=J**2*B(J)*TANH(J*K*D)*C2(M,J)+SUM2
      SUM3=J**2*B(J)*C1(M,J)+SUM3
    1 SUM4=J**2*B(J)*S2(M,J)*TANH(J*K*D)+SUM4
      SUM1=(U(M)-B0)/K+K*ETA(M)*SUM1-K*D*SUM2
      SUM3=V(M)/K+K*ETA(M)*SUM3-K*D*SUM4
      DFKN2=U(M)*SUM1+V(M)*SUM3
      RETURN
      END
      REAL FUNCTION DFNN2(M)
      REAL SUM1,SUM2,K,B(100),B0,ETA(101),PI,D
      INTEGER J,M,N1,N
      COMMON/BLK1/PI,K,D,ETA,N1
      COMMON/BLK2/B0,B,N
      SUM1=0.0
      SUM2=0.0
      DO 1 J=1,N
      SUM1=J**2*B(J)*S1(M,J)+SUM1
    1 SUM2=J**2*B(J)*C1(M,J)+SUM2
      SUM1=SUM1*U(M)*K**2
      SUM2=SUM2*V(M)*K**2
      DFNN2=SUM1+SUM2+1.0
      RETURN
      END
      REAL FUNCTION F1N1(M)
      REAL Q,ETA(101),B0,B(100),PI,K,D,SUM
      INTEGER J,M,N,N1
      COMMON/BLK1/PI,K,D,ETA,N1
      COMMON/BLK2/B0,B,N
      COMMON/BLK3/Q
      SUM=0.0
      DO 1 J=1,N
    1 SUM=SUM+B(J)*S1(M,J)
      F1N1=B0*ETA(M)+SUM+Q
      RETURN
      END
C
      SUBROUTINE SOLVE(K,IK,U,F,N)
C
C   THIS SUBROUTINE PERFORMS GAUSSIAN ELIMINATION WITH
C   PARTIAL PIVOTING ON A FULL N*N MATRIX
C
      REAL K(IK,*),F(*),U(*)
C
C   PIVOTING STAGE
C
      DO 1 I=1,N-1
      BIG=ABS(K(I,I))
      IHOLD=I
      DO 10 J=I+1,N
      IF(ABS(K(J,I)).GT.BIG) THEN
      BIG=ABS(K(J,I))
      IHOLD=J
      END IF
 10   CONTINUE
      IF(IHOLD.NE.I) THEN
      DO 12 J=I,N
      HOLD=K(I,J)
      K(I,J)=K(IHOLD,J)
      K(IHOLD,J)=HOLD
 12   CONTINUE
      HOLD=F(I)
      F(I)=F(IHOLD)
      F(IHOLD)=HOLD
      END IF
C
C   ELIMINATION STAGE
C
       DO 3 J=I+1,N
       FAC=K(J,I)/K(I,I)
       DO 4 L=I,N
 4     K(J,L)=K(J,L)-K(I,L)*FAC
       F(J)=F(J)-F(I)*FAC
 3     CONTINUE
 1     CONTINUE
C
C   BACK SUBSTITUTION STAGE
C
       DO 9 I=N,1,-1
       SUM=0.0
       DO 6 L=I+1,N
 6     SUM=SUM+K(I,L)*U(L)
       U(I)=(F(I)-SUM)/K(I,I)
 9     CONTINUE
       RETURN
       END
