!*****************************************************************
!     cgaus.f                2012-12-04  
!     include two subroutine:
!      1) cfgaus.f--mult-integrate, complex
!      2) gauss_oset.f
!****************************************************************
!
!=========================================================
!--------------cfgaus.f------------------------------
!=========================================================
	SUBROUTINE FGAUS(N,JS,X,FS,F,S)
	DIMENSION JS(N),X(N)
	DIMENSION T(5),C(5),IS(2,11),D1(11),CC(11)
	DOUBLE PRECISION T,C,DN,UP,CC,D1,X
        complex*16 :: F,S,P
        complex*16,dimension(11) :: D2
	DATA T/-0.9061798459,-0.5384693101,0.0,
     *         0.5384693101,0.9061798459/
	DATA C/0.2369268851,0.4786286705,0.5688888889,
     *         0.4786286705,0.2369268851/
	M=1
	D1(N+1)=1.0
	D2(N+1)=1.0
10	DO 20 J=M,N
	  CALL FS(J,N,X,DN,UP)
	  D1(J)=0.5*(UP-DN)/JS(J)
	  CC(J)=D1(J)+DN
	  X(J)=D1(J)*T(1)+CC(J)
	  D2(J)=0.0
	  IS(1,J)=1
	  IS(2,J)=1
20	CONTINUE

	J=N
30	K=IS(1,J)
	IF (J.EQ.N) THEN
	  P=F(N,X)
	ELSE
	  P=1.0
	END IF
	D2(J)=D2(J+1)*D1(J+1)*P*C(K)+D2(J)
	IS(1,J)=IS(1,J)+1
	IF (IS(1,J).GT.5) THEN
	  IF (IS(2,J).GE.JS(J)) THEN
	    J=J-1
	    IF (J.EQ.0) THEN
	      S=D2(1)*D1(1)
	      RETURN
	    END IF
	    GOTO 30
	  END IF
	  IS(2,J)=IS(2,J)+1
	  CC(J)=CC(J)+D1(J)*2.0
	  IS(1,J)=1
	END IF
	K=IS(1,J)
	X(J)=D1(J)*T(K)+CC(J)
	IF (J.EQ.N) GOTO 30
	M=J+1
	GOTO 10
	END

!=========================================================
!--------------gauss_oset.f------------------------------
!=========================================================

!=======EXAMPLE =========================================
!      subroutine gpropagator(s,xm1,xm2,g)
!      IMPLICIT NONE

!        common/loop/s_xie
!	common/xmass/xm1_xie,xm2_xie
!	
!        COMPLEX*16, INTENT(IN) :: s
!	COMPLEX*16  :: g,s_xie
!        REAL*8 :: xm1,xm2,xm1_xie,xm2_xie,a1,b
!	integer :: n
!	COMPLEX*16, EXTERNAL :: kernel
      
!	a1=0.d0
!	b=630.d0
!        n=100

!	s_xie=s
!	xm1_xie=xm1
!	xm2_xie=xm2

!	call gaussc(kernel,a1,b,n,g)

!      RETURN
!      END subroutine gpropagator
!
!----------------------------------------------------
!
!      function kernel(q)
!      IMPLICIT NONE

!	common/loop/s_xie
!	common/xmass/xm1_xie,xm2_xie
	
!        COMPLEX*16  :: s_xie,s
!        REAL*8 :: pie,epsilon
!	real*8 :: q,omega,el,xm1_xie,xm2_xie,xm1,xm2
!	complex*16 :: kernel
   
!	pie=4.d0*atan(1.d0)
!	epsilon=0.2d0
      
!	s=s_xie
!	xm1=xm1_xie     ! mass of meson
!	xm2=xm2_xie     ! mass of baryon

!    	omega=sqrt(xm1**2+q**2)
!	el=sqrt(xm2**2+q**2)

!	kernel=4.d0*pie*q**2*xm2/(2.d0*omega)/el/(sqrt(s)-omega-el
!     f +(0.d0,1.d0)*epsilon)/(2.d0*pie)**3

!      RETURN
!      END function kernel

      SUBROUTINE gaussc(kernel,a,b,n,t)
      IMPLICIT NONE
      complex*16,dimension(2000)  :: f1
      real*8,dimension(2000) ::x
      complex*16, EXTERNAL :: kernel
      integer :: i,n,np
      real*8 :: a,b
      complex*16 :: t

      call DSG20r(a,b,n,x,np)
      do i=1,np
      f1(i)=kernel(x(i))
      end do

      call DRG20c(a,b,n,f1,t)
     
      return
      END
     
      SUBROUTINE gaussr(kernel,a,b,n,t)
      IMPLICIT NONE
      real*8,dimension(2000)  :: x,f1
      real*8, EXTERNAL :: kernel
      integer :: i,n,np
      real*8 :: a,b,t

      call DSG20r(a,b,n,x,np)
      do i=1,np
      f1(i)=kernel(x(i))
      end do

      call DRG20r(a,b,n,f1,t)
     
      return
      END

!*******************************************************************************
!*  Proposito:        Integracion numerica por Gauss (funciones reales y complejas)
!*
!*      SUBROUTINE SG20R(A,B,N,X,NP)----->Realiza la particion
!*                Entradas: A--->limite inferior
!*                          B--->limite superior
!*                          N--->Numero de puntos Gauss/20 (n=1 gauss 20 puntos
!*                                                        (n=2 gauss 40 puntos,..)
!*                Salidas:  X--->Particion de la variable (ha de ser dimensionada)
!*                          NP-->Numero de puntos gauss (20*n)
!*
!*      SUBROUTINE RG20C(A,B,N,CF,CRES)--->Integracion propiamente dicha
!*                Entradas: A,B,N-->Idem antes
!*                          CF----->Valor de la funcion a integrar en los
!*                                  puntos X, CF(i)=funcion(x(i))
!*                                  ha de ser dimensionada
!*                Salidas:  CRES--->Resultado de la integral
!*                                  En caso de integracion en varias variables,
!*                                  CRES pasa a ser la entrada CF de la siguiente
!*                                  integral.
!*
!*
!*      SUBROUTINE DSG20R(A,B,N,X,NP)
!*      SUBROUTINE DRG20C(A,B,N,CF,CRES)
!*               Idem doble precision
!*
!*
!*
!*      SUBROUTINE  RG20R(A,B,N,F,RES)
!*      SUBROUTINE DRG20R(A,B,N,F,RES)
!*                        Idem funciones reales
!*
!* El programa esta preparado para hacer hasta 2000 puntos gauss.
!* Para aumentarlo hay que cambiar las dimensiones de las variables X(2000)
!* y CF(2000) al valor deseado.
!    *******************************************************************
!               SUBROUTINAS DE INTEGRACION NUMERICA POR GAUSS (Double precision)
!    *******************************************************************

      SUBROUTINE DSG20R(A,B,N,X,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(10),X(2000)
      DATA Y/.9931285991d0,.9639719272d0,.9122344282d0,.8391169718d0,
     F .7463319064d0,.6360536807d0,.5108670019d0,.3737060887d0,
     F .2277858511d0,.0765265211d0/
      
	NP=20*N
      DINT=(B-A)/dble(N)
      DELT=DINT*0.5d0
      ORIG=A-DELT
      I1=-20
      DO 1 I=1,N
      ORIG=ORIG+DINT
      DORIG=ORIG+ORIG
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
      X(J1)=ORIG-DELT*Y(J)
 2    X(J2)=DORIG-X(J1)
 1    CONTINUE
      RETURN
      END

      SUBROUTINE DRG20C(A,B,N,CF,CRES)
      IMPLICIT REAL*8 (A,B,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      IMPLICIT complex*16 (C)
      DIMENSION W(10),CF(2000)
      DATA W/.0176140071d0,.0406014298d0,.0626720483d0,.0832767415d0,
     F .1019301198d0,.1181945319d0,.1316886384d0,.1420961093d0,
     f .1491729864d0,.1527533871d0/
      CR=(0.d0,0.d0)
      I1=-20
      DO 1 I=1,N
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
 2    CR=CR+W(J)*(CF(J1)+CF(J2))
 1    CONTINUE
      CRES=CR*0.5d0*(B-A)/dble(N)
      RETURN
      END

      SUBROUTINE DRG20R(A,B,N,CF,CRES)
      IMPLICIT REAL*8 (A,B,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)
      IMPLICIT real*8 (C)
      DIMENSION W(10),CF(2000)
      DATA W/.0176140071d0,.0406014298d0,.0626720483d0,.0832767415d0,
     F .1019301198d0,.1181945319d0,.1316886384d0,.1420961093d0,
     f .1491729864d0,.1527533871d0/
      CR=(0.d0,0.d0)
      I1=-20
      DO 1 I=1,N
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
 2    CR=CR+W(J)*(CF(J1)+CF(J2))
 1    CONTINUE
      CRES=CR*0.5d0*(B-A)/dble(N)
      RETURN
      END


!--------1---------2---------3---------4---------5---------6--------7--
!#NUMPAC#MINVB               REVISED ON 1984-11-30
!--------1---------2---------3---------4---------5---------6--------7--

      SUBROUTINE MINVB(A,KA,N,EPS,ILL)

      REAL*8 EPS,AM,AA

      COMPLEX*16 A(KA,N),S,W,P

      INTEGER*2 MX(1000)

      IF(N.LT.1.OR.N.GT.1000.OR.N.GT.KA.OR.EPS.LE.0.d0) GO TO 250

!-----LUDECOMPOSITION------------------------------------------

      NM1=N-1

      DO 90 J=1,N

      IF(J.EQ.1) GO TO 30

      JM1=J-1

      DO 20 I=1,JM1

      M=MX(I)

      S=A(M,J)

      A(M,J)=A(I,J)

      IF(I.EQ.1) GO TO 20

      IM1=I-1

      DO 10 K=1,IM1

   10 S=A(I,K)*A(K,J)+S

   20 A(I,J)=S

   30 AM=0.d0

      DO 60 I=J,N

      S=A(I,J)

      IF(J.EQ.1) GO TO 50

      DO 40 K=1,JM1

   40 S=A(I,K)*A(K,J)+S

      A(I,J)=S

   50 AA=CDABS(S)

      IF(AA.LE.AM) GO TO 60

      AM=AA

      M=I

   60 CONTINUE

      IF(AM.LT.EPS) GO TO 240

      MX(J)=M

      IF(M.EQ.J) GO TO 80

      DO 70 K=1,J

      W=A(M,K)

      A(M,K)=A(J,K)

   70 A(J,K)=W

   80 IF(J.EQ.N) GO TO 100

      JP1=J+1

      W=-A(J,J)

      DO 90 I=JP1,N

   90 A(I,J)=A(I,J)/W

  100 IF(N.LE.2) GO TO 130

!-----INPLACE INVERSION OF L-COMPONENT---------------------------

      DO 120 I=3,N

      IM1=I-1

      IM2=I-2

      DO 120 J=1,IM2

      S=A(I,J)

      JP1=J+1

      DO 110 K=JP1,IM1

  110 S=A(I,K)*A(K,J)+S

  120 A(I,J)=S

!-----INPLACE INVERSION OF U-COMPONENT----------------------------

  130 A(1,1)=(1.d0,0.d0)/A(1,1)

      IF(N.EQ.1) GO TO 230

      DO 150 J=2,N

      A(J,J)=(1.d0,0.d0)/A(J,J)

      P=-A(J,J)

      JM1=J-1

      DO 150 I=1,JM1

      S=0.d0

      DO 140 K=I,JM1

  140 S=A(I,K)*A(K,J)+S

  150 A(I,J)=S*P

!-----INPLACE MULTIPLICATION OF L AND U COMPONENT-----------------

      DO 190 J=1,NM1

      JP1=J+1

      DO 170 I=1,J

      S=A(I,J)

      DO 160 K=JP1,N

  160 S=A(I,K)*A(K,J)+S

  170 A(I,J)=S

      DO 190 I=JP1,N

      S=0.d0

      DO 180 K=I,N

  180 S=A(I,K)*A(K,J)+S

  190 A(I,J)=S

!------INTERCHANGE OF COLUMNS------------------------------------

      J=NM1

  200 M=MX(J)

      IF(M.EQ.J) GO TO 220

      DO 210 I=1,N

      W=A(I,M)

      A(I,M)=A(I,J)

  210 A(I,J)=W

  220 J=J-1

      IF(J.GE.1) GO TO 200

  230 ILL=0

      RETURN

  240 ILL=J

      RETURN

  250 ILL=30000

      RETURN

      END


