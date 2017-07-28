!****************************************************************************
!
!  cs.f90 
!
!  FUNCTIONS:
!  Calculate the numerical value of cross section
!  include the file 'diag.f90,diracmatrix.inc,'fgaus.f'
!  2015-10-21 dsigma/dQ^2dw
!****************************************************************************

       module parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(0:3,0:3,4,4) :: c_gsigma          ! Dirac sigma matrix
       complex*16,dimension(0:3,4,4) :: c_ga              ! Dirac Gamma matrix
       complex*16,dimension(4,4) :: c_ga5                     ! Gamma_5 
       real*8,dimension(0:3,0:3) :: g                         ! metrix tensor g,   xunm????
       real*8,dimension(4,4) :: unm,unm0                         ! unit identify matrix 
       real*8,dimension(0:3) :: xk,xkp,xp,xpp,xq,xqf,xpd,xpdc              ! the momentum of initial and final particals   
       complex*16,dimension(4,4) :: c_k,c_kp,c_p,c_pp,c_q,c_qf  
       common/ss/c_ga,c_ga5,g,c_gsigma,unm,unm0,xk,xkp,xp,xpp,xqf,xq, &
           s_kqf,xq0,xpd,xpdc,xenergy,start,finish,beta_cos,xpf0,nprecise, &
           Ddensity,w_sqf0,fca5p33,qfmin,ndiag,nnstep,nucleon,nff,mdelta
           ! ndelta =1 direct Delta , else:direct and crossed
     !************************************
       parameter(xmp=0.938272d0)             ! the mass of proton
       parameter(xmn=0.939)                  ! the mass of neutron 0.939565d0
       parameter(xmn2=xmn**2) 

       parameter(ap=0.942)                    ! these parameters is used to define the form factor F1,F2
       parameter(bp=4.61)
       parameter(xmiup=2.793)
       parameter(xmiun=-1.913)

       parameter(pi = acos(-1.d0))
       parameter(alpha=1.d0/137.d0)          ! this is the e**2/(4*pi)
       
       parameter(cu=cmplx(0.d0,1.d0))        ! the unit of imagic 
 
       parameter(gfermi=1.16637d-5)             ! fermi coupling constant
       parameter(ga=-1.267)                  ! define the formfactor Fa, Diag E(1.26)

       parameter(xmpi=0.138)                 ! the mass of pi mseon, define the formfactor Fp
       
       parameter(delta_s=-0.d0)              ! the parameter of FIT I, GLW93   -0.21, -0.15
       parameter(f1s_0=0.d0)                 ! 0.53
       parameter(f2s_0=-0.d0)                 ! -0.4

       parameter(xmw=80.398)                 ! the mass of w gauge boson
       parameter(xmz=91.1876)                ! the mass of z gauge boson
       
       parameter(gev_bar=3.89379304d-28)     ! unit is cm**2,1(GeV-2)=0.389379mb,1mb=10(-27)(cm2)

       parameter(xmd=1.232)                  ! the mass of Delta 1.232

       parameter(fstar=2.14)                 ! For the width of Delta, Ddelta function

       parameter(sin_theta_w2=0.231)       

       parameter(fpi= 0.0924)                     ! For the Diag E  0.0924


       parameter(x1ma2=1. ) !  1.012**2               ! the necleon axial mass 
       parameter(x2ma2=1.) !  1.012**2                ! the heavier resonances axial mass
       parameter(xma2=0.93**2) !  1.012**2               ! For the factor of Dig CD, xma=1.05, 1.012 
       parameter(xmv2=0.71)  !0.84             ! 0.71    For the Factor of Fig CD, xmv=0.84

     ! the excited resonances mass
       parameter(xmp33=1.232)                           ! the mass of Delta P33(1232)  
       parameter(xmp11=1.44)                            ! the mass of P11(1440)
       parameter(xmd13=1.52)                            ! the mass of D13(1530)
       parameter(xms11=1.535)                           ! the mass of S11(1353)

      ! parameter(fca5p33=1.d0)                          
       parameter(fa0p11=-0.47)                
       parameter(fca5d13=-2.14)!753)!2.14)                           
       parameter(fa0s11=-0.21)                        
      ! parameter(fa0p11=-0.52)                            ! the mass of P11(1440)
      ! parameter(fca5d13=-2.15)                            ! the mass of D13(1530)
      ! parameter(fa0s11=-0.23)                           ! the mass of S11(1353)

     ! for the width.f90
       parameter(xmsigma=0.475)                         ! the mass of Sigma,For Sigma-> Pi Pi in S-wave
       parameter(xmeta=0.548)                           ! the mass of eta
       parameter(xmrho=0.77549)                         ! the mass of rho


       end module
      include 'diag.f90'
    !  include 'diag2.f90'
    !  include 'deltaff_3.f90'
    !  include 'heliamp07.f90'
    !  include 'axialff.f90'
    !  include 'width.f90'

    !  include 'deltamedium.f90'

      !************************************************************* 
      function dxsec(Enu, W, Qsq, nucl, nd, dsigmaneutrino, dsigmaantineutrino)
        use parameter
        implicit real*8 (a,b,d-h,o-y)
        implicit complex*16 (c)
        implicit character(z)
        real*8,dimension(2000) :: x1,x2
        complex*16,dimension(2000) :: cf1,cf2
        character(len=20)::zfile
        common/q2w/tq2,tw
        call cpu_time(start)
        call diracmatrix()
        
        mdelta=0
        tq2=Qsq
        tw=W
        !      print *,'Selcet the diagram: tatal(0),N(1),Delta(2),Pi(3),D13(1520)(4),P11(1440)(5),S11(1535)(6),no resoances(7),Delat+N(8)'
        !      read(*,*),ndiag
        ndiag = nd
        print*, 'ndiag',ndiag
        
        !        print *,'input fca5p33'
        !        read(*,*),fca5p33
        fca5p33=1.
        
        !      print *,'the precise of the integral(5)'
        !      read(*,*),nprecise
        nprecise=5
        
        !      if(ndiag.lt.0.or.ndiag.gt.8) then
        !         print *,'ndiag should be 0 to 6'
        !         stop
        !      else if(ndiag.eq.2.or.ndiag.eq.3) then
        !         nucleon=1
        !        print *,'direct diagram(1) or both(0)'
        !        read(*,*),mdelta
        !      else 
        !        print *,'input the nucleon: Proton (1), Neutron (-1)'
        !        read(*,*),nucleon
        !      endif
        nucleon=nucl
        ! print*,'Please select the nucleon: proton(1) or neutrion(-1)',   nucleon
        
        !      print *,'input the name of file(9)'
        !      read(*,*),zfile
        ! if(nucleon==1) then
        !   open(unit=8,file='dcs_proton.txt')
        ! else
        !   open(unit=8,file='dcs_neutrino.txt')
        ! endif
        ! print*,'The output is the double differential cross section as the function of Q^2 and W' 
        !           'ndiag is ',ndiag, 'nucleon is',nucleon, 'nprecise ',nprecise,'fca5(0)',fca5p33
        
        nnstep=0
        nff=4      ! Form factors of the helicity amplitude
        
        !   if(nucleon==1) then
        !     nfile=ndiag+30
        !   else
        !    nfile=ndiag+20
        !   endif
        
        ! write(8,*) '# neutrion energy(GeV)  neutrino           antineutrino  cross section[10-42 cm2]'
        ! write(8,*) '#ndiag is ',ndiag, 'nucleon is',nucleon, 'nprecise ',nprecise,'fca5(0)',fca5p33
        ! write(8,*) '#nprecise=',nprecise,'namefile=ndiag + 20(neutron) or proton(30)','fca5(0)',fca5p33
        ! print *, '# neutrion energy(GeV)  neutrino           antineutrino  cross section[10-42 cm2]'
        
        
        !  do k=0,18
        xenergy=Enu!0.2+0.1*k
        qfmin=0.14
        s=xmn2+2.*xmn*xenergy 
        
        
        twmin=xmn
        twmax=sqrt(s)
        
        
        !! tq2=0.6
        !        print*,1.
        !        stop
        !        print*, acos(1.000001)
        !        stop
        !c        print*,q2max,wmax,wmin
        !        stop
        sf=alpha/((2.d0*pi)**3*16.d0*xenergy*xmn)*gev_bar*1d42 
        
        ! call DSG20r(twmin,twmax,1,x1,np1)
        
        ! do i=0,200!np1
        !! tw=xmn+(sqrt(s-tq2*(1.+xmn/2./xenergy))-xmn)/200.*i !x1(i)
        
        tq2=(s-tw**2)*(s-xmn2)/s/200.*i
        
        ! print*,i,s-tq2*(1.+xmn/2./xenergy),tw
        ! cycle
        !        tq2min=0.
        !        tq2max=(s-tw**2)*(s-xmn2)/s
        
        !        call DSG20r(tq2min,tq2max, 1,x2,np2)
        !           do j=1,np2
        ! tq2=0. ! x2(j)
        call dcross2(2,cdc)
        !             cf2(j)=cdc
        !           enddo
        
        !           call DRG20c(tq2min,tq2max,1,cf2,ctst2)
        !           cf1(i)=ctst2*2.*tw
        !         enddo
        
        !         call DRG20c(twmin,twmax,1,cf1,ctst1)
        !         ctcr=ctst1*sf
        ctcr=cdc*2*tw*sf
        !  write(8,*) xenergy,real(ctcr),aimag(ctcr),nnstep
        
        ! print *,  xenergy,tq2,tw,real(ctcr),aimag(ctcr),nnstep
        dsigmaneutrino = real(ctcr)
        dsigmaantineutrino = aimag(ctcr)
        call comp_time()
        
        
      END function dxsec
      
!*********************Four integration**************************
      SUBROUTINE dcross2(n,cdc)
        use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(n) :: x
       integer,dimension(n) :: js
       real*8,external :: fs
       complex*16,external :: ckernel
       
       do i=1,n 
        js(i)=nprecise
       enddo
    !    js(2)=8       ! phi_qf  
    !    js(1)=8       ! theta_kp
     !   JS(3)=8       ! xqf0 
     !   JS(4)=8       ! xqf0 

	CALL FGAUS(n,JS,X,FS,ckernel,cs)
        cdc=cs
       RETURN
       END SUBROUTINE


!*****************************************************
       FUNCTION ckernel(n,x)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(n) :: x
       common/q2w/tq2,tw


        xkp0=(xmn2-tq2-tw**2+2.*xmn*xenergy)/2./xmn
        ddd= 1.-tq2/2./xenergy/xkp0
        if(abs( ddd).gt.(1.d0)) then
           print*,'error: the costheta_kp is larger than 1!'
           print*,'costheta_kp=',1.-tq2/2./xenergy/xkp0
           print*,tq2,xenergy,xkp0,tq2,tw
          ! stop
          ckernel=0.
          return
    !    else if(abs(ddd).gt.1) then
    !      ddd=1.d0
        endif

        theta_kp=acos(ddd)
!        print*,sin(theta_kp)
        phi_qf=x(1) 
        theta_qf=x(2)

        xqf0=(xmn*(xenergy-xkp0)-tq2*0.5)/(xmn+xenergy-xkp0-xenergy*cos(theta_qf) &
            +xkp0*(sin(theta_kp)*sin(theta_qf)*cos(phi_qf)+cos(theta_kp)*cos(theta_qf)))

           beta_cos=sin(theta_kp)*sin(theta_qf)*cos(phi_qf) + cos(theta_kp)*cos(theta_qf) 
           s_kqf=xenergy*cos(theta_qf) - xkp0*beta_cos   


         if((xqf0.le.0.d0).or.(xqf0.gt.(xenergy-xkp0))) then  !.or.(abs(beta_cos).gt.1.d0)
            ckernel=0.d0

         else 
       !  print *,'time'
       !      stop
          call momentum(xkp0,theta_kp,xqf0,theta_qf,phi_qf,xmn,0.d0,0.d0,t1,t2)  

          call  amp_num(t1,t2,ctensor_lh)
         ! ndiag=0
         ! call amp_num2(t1,t2,ctensor_lh,c_lh_both_n)      

!c         ckernel=xqf0*sin(theta_qf)*xkp0/abs((xmn+xenergy-xqf0)-s_kqf)&!*4.d0*(2.d0*pi)**5*xenergy*xmn! & 
!c                   *ctensor_lh*gfermi**2/2.d0  &
!c                   *abs(1./4./xenergy/xmn/xkp0) !+tq2/4./xmn2/xenergy/xkp0**2)

         ckernel=xqf0*sin(theta_qf)*xkp0/abs((xmn+xenergy-xkp0)-s_kqf)&!*4.d0*(2.d0*pi)**5*xenergy*xmn! & 
                   *ctensor_lh*gfermi**2/2.d0  &
                   *abs(1./4./xenergy/xmn/xkp0) !+tq2/4./xmn2/xenergy/xkp0**2)
        
         nnstep=nnstep+1
         endif   
         ! print *,ctensor_lh,nnstep
         ! if(nnstep==10) stop
        return
       END FUNCTION
!***************************************************
	SUBROUTINE FS(J,N,X,DN,UP)
        use parameter
        implicit real*8 (a,b,d-h,o-z)
        implicit complex*16 (c)
        real*8,dimension(n) :: x

	IF (J.EQ.1) THEN
          DN=0.d0
          UP=2.d0*pi
	ELSE IF (J.EQ.2) THEN      
	  DN=0.0d0
	  UP=pi
	END IF
	RETURN
	END

  
!****************************
       subroutine comp_time()
       use parameter     
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)

        call cpu_time(finish)
        stime=finish-start
        sh=aint(stime/3600.)
        sm=aint(mod(stime,3600.)/60.)
        ss=mod(stime,60.)
        write(*,200)'The time of calculating is:',sh,'h ',    &
           sm,'m', ss,'s'     
        write(8,200)'#The time of calculating is:',sh,'h ',    &
           sm,'m', ss,'s'  
        print'("Time = ",f10.3," seconds.")',stime
200     Format(a,F5.0,a,F3.0,a,F7.4,a)

        if(ndiag==1.and.nucleon==1) then
         print *,'cs_compton_proton.txt'
        else if(ndiag==1.and.nucleon==-1)then
         print *,'cs_compton_neutron.txt'
        else if(ndiag==2) then
         print *,'cs_delta.txt'
        else if(ndiag==3) then
         print *,'cs_pie.txt'
        else if(ndiag==0.and.nucleon==1) then
         print *,'cs_total(p).txt'
        else if(ndiag==0.and.nucleon==-1) then
         print *,'cs_total(n).txt'    

        endif
        if(neutrino==1) then
           print *,'neutrino'
        else if(neutrino==-1) then
           print *,'antineutrino'
        endif
        Return
        end subroutine comp_time

!******************mult integration **************************************
!***********Definition fo the kinematical in the lab frame *******************
       SUBROUTINE momentum(xkp0,theta_kp,xqf0,theta_qf,phi_qf,xp0,theta_p,phi_p,t1,t2)            
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c) 
       real*8,external :: FMV

       xk(0)=xenergy
       xk(1)=0.d0
       xk(2)=0.d0
       xk(3)=xenergy
   
       spn=sqrt(xp0**2-xmn2)
       xp(0)=xp0
       xp(1)=spn*sin(theta_p)*cos(phi_p)
       xp(2)=spn*sin(theta_p)*sin(phi_p)
       xp(3)=spn*cos(theta_p)
        tpn=theta_p
       if(spn.lt.(-1d-6)) then
          print *,'spn is less than 0', xp0,xpn
          stop
       endif

       xqf(0)=xqf0
       xqf(1)=xqf0*sin(theta_qf)*cos(phi_qf)
       xqf(2)=xqf0*sin(theta_qf)*sin(phi_qf)
       xqf(3)=xqf0*cos(theta_qf)

       xkp(0)=xkp0
       xkp(1)=xkp0*sin(theta_kp)
       xkp(2)=0.d0
       xkp(3)=xkp0*cos(theta_kp)

      
       xpp=xk+xp-xkp-xqf
       xq=xk-xkp
       xpd=xp+xq                      !  P_detla   = p+q = pp+qf
       xpdc=xpp-xq                    !  P_delta_c =pp-q = p-qf

       t1=FMV(xqf,xqf) 
       t2=FMV(xq,xq)
       t3=FMV(xpp,xpp)

       if(xpp(0).lt.xmn) then
          print *, 'xpp(0) is less than xmn',xpp
          print *,xk
          print *,xkp
          print *,xp
          print *,xqf
          stop
       endif
       
       Return 
      END SUBROUTINE 














          



