!***********************************************************
!     2012-06-06                
!
!     This code includes the hardron part and light part of 
!   Diag ABCDE.     
!
!   The code of Diag CD is complicated.
!
!***********************************************************
      SUBROUTINE amp_num(t1,t2,c_lh_both)      
      use parameter
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      complex*16,dimension(0:3,0:3) :: c_tr_l,c_tr_l_anti
      complex*16,dimension(0:3,0:3,4,4) :: ch_ver,ch_ver_t,ch_ver_ab,ch_ver_abt,ch_ver_cd,ch_ver_cdt,ch_ver_e,ch_ver_et,      &
                         ch_ver_d13,ch_ver_d13t,ch_ver_p11,ch_ver_p11t,ch_ver_s11,ch_ver_s11t
      complex*16,dimension(4,4) :: cpm,cppm
     
       call c_gs(xp,c_p)
       call c_gs(xpp,c_pp)
       call c_gs(xk,c_k)
       call c_gs(xkp,c_kp)
       call c_gs(xq,c_q)
       call c_gs(xqf,c_qf)

          cpm =c_p  + xmn*unm
          cppm=c_pp + xmn*unm
     

       If(ndiag==0)  then
        call ver_ab(t1,t2,ch_ver_ab,ch_ver_abt)
        call ver_cd(t1,t2,ch_ver_cd,ch_ver_cdt) 
        call ver_e(ch_ver_e) 
           ch_ver_et=-ch_ver_e
        call ver_j32(-1,t2,ch_ver_d13,ch_ver_d13t)
        call ver_j12( 1,t2,ch_ver_p11,ch_ver_p11t)
        call ver_j12(-1,t2,ch_ver_s11,ch_ver_s11t)

       ch_ver  =ch_ver_cd + ch_ver_e   +ch_ver_ab  + ch_ver_s11 +ch_ver_p11  + ch_ver_d13   
       ch_ver_t=ch_ver_cdt+ ch_ver_et  +ch_ver_abt + ch_ver_s11t+ch_ver_p11t + ch_ver_d13t
       else if(ndiag==1) then
          call ver_ab(t1,t2,ch_ver,ch_ver_t)
       else if(ndiag==2) then
          call ver_cd(t1,t2,ch_ver,ch_ver_t)  
          nparity=1
         ! call ver_j32(nparity,t2,ch_ver,ch_ver_t)
       else if(ndiag==3) then
          call ver_e(ch_ver)
           ch_ver_t=-ch_ver
       else if(ndiag==4) then
          nparity=-1
         call ver_j32(nparity,t2,ch_ver,ch_ver_t)
       else if(ndiag==5) then
          nparity=1
         call ver_j12(nparity,t2,ch_ver,ch_ver_t)
       else if(ndiag==6) then
          nparity=-1
         call ver_j12(nparity,t2,ch_ver,ch_ver_t)
       else if(ndiag==7)  then
          call ver_ab(t1,t2,ch_ver_ab,ch_ver_abt)
          call ver_cd(t1,t2,ch_ver_cd,ch_ver_cdt) 
          call ver_e(ch_ver_e) 
            ch_ver_et=-ch_ver_e
          ch_ver  =ch_ver_cd  +ch_ver_ab   + ch_ver_e !+ ch_ver_s11!! +ch_ver_p11 ! + ch_ver_d13   
          ch_ver_t=ch_ver_cdt  +ch_ver_abt + ch_ver_et
       else if(ndiag==8)  then
          call ver_ab(t1,t2,ch_ver_ab,ch_ver_abt)
          call ver_cd(t1,t2,ch_ver_cd,ch_ver_cdt) 
         ! call ver_e(ch_ver_e) 
         !   ch_ver_et=-ch_ver_e
          ch_ver  =ch_ver_cd  +ch_ver_ab  ! + ch_ver_e !+ ch_ver_s11!! +ch_ver_p11 ! + ch_ver_d13   
          ch_ver_t=ch_ver_cdt  +ch_ver_abt !+ ch_ver_et
       endif

       c_lh=0.d0
       c_lh_anti=0.d0
       call tracelight(c_tr_l,c_tr_l_anti)
       do n_alpha=0,3
         do n_beta=0,3
           
           call c_mult2(cpm,ch_ver_t,cppm,ch_ver, n_alpha,n_beta,c_tr_h)
         
           c_lh=c_lh + c_tr_l(n_alpha,n_beta)*c_tr_h*g(n_alpha,n_alpha)*g(n_beta,n_beta)
           c_lh_anti=c_lh_anti + c_tr_l_anti(n_alpha,n_beta)*c_tr_h*g(n_alpha,n_alpha)*g(n_beta,n_beta)
           
         enddo
       enddo
       
    ! real part is for neutrino, image part is for antineutrino
           if(ncheck==520) then
               if( (abs(aimag(c_lh)).gt.1d-6).or.(abs(aimag(c_lh_anti)).gt.1d-6)  ) then
                  print *, 'The image part of c_lh or c_lh_anti is not zero'
                  print *, c_lh,c_lh_anti
                 stop
               endif
           endif
          
          c_lh_both=real(c_lh)+cu*real(c_lh_anti) 
    
           if(ncheck==520) then
               if( (real(c_lh_both).lt.0.d0).or.(aimag(c_lh_both).lt.0.d0)  ) then
                  print *, 'the LH tensor is less than zero'
                  print *, c_lh_both
                 stop
               endif
           endif
        
       RETURN
       END SUBROUTINE 

!***********************************************************************
      SUBROUTINE amp_num2(t1,t2,c_lh_both_p,c_lh_both_n)      
      use parameter
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      complex*16,dimension(0:3,0:3) :: c_tr_l,c_tr_l_anti 
      complex*16,dimension(0:3,0:3,4,4) :: ch_ver_p,ch_ver_t_p,ch_ver_ab_p,ch_ver_abt_p,    &
                                           ch_ver_n,ch_ver_t_n,ch_ver_ab_n,ch_ver_abt_n,    &
                                           ch_ver_cd,ch_ver_cdt,ch_ver_e,ch_ver_et,         &
                                           ch_ver_p11_p,ch_ver_p11t_p,ch_ver_p11_n,ch_ver_p11t_n, &
                                           ch_ver_d13_p,ch_ver_d13t_p,ch_ver_d13_n,ch_ver_d13t_n, &
                                           ch_ver_s11_p,ch_ver_s11t_p,ch_ver_s11_n,ch_ver_s11t_n
      complex*16,dimension(4,4) :: cpm,cppm,ckm,ckpm
         
       call c_gs(xp,c_p)
       call c_gs(xpp,c_pp)
       call c_gs(xk,c_k)
       call c_gs(xkp,c_kp)
       call c_gs(xq,c_q)
       call c_gs(xqf,c_qf)

          cpm =c_p  + xmn*unm
          cppm=c_pp + xmn*unm
      
       if(ndiag==7)  then
        nucleon=-1
        call ver_ab(t1,t2,ch_ver_ab_n,ch_ver_abt_n)
        nucleon=1
        call ver_ab(t1,t2,ch_ver_ab_p,ch_ver_abt_p)

        call ver_cd(t1,t2,ch_ver_cd,ch_ver_cdt) 
        call ver_e(ch_ver_e) 
             ch_ver_et=-ch_ver_e

        ch_ver_p  =ch_ver_cd + ch_ver_e   +ch_ver_ab_p  
        ch_ver_t_p=ch_ver_cdt+ ch_ver_et  +ch_ver_abt_p 

        ch_ver_n  =ch_ver_cd - ch_ver_e   +ch_ver_ab_n  
        ch_ver_t_n=ch_ver_cdt- ch_ver_et  +ch_ver_abt_n 
   
       else if(ndiag==0) then
        nucleon=-1
        call ver_ab(t1,t2,ch_ver_ab_n,ch_ver_abt_n)
        call ver_j32(-1,t2,ch_ver_d13_n,ch_ver_d13t_n)
        call ver_j12( 1,t2,ch_ver_p11_n,ch_ver_p11t_n)
        call ver_j12(-1,t2,ch_ver_s11_n,ch_ver_s11t_n)
        nucleon=1
        call ver_ab(t1,t2,ch_ver_ab_p,ch_ver_abt_p)
        call ver_j32(-1,t2,ch_ver_d13_p,ch_ver_d13t_p)
        call ver_j12( 1,t2,ch_ver_p11_p,ch_ver_p11t_p)
        call ver_j12(-1,t2,ch_ver_s11_p,ch_ver_s11t_p)

        call ver_cd(t1,t2,ch_ver_cd,ch_ver_cdt) 
        call ver_e(ch_ver_e) 
             ch_ver_et=-ch_ver_e

       ch_ver_p  =ch_ver_cd + ch_ver_e   +ch_ver_ab_p  +  ch_ver_d13_p  +ch_ver_p11_p  + ch_ver_s11_p
       ch_ver_t_p=ch_ver_cdt+ ch_ver_et  +ch_ver_abt_p +  ch_ver_d13t_p +ch_ver_p11t_p + ch_ver_s11t_p

       ch_ver_n  =ch_ver_cd - ch_ver_e   +ch_ver_ab_n  +  ch_ver_d13_n  +ch_ver_p11_n  + ch_ver_s11_n  
       ch_ver_t_n=ch_ver_cdt- ch_ver_et  +ch_ver_abt_n +  ch_ver_d13t_n +ch_ver_p11t_n + ch_ver_s11t_n 

       else
        print *,'Ndiag should be 0 or 7'
        stop
       endif


       c_lh_p=0.d0
       c_lh_n=0.d0
       c_lh_anti_p=0.d0
       c_lh_anti_n=0.d0
       call tracelight(c_tr_l,c_tr_l_anti)

       do n_alpha=0,3
         do n_beta=0,3
                 
           call c_mult2(cpm,ch_ver_t_p,cppm,ch_ver_p, n_alpha,n_beta,c_trh_p)
           call c_mult2(cpm,ch_ver_t_n,cppm,ch_ver_n, n_alpha,n_beta,c_trh_n)
   

           c_lh_p=c_lh_p + c_tr_l(n_alpha,n_beta)*c_trh_p*g(n_alpha,n_alpha)*g(n_beta,n_beta)
           c_lh_n=c_lh_n + c_tr_l(n_alpha,n_beta)*c_trh_n*g(n_alpha,n_alpha)*g(n_beta,n_beta)

           c_lh_anti_p=c_lh_anti_p + c_tr_l_anti(n_alpha,n_beta)*c_trh_p*g(n_alpha,n_alpha)*g(n_beta,n_beta)
           c_lh_anti_n=c_lh_anti_n + c_tr_l_anti(n_alpha,n_beta)*c_trh_n*g(n_alpha,n_alpha)*g(n_beta,n_beta)

          enddo
       enddo
   
           if(ncheck==520) then
               if( (abs(aimag(c_lh_p)).gt.1d-6).or.(abs(aimag(c_lh_anti_p)).gt.1d-6).or.  &
                      (abs(aimag(c_lh_n)).gt.1d-6).or.(abs(aimag(c_lh_anti_n)).gt.1d-6)) then
                  print *, 'The image part of c_lh or c_lh_anti is not zero'
                  print *, 'c_lh_p,c_lh_anti_p,c_lh_n,c_lh_anti_n'
                  print *, c_lh_p,c_lh_anti_p,c_lh_n,c_lh_anti_n
                 stop
               endif
           endif

          c_lh_both_p=real(c_lh_p)+cu*real(c_lh_anti_p) 
          c_lh_both_n=real(c_lh_n)+cu*real(c_lh_anti_n) 

           if(ncheck==520) then
               if( (real(c_lh_both_p).lt.0.d0).or.(aimag(c_lh_both_p).lt.0.d0).or.  &
                    (real(c_lh_both_n).lt.0.d0).or.(aimag(c_lh_both_n).lt.0.d0)    ) then
                  print *, 'the LH tensor is less than zero'
                  print *, 'c_lh_both_p,c_lh_both_n'
                  print *, c_lh_both_p,c_lh_both_n
                 stop
               endif
           endif
 
       RETURN
       END SUBROUTINE 

!**********************************************************
!
!==============Light part==================================
!******************************************************
       SUBROUTINE tracelight(c_tr_l,c_tr_l_anti)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c) 


       complex*16,dimension(0:3,0:3) :: c_tr_l,c_tr_l_anti 
        
       do n_alpha=0,3
         do n_beta=0,3
           c_tr_l(n_alpha,n_beta)=0.d0
              stl=0.d0
           do n_sigma=0,3
             do n_delta=0,3
              stl=stl - Lvten(n_alpha,n_beta,n_sigma,n_delta)     &
                    *xk(n_sigma)*xkp(n_delta)*g(n_sigma,n_sigma)*g(n_delta,n_delta)
              enddo
           enddo
         !  if(neutrino==1) then
           c_tr_l(n_alpha,n_beta)=8.d0*(xk(n_alpha)*xkp(n_beta)+xk(n_beta)*xkp(n_alpha)  &
              -FMV(xk,xkp)*g(n_alpha,n_beta)+cu*stl )
         !  else if(neutrino==-1) then
           c_tr_l_anti(n_alpha,n_beta)=8.d0*(xk(n_alpha)*xkp(n_beta)+xk(n_beta)*xkp(n_alpha)  &
              -FMV(xk,xkp)*g(n_alpha,n_beta)-cu*stl )
         !  else 
         !     print *,'The neutrino only can be 1 and -1'
         !     stop
         !  endif
           enddo
         enddo 

        RETURN
        END SUBROUTINE 

!==============Hardon part Diag A and B==================================
!****************The form factor of EM***********************************
      Function smomentum(sp)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: sp
      
       smomentum=0.d0       
       do i=0,3
         smomentum=smomentum + sp(i)*sp(i)*g(i,i)     
       enddo
  
      Return
      END FUNCTION 
       SUBROUTINE ff_em(t,f1,f2)     
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)


       gd=(1.d0-t/xmv2)**(-2)
       gmp=xmiup*gd
       gmn=xmiun*gd

       tao=t/(4.d0*xmn2)
       gep=gd
       gen=xmiun*(ap*tao/(1.d0-bp*tao))*gd       

       if(nucleon==1) then
          ge=gep
          gm=gmp
         else
          ge=gen
          gm=gmn
       endif

       f1=(ge-tao*gm)/(1.d0-tao)
       f2=(gm-ge)/(1.d0-tao)

       Return 
       END SUBROUTINE
!******************The form factor of NC*********************************
       SUBROUTINE ff_nc(t,fv1,fv2,fva)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
     
       tao=t/(4.d0*xmn2)    

       gd=1.d0/(1.d0-t/xmv2)**2
       gmp=xmiup*gd
       gmn=xmiun*gd

       gep=gd
       gen=xmiun*(ap*tao/(1.d0-bp*tao))*gd       

       f1_p=(gep-tao*gmp)/(1.d0-tao)
       f2_p=(gmp-gep)/(1.d0-tao)

       f1_n=(gen-tao*gmn)/(1.d0-tao)
       f2_n=(gmn-gen)/(1.d0-tao)

       fa=ga/(1.d0-t/x1ma2)**2

       f1_s=-f1s_0*t /(1.d0-tao)/(1.d0-t/xmv2)**2
       f2_s=f2s_0/(1.d0-tao)/(1.d0-t/xmv2)**2
       fa_s=delta_s/(1.d0-t/x1ma2)**2

       stheta=1.d0-4.d0*sin_theta_w2  
         
       if(nucleon==1) then
           fv1=(stheta*f1_p - f1_n - f1_s)/2.d0
           fv2=(stheta*f2_p - f2_n - f2_s)/2.d0  
           fva=(fa + fa_s)/2.d0         
         else if (nucleon==-1) then       
           fv1=(stheta*f1_n - f1_p - f1_s)/2.d0
           fv2=(stheta*f2_n - f2_p - f2_s)/2.d0 
           fva=(-fa + fa_s)/2.d0 
       endif

       Return 
       END SUBROUTINE
!**************************************************************
       SUBROUTINE ver_ab(t1,t2,ch_vertex1,ch_vertex2)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(0:3,0:3,4,4) :: ch_vertex1,ch_vertex2
       complex*16,dimension(0:3,4,4) :: cA_nc,cA_em,cA_nc_con,cA_em_con !cl_vertex
       complex*16,dimension(4,4) :: cA_nc2,cA_em2,cA_nc_con2,cA_em_con2
       complex*16,dimension(4,4) :: cpqm,cppqm,c_ga_n1,c_ga_n2,c_ga_con1,c_ga_con2,c_ga0,c_p0,c_q0

       

        call ff_em(t1,f1,f2)
        call ff_nc(t2,ft1,ft2,fta)
               
 
        do n1=0,3
          do i=1,4
           do j=1,4
              c_sigma_q = 0.d0
              c_sigma_qf = 0.d0
             do n2=0,3
              c_sigma_q = c_sigma_q + c_gsigma(n1,n2,i,j)*xq(n2)*g(n2,n2) 
              c_sigma_qf = c_sigma_qf + c_gsigma(n1,n2,i,j)*xqf(n2)*g(n2,n2) 
             enddo      

             c_ga_5 = 0.d0
             do k=1,4
               c_ga_5 =c_ga_5 + c_ga(n1,i,k)*c_ga5(k,j)
             enddo
   
    ! the Fig A and Fig B
           cA_nc(n1,i,j)=ft1*c_ga(n1,i,j)                &          
                       + ft2*cu*c_sigma_q/(2.d0*xmn)     &
                       + c_ga_5*fta                                  ! for the vertex of hadron

           cA_em(n1,i,j)= f1*c_ga(n1,i,j) - f2*cu*c_sigma_qf/(2.d0*xmn)
           
           cA_nc_con(n1,i,j)=ft1*c_ga(n1,i,j)           &            ! for the conj vertex
                       - ft2*cu*c_sigma_q/(2.d0*xmn)    &           
                       + c_ga_5*fta                                 

           cA_em_con(n1,i,j)= f1*c_ga(n1,i,j)   + f2*cu*c_sigma_qf/(2.d0*xmn)        

           enddo
         enddo
      enddo  


      cpqm =c_p + c_q + xmn*unm
      cppqm=c_pp - c_q + xmn*unm

      call d_prog(xp+xq,dpq)
      call d_prog(xpp-xq,dppq)

     ! if(nucleon==1.and.w_sqf0.lt.0.1) then
     !       dpq=0.
     !       dppq=0.
     ! endif
      
      do n_miu=0,3
        do n_alpha=0,3
          call matrix_32(cA_em,n_miu,cA_em2)
          call matrix_32(cA_nc,n_alpha,cA_nc2) 

          call matrix_32(cA_em_con,n_miu,cA_em_con2)
          call matrix_32(cA_nc_con,n_alpha,cA_nc_con2) 

          call mult_matrix(cA_em2,cpqm,cA_nc2,c_ga_n1)  
          call mult_matrix(cA_nc2,cppqm,cA_em2,c_ga_n2)

          call mult_matrix(cA_nc_con2,cpqm,cA_em_con2,c_ga_con1) 
          call mult_matrix(cA_em_con2,cppqm,cA_nc_con2,c_ga_con2)

                 
          do i=1,4
            do j=1,4
              ch_vertex1(n_miu,n_alpha,i,j)=c_ga_n1(i,j)*dpq + c_ga_n2(i,j)*dppq
              ch_vertex2(n_miu,n_alpha,i,j)=c_ga_con1(i,j)*dpq +c_ga_con2(i,j)*dppq
            enddo
          enddo 

        enddo
     enddo

       RETURN
       END SUBROUTINE
!********* P -> DiracSlash(P)*****************************
       SUBROUTINE c_gs(p,cp)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c) 
       real*8,dimension(0:3) :: p
       complex*16,dimension(4,4) :: cp 

       do i=1,4
         do j=1,4
           cp(i,j)=0.d0
           do n=0,3
            cp(i,j)=cp(i,j)+p(n)*c_ga(n,i,j)*g(n,n)
           enddo
         enddo
       enddo
        
        return
       END SUBROUTINE

!********Trace of hadron tensor***************************************
      SUBROUTINE c_mult2(c_a,c_b,c_c,c_d,n_alpha,n_beta,ch2) 
      use parameter
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      parameter(n=4)                                   ! n is the dimension of the matrixs
      complex*16,dimension(n,n) :: c_a,c_c,c_m1,c_m2,c_m 
      complex*16,dimension(0:3,0:3,n,n) :: c_b,c_d
      
      ch=0.d0

      do n_miu=0,3
      do n_niu=0,3

      do i=1,n
        do j=1,n
            c_m1(i,j)=0.d0
            c_m2(i,j)=0.d0
          do k=1,n
            c_m1(i,j)=c_m1(i,j) + c_a(i,k)*c_b(n_miu,n_alpha,k,j)
            c_m2(i,j)=c_m2(i,j) + c_c(i,k)*c_d(n_niu,n_beta,k,j)
          enddo
        enddo
      enddo

      c_m=matmul(c_m1,c_m2)
      
      ctr=0.d0
      do i=1,n
        ctr=ctr + c_m(i,i)
      enddo

      ch=ch+ctr*g(n_miu,n_niu)
      enddo
      enddo
 
      ch2=-ch/2.d0

      RETURN      
      END SUBROUTINE
!********Trace of four lepton tensor***************************************
      SUBROUTINE c_mult1(c_a,c_b,c_c,c_d,n_miu,n_niu,ctr) 

      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      parameter(n=4)                                   ! n is the dimension of the matrixs
      complex*16,dimension(n,n) :: c_a,c_c,c_m1,c_m2,c_m 
      complex*16,dimension(0:3,n,n) :: c_b,c_d
      
      do i=1,n
        do j=1,n
            c_m1(i,j)=0.d0
            c_m2(i,j)=0.d0
          do k=1,n
            c_m1(i,j)=c_m1(i,j) + c_a(i,k)*c_b(n_miu,k,j)
            c_m2(i,j)=c_m2(i,j) + c_c(i,k)*c_d(n_niu,k,j)
          enddo
        enddo
      enddo

      c_m=matmul(c_m1,c_m2)
      
      ctr=0.d0
      do i=1,n
        ctr=ctr + c_m(i,i)
      enddo

      RETURN      
      END SUBROUTINE


!********Trace of four matrix multiply***************************************
      SUBROUTINE mult_matrix(c_a,c_b,c_c,c_m) 
      implicit real*8 (a,b,d-h,o-z)
      implicit complex*16 (c)
      parameter(n=4)                                   ! n is the dimension of the matrixs
      complex*16,dimension(n,n) :: c_a,c_b,c_c,c_m1,c_m 
         c_m1=matmul(c_a,c_b)
         c_m=matmul(c_m1,c_c)
        
      RETURN      
      END SUBROUTINE

!*********Change the dimension of matrix from 3 to 2****************
       SUBROUTINE matrix_32(cmatrix_3,n,cmatrix_2)
       
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(0:3,4,4) :: cmatrix_3 
       complex*16,dimension(4,4) :: cmatrix_2 

       do i=1,4
         do j=1,4
            cmatrix_2(i,j)=cmatrix_3(n,i,j)
         enddo
       enddo
        
        return
       END SUBROUTINE   
!*********Change the dimension of matrix from 3 to 2****************
       SUBROUTINE matrix_42(cmatrix_4,n1,n2,cmatrix_2)
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(0:3,0:3,4,4) :: cmatrix_4 
       complex*16,dimension(4,4) :: cmatrix_2 

       do i=1,4
         do j=1,4
            cmatrix_2(i,j)=cmatrix_4(n1,n2,i,j)
         enddo
       enddo
        
        return
       END SUBROUTINE  

!****************D(Q)=1/(Q**2-M**2)**********
       SUBROUTINE d_prog(am,dd)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: am
       real*8,external :: smomentum
       
       temp=smomentum(am)
       dd=1.d0/(temp - xmn2)

       RETURN
       END SUBROUTINE
!**************************************************************
!
!=============Hardon part Diag C and D=========================
!
!**************************************************************
  SUBROUTINE Factor_C(t1,t2,fcv,fcvt,fcat)   ! the factor with ~t is for NC,  without ~t is for EM
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(3:6) :: fcat
       real*8,dimension(3:5) :: fcv,fcvt                        ! simplify, some factors are zero

       if(nff==1) then
   
       sma2=1.05d0**2
       sfca50=1.2
       sfcv30=1.95
       
       fcv(3)=sfcv30   !/(1.d0- t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))    ! fcv30=1.95   xmv=0.84GeV
       fcv(4)=-xmn/xmd*fcv(3)
       fcv(5)=0.
        
       fcat(3)=0.d0
       fcat(5)=sfca50/(1.d0- t2/sma2)**2/(1.d0- t2/(3.*sma2))    ! xma=1.05 GeV, fca50=1.2
       fcat(4)=-fcat(5)/4.d0
       fcat(6)=fcat(5)*xmn2/(xmpi**2-t2)

       fcvt(3)=(1.d0-2.d0*sin_theta_w2)*sfcv30/(1.d0-t2/xmv2)**2/(1.d0-t2/(4.d0*xmv2))
       fcvt(4)=-xmn/xmd*fcvt(3)
       fcvt(5)=0.d0
        
        else if(nff==2) then
     ! Hills form factor   
         sma2=1.05d0**2
         sfca50=1.2
       
       fcv(3)=2.13  !/(1.d0-t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))
       fcv(4)=-1.51 !/(1.d0-t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))
       fcv(5)=0.48  !/(1.d0-t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))
        
       fcat(3)=0.d0
       fcat(5)=sfca50/(1.d0- t2/sma2)**2/(1.d0- t2/(3.*sma2))  
       fcat(4)=-fcat(5)/4.d0
       fcat(6)=fcat(5)*xmn2/(xmpi**2-t2)

       fcvt(3)=2.13d0/(1.d0-t2/xmv2)**2/(1.d0-t2/(4.d0*xmv2))*(1.d0-2.d0*sin_theta_w2)
       fcvt(4)=-1.51d0/(1.d0-t2/xmv2)**2/(1.d0-t2/(4.d0*xmv2))*(1.d0-2.d0*sin_theta_w2)
       fcvt(5)=0.48d0/(1.d0-t2/xmv2)**2/(1.d0-t2/(0.776d0*xmv2))*(1.d0-2.d0*sin_theta_w2)

        else if(nff==3) then
       fca50=fca5p33
       fcv(3)=2.13  !/(1.d0-t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))
       fcv(4)=-1.51 !/(1.d0-t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))
       fcv(5)=0.48  !/(1.d0-t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))
        
       fcat(3)=0.d0
       fcat(5)=fca50/(1.d0- t2/xma2)**2  
       fcat(4)=-fcat(5)/4.d0
       fcat(6)=fcat(5)*xmn2/(xmpi**2-t2)

       fcvt(3)=2.13d0/(1.d0-t2/xmv2)**2/(1.d0-t2/(4.d0*xmv2))*(1.d0-2.d0*sin_theta_w2)
       fcvt(4)=-1.51d0/(1.d0-t2/xmv2)**2/(1.d0-t2/(4.d0*xmv2))*(1.d0-2.d0*sin_theta_w2)
       fcvt(5)=0.48d0/(1.d0-t2/xmv2)**2/(1.d0-t2/(0.776d0*xmv2))*(1.d0-2.d0*sin_theta_w2)


        else if(nff==4) then
          fca50=fca5p33
          q2=0.d0
          call deltaff(q2,sc3v,sc4v,sc5v,sc3a,sc4a,sc5a,sc6a)
          fcv(3)=sc3v
          fcv(4)=sc4v
          fcv(5)=sc5v
          call deltaff(t2,s3v,s4v,s5v,s3a,s4a,s5a,s6a)
          fcvt(3)=(1.d0-2.d0*sin_theta_w2)*s3v
          fcvt(4)=(1.d0-2.d0*sin_theta_w2)*s4v
          fcvt(5)=(1.d0-2.d0*sin_theta_w2)*s5v


         fcat(3)=0.d0
         fcat(5)=fca50/(1.d0- t2/xma2)**2  
         fcat(4)=-fcat(5)/4.d0
         fcat(6)=fcat(5)*xmn2/(xmpi**2-t2)
     else if(nff==5) then
     ! Hills form factor   
         smv2=0.8**2           ! 
         sma2=1.d0
         sfca50=1.2
         sfcv30=2.0 
       
       fcv(3)=sfcv30   !/(1.d0- t1/xmv**2)**2/(1.d0-t1/(4.d0*xmv**2))    ! fcv30=1.95   xmv=0.84GeV
       fcv(4)=-xmn/xmd*fcv(3)
       fcv(5)=0.d0
        
       fcat(3)=0.d0
       fcat(5)=sfca50/(1.d0- t2/sma2)**2 !/(1.d0- t2/(3.*xma2))  
       fcat(4)=-fcat(5)/4.d0
       fcat(6)=fcat(5)*xmn2/(xmpi**2-t2)

       fcvt(3)=(1.d0-2.d0*sin_theta_w2)*sfcv30/(1.d0-t2/xmv2)**2!/(1.d0-t2/(4.d0*xmv**2))
       fcvt(4)=-xmn/xmd*fcvt(3)
       fcvt(5)=0.d0

        else 
          print *, 'nff should be 1 to 4'
          stop
        endif

         RETURN
       END SUBROUTINE 
!**************************************************************
       SUBROUTINE ver_cd(t1,t2,ch_ver_cd,ch_ver_cdt)          ! Gamma_delta for the Fia C and Fig D
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(0:3,0:3,4,4) :: ch_ver_cd,ch_ver_cdt,cAC_emt,cAC_nc,        &
                               cAD_em,cAD_nct,cAC_em,cAD_nc
       complex*16,dimension(4,4) :: clpq,clppq,cg_c,cg_d,cgtem_c,cgtem_d,               &
                               cacemt,cacnc,cadnct,cadem

        cdeltapq=  cDdelta(xpd)         ! xpd=xp+xq
        cdeltappq= cDdelta(xpdc)       ! xpdc=xpp-xq
   ! in_medium affection of Delta 
        
       !call medium(xpd,densi,cdeltapq)
       !call medium(xpdc,densi,cdeltappq)


       call AEM(t1,t2,xpp,xqf,cAC_em)
       call conj(cac_em,cac_emt)
       call ANC(t1,t2,xp,xq,cAC_nc)

       call ANC(t1,t2,xpp,-xq,cAD_nc)          
       call conj(cad_nc,cad_nct)
       call AEM(t1,t2,xp,-xqf,cAD_em)     

        do nmiu=0,3
          do nalpha=0,3

           cgtem_c=unm0
           cgtem_d=unm0


           do ndelta=0,3
             do nsigma=0,3

         ! For Fig C    
              call matrix_42(cAC_emt,ndelta,nmiu,cacemt)
              call matrix_42(cAC_nc, nsigma,nalpha,cacnc)

              call Lamda(ndelta,nsigma,xpd,clpq)

              call Mult_matrix(cacemt,clpq,cacnc,cg_c)
               
              cgtem_c =cgtem_c  + cg_c *g(ndelta,ndelta)*g(nsigma,nsigma)

         ! For Fig D
              call matrix_42(cAD_nct,ndelta,nalpha,cadnct)
              call matrix_42(cAD_em,nsigma,nmiu,cadem)

              call Lamda(ndelta,nsigma,xpdc,clppq)

              call Mult_matrix(cadnct,clppq,cadem,cg_d)     

              cgtem_d =cgtem_d  + cg_d*g(ndelta,ndelta)*g(nsigma,nsigma)
           ENDDO
          ENDDO   
             if(ndiag==2.and.mDelta==1) cdeltappq=0.
          DO i=1,4
           DO j=1,4
             ch_ver_cd(nmiu,nalpha,i,j) =cgtem_c(i,j) *cdeltapq +cgtem_d(i,j)*cdeltappq 
           enddo
          enddo
         ENDDO
        ENDDO

      call conj(ch_ver_cd,ch_ver_cdt)
       Return
      END SUBROUTINE
!*******************************************************
       SUBROUTINE ANC(t1,t2,sp,sq,canc)       
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: sp,sq,spd
       real*8,dimension(3:6) :: fcat
       real*8,dimension(3:5) :: fcv,fcvt
       complex*16,dimension(4,4) :: cq,cnc1,cnc2
       complex*16,dimension(0:3,0:3,4,4) :: CANC
       real*8,external:: FMV
     
       spd=sp+sq
       call c_gs(sq,cq)
 !      t1=0.d0
 !      t2=FMV(sq,sq)
       call Factor_c(t1,t2,fcv,fcvt,fcat)
    
       do n1=0,3
         do n2=0,3
          
        cs2=fcvt(4)/xmn2*(g(n1,n2)*FMV(sq,spd)-sq(n1)*spd(n2))      &
               +fcvt(5)/xmn2*(g(n1,n2)*FMV(sq,sp)-sq(n1)*sp(n2)) 

       ctemp1=fcat(4)/xmn2*(g(n1,n2)*FMV(sq,spd) - sq(n1)*spd(n2))             &
                      +fcat(5)*g(n1,n2)                                        &
                      +fcat(6)/xmn2*sq(n1)*sq(n2)                           
                         
       DO i=1,4
         DO j=1,4     

          do j1=1,4
             cs1=fcvt(3)/xmn*(g(n1,n2)*cq(i,j1)-sq(n1)*c_ga(n2,i,j1)) 
             cnc1(i,j1) =cs1+cs2*unm(i,j1) 
          enddo

            cnc2(i,j) =0.d0
          do k=1,4
            cnc2(i,j) =cnc2(i,j) +cnc1(i,k) *c_ga5(k,j)
          enddo  

            CANC(n1,n2,i,j) =cnc2(i,j) +ctemp1*unm(i,j)
        ENDDO 
      ENDDO

      enddo
      enddo
      
      RETURN
      END SUBROUTINE
!***************************************************************************************
       SUBROUTINE AEM(t1,t2,sp,sq,caem)       
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: sp,sq,spd
       real*8,dimension(3:6) :: fcat
       real*8,dimension(3:5) :: fcv,fcvt
       complex*16,dimension(4,4) :: cq,cem
       complex*16,dimension(0:3,0:3,4,4) :: CAEM
       real*8,external :: FMV

       spd=sp+sq
       call c_gs(sq,cq)
!       t1=0.d0
!       t2=FMV(xq,xq)
       call Factor_c(t1,t2,fcv,fcvt,fcat)                

       do n1=0,3
        do n2=0,3
           cs2=fcv(4)/xmn2*(g(n1,n2)*FMV(sq,spd)-sq(n1)*spd(n2))     &
               +fcv(5)/xmn2*(g(n1,n2)*FMV(sq,sp)-sq(n1)*sp(n2)) 

         DO i=1,4
           DO j=1,4    

             do j2=1,4
              cs1=fcv(3)/xmn*(g(n1,n2)*cq(i,j2)-sq(n1)*c_ga(n2,i,j2)) 
              cem(i,j2)=cs1+cs2*unm(i,j2)    
            enddo
            CAEM(n1,n2,i,j)=0.d0
            
           do k2=1,4
              CAEM(n1,n2,i,j)=CAEM(n1,n2,i,j)+cem(i,k2)*c_ga5(k2,j)
           enddo  
        enddo
       ENDDO 
     
       enddo
       enddo

       Return
       END SUBROUTINE
!*****************************************************************
       subroutine conj(cm,conjm)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(4,4) :: cg0,cm0,cm1,cm2
       complex*16,dimension(0:3,0:3,4,4) :: cm,conjm

       call matrix_32(c_ga,0,cg0)

       DO n1=0,3
        Do n2=0,3

         call matrix_42(cm,n1,n2,cm0)

         DO i=1,4
          DO j=1,4
             cm1(j,i) =conjg(cm0(i,j))
          enddo
        enddo
       
        cm2=matmul(matmul(cg0,cm1),cg0)
   
         DO i=1,4
          DO j=1,4
             conjm(n1,n2,i,j) =cm2(i,j)
          enddo
        enddo
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE



!**************************************************************
       SUBROUTINE Lamda(ns,nd,ppd,clam)       !  ns->n_sigma,nd->n_delta
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(4,4) :: cppd,clam,cpmd,csbar
       real*8,dimension(0:3) :: ppd
       
       call c_gs(ppd,cppd)
       cpmd=cppd+xmd*unm
       do i=1,4
         do j=1,4 
             csdij=0.d0
           do k=1,4
             csdij=csdij+c_ga(ns,i,k)*c_ga(nd,k,j)
           enddo   
     
        csbar(i,j)=g(ns,nd)*unm(i,j)-2.d0/3.d0*ppd(ns)*ppd(nd)/xmd**2*unm(i,j)     &
                 +(ppd(ns)*c_ga(nd,i,j)-ppd(nd)*c_ga(ns,i,j))/(3.d0*xmd)  &
                 -csdij/3.d0 
          enddo
       enddo

       do i=1,4
         do j=1,4
          clam(i,j)=0.d0
           do k=1,4
             clam(i,j)=clam(i,j)-cpmd(i,k)*csbar(k,j)
           enddo
          enddo
       enddo

       Return
      END SUBROUTINE


!***************D_delta FUNCTION****************************************
       Function cDdelta(spp)       
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: spp
       real*8,external :: Flam

      ! spp2=smomentum(spp)
       spp2=spp(0)**2-spp(1)**2-spp(2)**2-spp(3)**2
       if((spp2-(xmn+xmpi)**2).GT.0) then
       slam=Flam(spp2,xmpi**2,xmn2)
       wd=1.d0/(6.d0*pi)*(fstar/xmpi)**2 *xmn/spp2**2  *                                 &
           (sqrt(slam)/(2.d0))**3  !*Ftheta(spsqr-(xmn+xmpi)**2)            !  csw-xmn-xmpi
        else 
           wd=0.d0
        endif
      cDdelta=1.d0/(spp2-xmd**2+cu*wd*xmd)                    !1.d0/(csw+xmd)/(csw-xmd+cu*cwd/2.d0)
   !***************++hills
    !    IF ((spp2-(xmn+xmpi)**2).GT.0) then
    !    wd=0.12*( sqrt(Flam(spp2,xmpi**2,xmn2)) *xmd/(sqrt(spp2)*sqrt(Flam(xmd**2,xmpi**2,xmn2))))**3
    !   else 
    !   wd=0.
    !   endif

     !  cDdelta=1.d0/(spp2-xmd**2+cu*wd*xmd)
       Return
      END FUNCTION



!******************************************************
       Function Flam(sx,sy,sz)       
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
             
       Flam=sx**2 + sy**2 + sz**2 -2.d0*(sx*sy + sy*sz + sx*sz)
             
       Return
       END FUNCTION
!*****************************************************
       Function FMV(xp1,xp2)       
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: xp1,xp2
       
       FMV=0.d0
       do n=0,3
         FMV=FMV+xp1(n)*xp2(n)*g(n,n)
       enddo
           
       Return
      END FUNCTION


!=============Hardon part Diag E=========================
!
!**************************************************************
!*************************************************************
!*************************************************************
        FUNCTION Lvten(n0,n1,n2,n3)
        implicit real*8 (a,b,d-h,o-z)
        implicit complex*16 (c)
        integer,dimension(0:3) :: nvector
   
        nvector(0)=n0
        nvector(1)=n1
        nvector(2)=n2
        nvector(3)=n3

       do i=0,2
          do j=i+1,3
           if(nvector(i)==nvector(j)) then
             lvten=0
             goto 20
           endif
         enddo
       enddo    

        num=0
        do i=0,2
          do j=i+1,3       
            if (nvector(i).GT.nvector(j)) then
              ntem=nvector(j)
              nvector(j)=nvector(i)
              nvector(i)=ntem
              num=num+1
            endif
          enddo
        enddo

        if(mod(num,2)==0) then
           Lvten=1
        else 
           Lvten=-1
        endif
      

20        Return
        END FUNCTION 

!***************************************************************      
       SUBROUTINE ver_e(ch_ver_e) 
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(0:3,0:3,4,4) :: ch_ver_e
       complex*16,dimension(4,4) :: csp5,csppp
       real*8,dimension(0:3,0:3) :: st
       integer,external :: Lvten

       call deltapi(xpp-xp,sdel)
       csfactor=-cu*nucleon*(-ga)/(8.*pi**2*fpi**2)*(0.5-2.*sin_theta_w2)*sdel !4.*sin_theta_w2   (1.-0.23122)

       csppp=c_pp-c_p

       csp5=matmul(csppp,c_ga5)

       do nmiu=0,3
         do nalpha=0,3
           st(nmiu,nalpha)=0.d0
           do nsigma=0,3
             do ndelta=0,3
               st(nmiu,nalpha)=st(nmiu,nalpha) + Lvten(ndelta,nsigma,nmiu,nalpha)     &
                    *xqf(ndelta)*xq(nsigma)
              enddo
           enddo
            
           do i=1,4
             do j=1,4
               ch_ver_e(nmiu,nalpha,i,j)=st(nmiu,nalpha)*csp5(i,j)
             enddo
           enddo

         enddo
       enddo

       ch_ver_e=ch_ver_e*csfactor
    !   ch_ver_et=-ch_ver_e  

       Return
       END SUBROUTINE
!*************************************************
       SUBROUTINE Deltapi(spp,sdel)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: spp
        
 
       temp=smomentum(spp)

       sdel=1.d0/(temp - xmpi**2)
      
       RETURN
       END SUBROUTINE 
           
!************************************************************
!    This subroutine is for the excited resonances contribution
!    include P11(1440) D13(1520)  S11(1535)
!
!
!************************************************************

!=============Hardon part P11(1440) and S11(1535)=========================
       SUBROUTINE ver_j12(nparity,t1,ch_verj12,ch_verj12_t)    ! nparity=1 for P11(1440), nparity=-1 for S11(1535)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(2,6) :: fem,fvem,fvnc
       real*8,dimension(2) :: fanc
       complex*16,dimension(0:3,0:3,4,4) :: ch_verj12,ch_verj12_t
       complex*16,dimension(0:3,4,4) :: cA_nc,cA_em
       complex*16,dimension(4,4) :: cA_nc2,cA_em2,cpqm,cppqm,c_ga_n1,c_ga_n2


       cpqm =c_p + c_q + xmn*unm
       cppqm=c_pp- c_q + xmn*unm

       if(nparity==1) then
         call ffP11(t1,fem,fvem,fvnc,fanc)
         nexcit=2
       else if(nparity==-1) then
         call ffS11(t1,fem,fvem,fvnc,fanc)
         nexcit=4
       else 
         print*,'nparity should be 1 or -1. Here is Form Factor of J=1/2'
         stop
       endif 

       ii=-(nucleon-3)/2       ! For nuetron, nucleon=-1, i=2, for proton, nucleon=1, i=1
       f1=fem(ii,1)
       f2=fem(ii,2)
       fa=0.d0
       ft1=fvnc(ii,1)
       ft2=fvnc(ii,2)
       fta=fanc(ii)  

       call vertex12(nparity,f1,f2,fa,-xqf,cA_em)
       call vertex12(nparity,ft1,ft2,fta,xq,cA_nc)
       
       call propagator(nexcit,xpd,cdpq)
       call propagator(nexcit,xpdc,cdppq)
        
       do n_miu=0,3
         do n_alpha=0,3
          call matrix_32(cA_em,n_miu,cA_em2)
          call matrix_32(cA_nc,n_alpha,cA_nc2) 

          call mult_matrix(cA_em2,cpqm,cA_nc2,c_ga_n1)  

          call mult_matrix(cA_nc2,cppqm,cA_em2,c_ga_n2)
                 
          do i=1,4
            do j=1,4
              ch_verj12(n_miu,n_alpha,i,j)=c_ga_n1(i,j)*cdpq + c_ga_n2(i,j)*cdppq
            enddo
          enddo 

        enddo
       enddo
         
        call conj(ch_verj12,ch_verj12_t)

       RETURN
       END SUBROUTINE

!**************V-A(1/2)  vertex12***************************
       SUBROUTINE vertex12(nparity,f1,f2,fa,sq,ver12)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: sq
       complex*16,dimension(4,4) :: csq
       complex*16,dimension(0:3,4,4) :: ver12,ver12p 

      call c_gs(sq,csq)
      do mu=0,3
          do i=1,4
           do j=1,4
              c_sigma_q = 0.d0
             do nu=0,3
              c_sigma_q = c_sigma_q + c_gsigma(mu,nu,i,j)*sq(nu)*g(nu,nu)      ! sigma(mu,nu)*q(nu)
             enddo      

             c_ga_5 = 0.d0
             do k=1,4
               c_ga_5 =c_ga_5 + c_ga(mu,i,k)*c_ga5(k,j)                         ! ga(mu)*ga5
             enddo

             ver12p(mu,i,j)=f1*( csq(i,j)*sq(mu)-FMV(sq,sq)*c_ga(mu,i,j) ) /(4.d0*xmn2)                &          
                       + f2*cu*c_sigma_q/(2.d0*xmn)  + c_ga_5*fa      
          enddo

          do j2=1,4       

             if(nparity==1) then  
               ver12(mu,i,j2)=ver12p(mu,i,j2)                                       ! For the positive parity P11(1440)
             else if(nparity==-1) then
               ver12(mu,i,j2)=0.d0
               do k2=1,4
                ver12(mu,i,j2)= ver12(mu,i,j2) + ver12p(mu,i,k2)*c_ga5(k2,j2)        ! For the negative parity S11(1535)  
               enddo
             endif
  
           enddo
         enddo
       enddo  

       END SUBROUTINE
!**************V-A(3/2)  vertex32***************************
       SUBROUTINE vertex32(nparity,fcv,fca,sp,sq,cver32)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: sq,sp,spq
       real*8,dimension(3:6) :: fcv,fca
       complex*16,dimension(4,4) :: csq,cvector0,caxial0
       complex*16,dimension(0:3,0:3,4,4) :: cver32,caxial,cvector
       real*8,external :: FMV

       fcv(6)=0.d0 

       call c_gs(sq,csq)
       spq=sp+sq

       do n1=0,3
        do n2=0,3
          cvector1=fcv(4)/xmn2*(g(n1,n2)*FMV(sq,spq)-sq(n1)*spq(n2))              &   ! the vector part
              + fcv(5)/xmn2*(g(n1,n2)*FMV(sq,sp)-sq(n1)*sp(n2))  !+fcv(6)*g(n1,n2)

          caxial1= fca(4)/xmn2*(g(n1,n2)*FMV(sq,spq) - sq(n1)*spq(n2))          &   ! the minus axial part
                      +fca(5)*g(n1,n2) +fca(6)/xmn2*sq(n1)*sq(n2)  

          DO i=1,4   
            do j2=1,4
             cvector2=fcv(3)/xmn*(g(n1,n2)*csq(i,j2)-sq(n1)*c_ga(n2,i,j2))
             cvector0(i,j2)=cvector2 + cvector1*unm(i,j2)

            ! caxial2=fca(3)/xmn*(g(n1,n2)*csq(i,j2)-sq(n1)*c_ga(n2,i,j2))   ! fca(3)=0
            ! caxial0(i,j2)=caxial2 + caxial1*unm(i,j2) 
              caxial0(i,j2)=caxial1*unm(i,j2)    
            enddo

           if(nparity==1) then          ! the parity of P33(1232) is positive
             DO j=1,4 
               cvector(n1,n2,i,j)=0.d0
               do k2=1,4
                cvector(n1,n2,i,j)=cvector(n1,n2,i,j)+ cvector0(i,k2)*c_ga5(k2,j)
               enddo  
                cver32(n1,n2,i,j)= cvector(n1,n2,i,j) + caxial0(i,j) 
             enddo
           else if(nparity==-1) then   ! the parity of D13(1520) is negative
             do j=1,4
               caxial(n1,n2,i,j)=0.d0
              do k2=1,4
                caxial(n1,n2,i,j)=caxial(n1,n2,i,j)+ caxial0(i,k2)*c_ga5(k2,j)
              enddo  
                cver32(n1,n2,i,j)= caxial(n1,n2,i,j) + cvector0(i,j) 
            enddo
           else
              print *,'parity of the P33(1232) or D13(1520) is not correct'
              stop
           endif
        
       ENDDO 
     
       enddo
       enddo

       END SUBROUTINE
!***********For D13(1520)*****************************************
       SUBROUTINE ver_j32(nparity,t1,ch_verj32,ch_verj32_t)    ! nparity=1 for P33(1232), nparity=-1 for D13(1520)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       complex*16,dimension(0:3,0:3,4,4) :: ch_verj32,ch_verj32_t,cAem,cAem_0,cAem_t,cAnc,cAnc_0,cAnc_t
       complex*16,dimension(4,4) :: clpq,clppq,cg_c,cg_d,cgtem_c,cgtem_d,ca_emt,ca_nc,ca_nct,ca_em
       real*8,dimension(2,6) :: fem,fvem,fvnc,fanc
       real*8,dimension(3:6) :: fcv,fcat,fcvt,fca

       if(nparity==1) then
         call ffP33(t1,fem,fvem,fvnc,fanc)
         nexcit=1
       else if(nparity==-1) then
         call ffD13(t1,fem,fvem,fvnc,fanc)
         nexcit=3
       else 
         print*,'nparity should be 1 or -1. Here is Form Factor of J=3/2'
         stop
       endif 
     
       if(nexcit==1) then
          ii=1
       else if(nexcit==3) then
          ii=-(nucleon-3)/2       ! For neuturon, nucleon=-1, i=2, for proton, nucleon=1, i=1
       else
          print *,'nexcit should be 1 or 3. Here is ver_j32'
          stop
       endif 
        
       do j=3,6
         fcv(j)=fem(ii,j)
         fcvt(j)=fvnc(ii,j)
         fcat(j)=fanc(ii,j)
         fca(j)=0.d0
        ! print *
       enddo
       
      ! for the direct diagram
       call vertex32(nparity,fcv,fca,xp,-xqf,cAem)
       call vertex32(nparity,fcv,fca,xpp,xqf,cAem_0)
       call conj(cAem_0,cAem_t)
      ! For the crossed diagram
       call vertex32(nparity,fcvt,fcat,xp,xq,cAnc)
       call vertex32(nparity,fcvt,fcat,xpp,-xq,cAnc_0)
       call conj(cAnc_0,cAnc_t)
      
       call propagator(nexcit,xpd,cdpq)
       call propagator(nexcit,xpdc,cdppq)   
   
       do nmiu=0,3
        do nalpha=0,3

          cgtem_d=unm0        ! For the direct diagram
          cgtem_c=unm0        ! For the crossed diagram

          do ndelta=0,3
            do nsigma=0,3
          ! For the direct diagram
            call matrix_42(cAem_t,ndelta,nmiu,cA_emt)
            call matrix_42(cAnc,nsigma,nalpha,cA_nc)

            call Lamda(ndelta,nsigma,xpd,clpq)

            call Mult_matrix(cA_emt,clpq,cA_nc,cg_d)
               
            cgtem_d =cgtem_d  + cg_d *g(ndelta,ndelta)*g(nsigma,nsigma)

          ! For the crossed diagram 
            call matrix_42(cAnc_t,ndelta,nalpha,cA_nct)
            call matrix_42(cAem,nsigma,nmiu,cA_em)


            call Lamda(ndelta,nsigma,xpdc,clppq)

            call Mult_matrix(cA_nct,clppq,cA_em,cg_c)     

            cgtem_c = cgtem_c  + cg_c*g(ndelta,ndelta)*g(nsigma,nsigma)
           ENDDO	
          ENDDO   

          DO i=1,4
           DO j=1,4
             ch_verj32(nmiu,nalpha,i,j) =cgtem_d(i,j) *cdpq +cgtem_c(i,j)*cdppq
           enddo
         enddo
            
         ENDDO
        ENDDO

       call conj(ch_verj32,ch_verj32_t)

         Return
      END SUBROUTINE
!***************The propagator of the exicited resonances****************************************
       SUBROUTINE propagator(nexcit,sp,cprop)       
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: sp
       real*8,external :: FMV
       complex*16,external :: cDdelta
       
       select case(nexcit)
       case(0)
         smr=xmn
       case(1) 
         smr=xmp33
       case(2)
         smr=xmp11
       case(3)
         smr=xmd13
       case(4)
         smr=xms11
       case default
         print*,'nexcit should be 0 to 4',nexcit
         stop   
       end select

       sp2=FMV(sp,sp)

       if(nexcit==1) then
        ! call medium(sp2,cprop)         ! for the Delta, we use the deltamedium.f90
          cprop=cDdelta(sp)
       else
         call width(nexcit,sp2,wid)
         cprop=1.d0/(sp2-smr**2+cu*wid*smr)
       endif   

      END SUBROUTINE
!************************************************************
!    width.f90                     2012-11-27
!    Give the width of the excited baryonic resances:  
!         P11(1440),D13(1520),S11(1535)
!
!          wid: the total witdth
!        pwid*: the partial width
!           *0: the width when the invariance mass equal the resonance mass 
!
!************************************************************
       Subroutine width(nexcit,sp2,wid) 
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(4,5) :: widd
       real*8,external :: flam
      ! xmsigma=0.475           ! For Sigma-> Pi Pi in S-wave
      ! xmeta=0.548


       select case(nexcit)  ! nexcit=0,1,2,3,4 for N,Delta,P11(1440),D13(1520),S11(1535)
     ! N or P
       case(0)
         wid=0.d0
     ! P33(1232)
       case(1)   
         if((sp2-(xmn+xmpi)**2).GT.0) then
           !fstar=2.14
           slam=flam(sp2,xmpi**2,xmn2)
            wd=1.d0/(6.d0*pi)*(fstar/xmpi)**2 *xmn/sp2**2  *                                 &
           (sqrt(slam)/(2.d0))**3  
        else 
           wd=0.d0
        endif

        smr=xmd
        pwidnpi0=0.117
        call pwidth1(sp2,smr,xmn,xmpi,1,pwidnpi0,pwidnpi)                           ! decay mode N-Pi
        wid=wd
      !  widd(nexcit,1)=pwidnpi 

      !   print *,'the propagator should be obtained by deltamedium.f90'
      !   stop       
     ! P11(1440) 
       case(2)
         smr=xmp11
         pwidnpi0=0.3*0.65  
         call pwidth1(sp2,smr,xmn,xmpi,1,pwidnpi0,pwidnpi)                           ! decay mode N-Pi
         pwiddpi0=0.3*0.2
         pwiddnpi0=0.117            ! the width of Delta 
         call pwidth2(sp2,smr,xmd,xmpi,1,pwiddpi0,xmn,xmpi,1,pwiddnpi0,pwiddpi)      ! decay mode Delta-pi
         pwidns0=0.3*0.15
         pwidspi0=0.6               ! the width of Sigma->Pi Pi
         call pwidth2(sp2,smr,xmsigma,xmn,0,pwidns0,xmpi,xmpi,0,pwidspi0,pwidns)     ! decay mode N-Sigma
       !  widd(nexcit,1)=pwidnpi
       !  widd(nexcit,2)=pwiddpi
       !  widd(nexcit,3)=pwidns
        
         wid= pwidnpi +pwiddpi + pwidns
         !print *,pwiddpi,pwidnpi,pwidns,j1,i
       !  stop
     ! D13(1520)
      case(3)   
         smr=xmd13        
         pwidnpi0=0.115*0.6
         call pwidth1(sp2,smr,xmn,xmpi,2,pwidnpi0,pwidnpi)                             ! n_pi
         pwiddpi00=0.115*0.15
         pwiddpi20=0.115*0.125     
         pwiddnpi0=0.117
         call pwidth2(sp2,smr,xmd,xmpi,0,pwiddpi00,xmn,xmpi,1,pwiddnpi0,pwiddpi0)    ! decay mode Delta-pi L=0
         call pwidth2(sp2,smr,xmd,xmpi,2,pwiddpi20,xmn,xmpi,1,pwiddnpi0,pwiddpi2)    ! decay mode Delta-pi L=2
         pwidnr00=0.115*0.09     
         pwidnr20=0.115*0.035
         pwidrpi0=0.1491      ! the width of Rho->pi pi P-wave 
         call pwidth2(sp2,smr,xmrho,xmn,0,pwidnr00,xmpi,xmpi,1,pwidrpi0,pwidnr0)       ! decay mode N-Rho L=0
         call pwidth2(sp2,smr,xmrho,xmn,2,pwidnr20,xmpi,xmpi,1,pwidrpi0,pwidnr2)       ! decay mode N-Rho L=2
            
        ! widd(nexcit,1)=pwidnpi
        ! widd(nexcit,2)=pwiddpi0
        ! widd(nexcit,3)=pwiddpi2
        ! widd(nexcit,4)=pwidnr0
        ! widd(nexcit,5)=pwidnr2

         wid= pwidnpi + pwiddpi0 + pwiddpi2 + pwidnr0 + pwidnr2
     ! S11(1535)
        case(4)
        smr=xms11
        pwidnpi0=0.15*0.45
        call pwidth1(sp2,smr,xmn,xmpi,0,pwidnpi0,pwidnpi)                             ! N-Pi
        pwidne0=0.15*0.42
        call pwidth1(sp2,smr,xmn,xmeta,0,pwidne0,pwidne)                              ! N-eta
        pwidnspi0=0.15*0.08
        pwidns0=0.d0   ! call the width of P11(1440)
        call pwidth2(sp2,smr,xmp11,xmpi,0,pwidnspi0,xmn,xmpi,1,pwidns0,pwidnspi)      ! N*(1440)-Pi   
        pwiddpi0=0.15*0.01
        pwiddnpi0=0.117
        call pwidth2(sp2,smr,xmd,xmpi,2,pwiddpi0,xmn,xmpi,1,pwiddnpi0,pwiddpi)     ! Delta-Pi L=2
        pwidnr0=0.15*0.02
        pwidrpi0=0.1491
        call pwidth2(sp2,smr,xmrho,xmn,0,pwidnr0,xmpi,xmpi,1,pwidrpi0,pwidnr)       ! decay mode N-Rho L=2
         pwidns0=0.15*0.02
         pwidspi0=0.6               ! the width of Sigma->Pi Pi
         call pwidth2(sp2,smr,xmsigma,xmn,0,pwidns0,xmpi,xmpi,0,pwidspi0,pwidns)     ! decay mode N-Sigma       

         !widd(nexcit,1)=pwidnpi
         !widd(nexcit,2)=pwidne
         !widd(nexcit,3)=pwidnspi
         !widd(nexcit,4)=pwiddpi
         !widd(nexcit,5)=pwidnr

        wid=pwidnpi + pwidne  + pwiddpi + pwidnspi + pwidnr + pwidns      
        END select 
       END SUBROUTINE
!****************width of the decay mode******************************
       SUBROUTINE  pwidth1(sp2,smr,sma,smb,l,pwid0,pwid)   ! Both of the particle A and B are stable
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,external :: pcm     
       If(sp2.LE.(sma+smb)**2) then
         pwid=0.d0  
       else
         sw=sqrt(sp2) 
         wpcm=pcm(sw,sma,smb)
         rpcm=pcm(smr,sma,smb) 
         rho=wpcm**(2*l+1)/sw
         rho0=rpcm**(2*l+1)/smr  
         pwid=pwid0*rho/rho0
       endif

       END SUBROUTINE
!*********** the momentum in the center mass frame ******************
       FUNCTION Pcm(sw,sma,smb)
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)

       pcm=sqrt( (sw**2-sma**2-smb**2)**2-4.d0*sma**2*smb**2 )/2.d0/sw
       RETURN

       END FUNCTION

!****************width of the decay mode******************************
       SUBROUTINE  pwidth2(sp2,smr,sma,smb,l,pwid1,smc,smd,l2,pwid2,pwid)   ! Both of the particle A and B are stable
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,external :: pcm   
        
       If(sp2.LE.(smc+smd+smb)**2) then
         pwid=0.d0 
       else 
         sw=sqrt(sp2)
         call rho_width(sw,sma,smb,l,smc,smd,l2,pwid2,rho)
         call rho_width(smr,sma,smb,l,smc,smd,l2,pwid2,rho0)
         pwid=pwid1*rho/rho0  
        ! print *,rho,rho0,pwid1,pwid,xw,sma,xmb  
        ! stop  
       endif
       end subroutine
!********************************************
       SUBROUTINE rho_width(sw,sma,smb,l,smc,smd,l2,pwid2,rho)  
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(200) :: x,f1
       real*8,external :: pcm   
        pi=acos(-1.d0)  
      
         pmin=0.d0
         pmax=pcm(sw,smc+smd,smb)
        n=1
        np=20*n
        call DSG20r(pmin,pmax,n,x,np)
        do i=1,np
         f1(i)=fkernel_rho(x(i))
        ! print *,f1(i)
        end do
        ! stop
        call DRG20r(pmin,pmax,n,f1,rho)

       contains

         FUNCTION fkernel_rho(pab)
         implicit real*8 (a,b,d-h,o-z)
         implicit complex*16 (c)
        ! real*8,dimension(4,5) :: swidd
         swa=sqrt(sw**2+smb**2-2.d0*sw*sqrt(pab**2+smb**2) ) 
         if(abs(sma-1.44).gt.1d-5) then
           call pwidth1(swa**2,sma,smc,smd,l2,pwid2,pwida)
          ! print *,swa,sma,smc,smd,l2,pwid2,pwida
          ! stop
         else 
           call width(2,swa**2,pwida) 
         endif
 
         fkernel_rho=2.d0/pi*pab**(2*l+2)/sqrt(pab**2+smb**2)*sma*pwida  &
             /( (swa**2-sma**2)**2 + sma**2*pwida**2 ) 
     !    print*,xwa,sma,xmc,xmdd,l2,pwid2,pwida,fkernel_rho,pab
     !    stop
         return
         end function

       END SUBROUTINE
!*********************************

!********************************************************
!    deltamedium.f90           2012-09-17
!
!
!
!********************************************************
       SUBROUTINE medium(spd,cpropagator)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(0:3) :: spd
       real*8,external :: fct,FMV

       rho=Ddensity 

       rho0=0.17 *0.197327**3 ! / fm_gev**3     ! nuclear matter density (in fm^-3)
       gprime=0.63
       crho=2.
       flamrho2=2.5**2            ! Lambda_rho **2

       pf=(3.*pi**2*rho/2.)**(.3333333)   ! Fermi momentum from the density
       rat=rho/rho0

       sp2=FMV(spd,spd)
       spn2=spd(0)**2-sp2


!     Calculation of the pauli-blocked width
       

        ! if (gamdfree.eq.0.) then 
         if (sp2.le.(xmn+xmpi)**2) then
            gamdpb=0.
            gamdfree=0. 
         else 
!           Approximation from Nieves et al. NPA 554(93)554     
            qcm=sqrt( (sp2-(xmpi+xmn)**2)*(sp2-(xmpi-xmn)**2)   /  &     ! the centre momentume of the Delta rest frame
               (4.*sp2) )  
            gamdfree=gamd(sp2,qcm)

            if(qcm.le.0.) then
              print *,qcm
              stop
            endif 
           
               r=qcm/pf
               if (r.gt.1) then
                  f=1.+(-2./5./r**2+9./35./r**4-2./21./r**6)
               else
                  f=34./35.*r-22./105*r**3
               end if
              if(f.gt.1) then
               print *,'f gt 1'
               stop
              endif
              gamdpb=gamdfree*f
         end if

!     sig_re=0.  PRC83,045502(2011)
       !  sq2=FMV(xq,xq)
       !  sqn2=xq(0)**2-sq2 

!        Calculation of the delta selfenergy
!        Real part
       !  if(nmedium==2) then
       !      sig_re=-0.053*rat
       !  else     
       !      frho=(flamrho2-xmrho2)/(flamrho2-sq2)          
       !      vt=sqn2/(sq2-xmrho2)*crho*frho**2+gprime         
       !      sig_re=-0.053*rat +4./9.*(fstar/xmpi)**2*rho*vt
       ! endif

!        Imaginary part: using Oset, Salcedo, NPA 468(87)631
!        Using eq. (3.5) to relate the energy of the delta with the pion energy used 
!        in the parametrization

!        Prescription for the effective pion energy
      !   ome=ppi(0)   ! because it coincides with the energy of the outgoing pion in the main program
         if (sp2.le.(xmn+xmpi)**2) then
             sig_im=0.
         else 
             ome=(sp2-xmn2-(3.*pf**2/10./xmn)**2-xmpi**2) /   &
                  (3.*pf**2/10./xmn+xmn)           

!           The parameterization is valid for  85 MeV < tpi < 315. above  we take a contant values
            if (ome.ge.(xmpi+0.315)) ome=xmpi+0.315
 
!           The parameterization of Oset, Salcedo, with ca3 extrapolated to zero at low kin. energies
            fcq=fct(-5.19d0,15.35d0,2.06d0,ome)/1000.
            fca2=fct(1.06d0,-6.64d0,22.66d0,ome)/1000.
            fca3=fct(-13.46d0,46.17d0,-20.34d0,ome)/1000.
            falpha=fct(0.382d0,-1.322d0,1.466d0,ome)
            fbeta=fct(-0.038d0,0.204d0,0.613d0,ome)
            fgamma=2.*fbeta

            if (ome.le.(xmpi+0.085)) fca3=fct(-13.46d0,46.17d0,-20.34d0,(xmpi+0.085d0))/0.085*(ome-xmpi)/1000.
         
            sig_im=-(fcq*rat**falpha+fca2*rat**fbeta+fca3*rat**fgamma) 

         endif

       If(nmedium==1) then
          cpropagator=1./(sp2-xmd**2+cu*xmd*(gamdpb-2.*sig_im))
       else if(nmedium==0) then
          cpropagator=1./(sp2-xmd**2+cu*xmd*gamdfree)  
       endif


       end subroutine    

       function fct(a,b,fc,ome)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)

         x=ome/0.138-1.           
         fct=a*x**2+b*x+fc
       end function


!   The width of free
       function gamd(s,sqcm)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
         
        ! if (s.le.(xmn+xmpi)**2) then
        !    gamd=0.
        ! else
            gamd=1./(6.*pi)*(fstar/xmpi)**2*xmn*sqcm**3/sqrt(s)
        ! end if

       end function gamd
       subroutine deltaff(q2,c3v,c4v,c5v,c3a,c4a,c5a,c6a)
!      N-Delta(1232) transition form factors 
!      using MAID2007 in the vector sector
!      IMPORTANT: these are the CC form factors: CiV and CiA in Tables 5.5 and 5.6 in Tina's thesis 
!      Rules to obtain EM and NC form factors can be found in these tables

         implicit none
         real*8,intent(in)::q2  ! q2=t <0 4-momentum transfer 
         real*8,intent(out)::c3v,c4v,c5v  ! vector form factors (dimensionless)
         real*8,intent(out)::c3a,c4a,c5a,c6a  ! axial form factors (dimensionless)
         
         real*8,parameter::pi=4.*atan(1.)
         real*8,parameter::alpha=1./137.    ! fine structure constant
         real*8,parameter::md=1.232         !  delta mass
         real*8,parameter::mn=.939          ! nucleon mass
         real*8,parameter::mpi=.139         ! pion mass
         
         real*8,parameter::Ma=0.93    ! N-Delta axial mass

         real*8 ::a12,a32,s12         ! helicity amplitudes
         real*8 ::r,kr,mpw2,mmw2

         real*8 ::Fd  ! dipole

!        vector form factors          
         
         call heliamp07('P33(1232)','p',q2,a12,a32,s12)  ! nucleon label 'p' irrelevant in this case
        ! a12=a12*1.e-3
        ! a32=a32*1.e-3
        ! s12=s12*1.e-3     ! (-) gives the solution like in Tina's thesis
           

         kr=(md**2-mn**2)/2./md
         mpw2=(mn+md)**2
         mmw2=(mn-md)**2

         r=sqrt(6.*kr/pi/alpha*mn*md/(mmw2-q2))

         c3v=-r*mn*md/(mpw2-q2)*(a12+a32/sqrt(3.))

         c4v=-2.*r*mn**2/(mpw2-q2)/(mmw2-q2)*(mn*md*a12 &
              & -(md**2+mn**2-mn*md-q2)*a32/sqrt(3.) &
             ! & - sqrt(2.)*md**2*(md**2-mn**2-q2)/sqrt((mpw2-q2)*(mmw2-q2))*s12)
              & + sqrt(2.)*md**2*(md**2-mn**2-q2)/sqrt((mpw2-q2)*(mmw2-q2))*s12)

         c5v=-2.*r*(md*mn)**2/(mpw2-q2)/(mmw2-q2)*(-a12+a32/sqrt(3.) &
               & -sqrt(2.)*(md**2-mn**2+q2)/sqrt((mpw2-q2)*(mmw2-q2))*s12)
             !  & +sqrt(2.)*(md**2-mn**2+q2)/sqrt((mpw2-q2)*(mmw2-q2))*s12)

!        axial form factors
!        simple dipole for c5a and adler

         Fd=(1.-q2/(Ma**2))**(-2)
         c5a=1.2*Fd        
         c4a=-c5a/4.
         c3a=0.
         c6a=c5a*mn**2/(mpi**2-q2)

       end subroutine deltaff
!===================================================================================
!----------------------Form Factor For P11(1440)------------------------------------   
!===================================================================================
       SUBROUTINE ffp11(t1,fem,fvem,fvnc,fanc)
       use parameter
       implicit real*8 (a,b,d-h,o-y)
       implicit complex*16 (c)
       real*8,dimension(2,6) :: fem,fvem,fvnc
       real*8,dimension(2) :: fanc

         !xma=1.d0    
         smr=xmp11 !1.462
         fa0=fa0p11 !-0.52

       rk=sqrt(xmn*(smr**2-xmn2)/pi/alpha)
       xmp2=(smr+xmn)**2
       xmm2=(smr-xmn)**2 
       Q2=-t1
       xmpq2=(smr+xmn)**2 + Q2
       xmmq2=(smr-xmn)**2 + Q2

     !  a12=a12*1.e-3
     !  a32=a32*1.e-3
     !  s12=s12*1.e-3     ! (-) gives the solution like in Tina's thesis

      ! EM form factor 
       tem=0.d0
       call heliamp07('P11(1440)','p',tem,a12,a32,s12) 
       fem(1,1)=sqrt(8.d0)*xmn2*rk/xmp2/sqrt(xmm2)*( a12 - sqrt(8.d0/xmp2/xmm2)*smr*(smr+xmn)*s12 )
       fem(1,2)=sqrt(2.d0)*xmn*rk/xmp2/sqrt(xmm2)* (smr+xmn)*a12 

       call heliamp07('P11(1440)','n',tem,a12,a32,s12) 
       fem(2,1)=sqrt(8.d0)*xmn2*rk/xmp2/sqrt(xmm2)*( a12 - sqrt(8.d0/xmp2/xmm2)*smr*(smr+xmn)*s12 )
       fem(2,2)=sqrt(2.d0)*xmn*rk/xmp2/sqrt(xmm2)* (smr+xmn)*a12 

      ! NC form factor 
       call heliamp07('P11(1440)','p',t1,a12,a32,s12) 
       fvem(1,1)=sqrt(8.d0)*xmn2*rk/xmpq2/sqrt(xmmq2)*( a12 - sqrt(8.d0/xmpq2/xmmq2)*smr*(smr+xmn)*s12 )
       fvem(1,2)=sqrt(2.d0)*xmn*rk/xmpq2/sqrt(xmmq2)* ( (smr+xmn)*a12 + sqrt(8.d0/xmpq2/xmmq2)*smr*Q2*s12 )

       call heliamp07('P11(1440)','n',t1,a12,a32,s12) 
       fvem(2,1)=sqrt(8.d0)*xmn2*rk/xmpq2/sqrt(xmmq2)*( a12 - sqrt(8.d0/xmpq2/xmmq2)*smr*(smr+xmn)*s12 )
       fvem(2,2)=sqrt(2.d0)*xmn*rk/xmpq2/sqrt(xmmq2)* ( (smr+xmn)*a12 + sqrt(8.d0/xmpq2/xmmq2)*smr*Q2*s12 )

       call ffem_nc(2,fvem,fvnc)
             
       fanc(1)=fa0*(1.+q2/x2ma2)**(-2) *0.5
       fanc(2)=-fanc(1)


       END SUBROUTINE
!===================================================================================
!----------------------Form Factor For S11(1535)------------------------------------   
!===================================================================================
       SUBROUTINE ffS11(t1,fem,fvem,fvnc,fanc)
       use parameter
       implicit real*8 (a,b,d-h,o-y)
       implicit complex*16 (c)
       real*8,dimension(2,6) :: fem,fvem,fvnc
       real*8,dimension(2) :: fanc

       !  xma=1.d0    
         smr=xms11  !1.535
         fa0=fa0s11 ! -0.23

       rk=sqrt(xmn*(smr**2-xmn2)/pi/alpha)
       xmp2=(smr+xmn)**2
       xmm2=(smr-xmn)**2 
       Q2=-t1
       xmpq2=(smr+xmn)**2 + Q2
       xmmq2=(smr-xmn)**2 + Q2


      ! EM form factor 
       tem=0.d0
       call heliamp07('S11(1535)','p',tem,a12,a32,s12) 
       fem(1,1)=sqrt(8.d0)*xmn2*rk/xmm2/sqrt(xmp2)*( a12 + sqrt(8.d0/xmp2/xmm2)*smr*(smr-xmn)*s12 )
       fem(1,2)=sqrt(2.d0)*xmn*rk/xmm2/sqrt(xmp2)* (smr-xmn)*a12 

       call heliamp07('S11(1535)','n',tem,a12,a32,s12) 
       fem(2,1)=sqrt(8.d0)*xmn2*rk/xmm2/sqrt(xmp2)*( a12 + sqrt(8.d0/xmp2/xmm2)*smr*(smr-xmn)*s12 )
       fem(2,2)=sqrt(2.d0)*xmn*rk/xmm2/sqrt(xmp2)* (smr-xmn)*a12 

      ! NC form factor 
       
       call heliamp07('S11(1535)','p',t1,a12,a32,s12) 
       fvem(1,1)=sqrt(8.d0)*xmn2*rk/xmmq2/sqrt(xmpq2)*( a12 + sqrt(8.d0/xmpq2/xmmq2)*smr*(smr-xmn)*s12 )
       fvem(1,2)=sqrt(2.d0)*xmn*rk/xmmq2/sqrt(xmpq2)* ( (smr-xmn)*a12 - sqrt(8.d0/xmpq2/xmmq2)*smr*Q2*s12 )

       call heliamp07('S11(1535)','n',t1,a12,a32,s12) 
       fvem(2,1)=sqrt(8.d0)*xmn2*rk/xmpq2/sqrt(xmmq2)*( a12 + sqrt(8.d0/xmpq2/xmmq2)*smr*(smr-xmn)*s12 )
       fvem(2,2)=sqrt(2.d0)*xmn*rk/xmpq2/sqrt(xmmq2)* ( (smr-xmn)*a12 - sqrt(8.d0/xmpq2/xmmq2)*smr*Q2*s12 )
        
       call ffem_nc(2,fvem,fvnc)
     
       fanc(1)=fa0*(1.+q2/x2ma2)**(-2) *0.5
       fanc(2)=-fanc(1)

       END SUBROUTINE
!================================================================================
!-------------------  Transfer from EM FF to NC FF------------------------------
!================================================================================
       Subroutine ffem_nc(n,fvem,fvnc)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(2,6) :: fvem,fvnc

       sw=0.5d0-2.d0*sin_theta_w2
       
       if(n==2) then
         m1=1
         m2=2
       else if(n==3) then          
         m1=3
         m2=5
       endif
       !  print *,m1,m2,sw
       !   stop
       do i=m1,m2
         fvnc(1,i)=sw*fvem(1,i)-0.5*fvem(2,i)
         fvnc(2,i)=sw*fvem(2,i)-0.5*fvem(1,i)
       enddo
         
       END SUBROUTINE    
!===================================================================================
!----------------------Form Factor For D13(1520)------------------------------------   
!===================================================================================

       SUBROUTINE ffD13(t1,fem,fvem,fvnc,fanc)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(2,6) :: fem,fvem,fvnc,fanc

       fca50=fca5d13 !-2.15
       smr=xmd13 !1.524
     !  x2ma2=1.d0

       rk=sqrt(xmn*(smr**2-xmn2)/pi/alpha)
       xmp2=(smr+xmn)**2
       xmm2=(smr-xmn)**2 
       Q2=-t1
       xmpq2=(smr+xmn)**2 + Q2
       xmmq2=(smr-xmn)**2 + Q2


      ! EM form factor 
       tem=0.d0
       call heliamp07('D13(1520)','p',tem,a12,a32,s12) 
        
       fem(1,3)=rk/sqrt(xmp2)*xmn*smr/xmm2*(sqrt(3.d0)*a12-a32)
       fem(1,4)=2.d0*rk/xmp2/sqrt(xmp2)*xmn2/xmm2 *( -sqrt(3.d0)*smr*xmn*a12 &
            +(smr**2+smr*xmn+xmn2)*a32 - sqrt(6.d0)*smr**2*(smr**2-xmn2)/sqrt(xmp2*xmm2)*s12 )
       fem(1,5)=2.d0*rk/xmp2/sqrt(xmp2)*smr**2*xmn2/xmm2 *( -sqrt(3.d0)*a12 &
            - a32 + sqrt(6.d0)*(smr**2-xmn2)/sqrt(xmp2*xmm2)*s12 )
       
       call heliamp07('D13(1520)','n',tem,a12,a32,s12)
        
       fem(2,3)=rk/sqrt(xmp2)*xmn*smr/xmm2*(sqrt(3.d0)*a12-a32)
       fem(2,4)=2.d0*rk/xmp2/sqrt(xmp2)*xmn2/xmm2 *( -sqrt(3.d0)*smr*xmn*a12 &
            +(smr**2+smr*xmn+xmn2)*a32 - sqrt(6.d0)*smr**2*(smr**2-xmn2)/sqrt(xmp2*xmm2)*s12 )
       fem(2,5)=2.d0*rk/xmp2/sqrt(xmp2)*smr**2*xmn2/xmm2 *( -sqrt(3.d0)*a12 &
            - a32 + sqrt(6.d0)*(smr**2-xmn2)/sqrt(xmp2*xmm2)*s12 )

      ! NC form factor
       tem=t1
       call heliamp07('D13(1520)','p',tem,a12,a32,s12) 
       fvem(1,3)=rk/sqrt(xmpq2)*xmn*smr/xmmq2*(sqrt(3.d0)*a12-a32)
       fvem(1,4)=2.d0*rk/xmpq2/sqrt(xmpq2)*xmn2/xmmq2 *( -sqrt(3.d0)*smr*xmn*a12 &
            +(smr**2+smr*xmn+xmn2+q2)*a32 - sqrt(6.d0)*smr**2*(smr**2-xmn2+q2)/sqrt(xmpq2*xmmq2)*s12 )
       fvem(1,5)=2.d0*rk/xmpq2/sqrt(xmpq2)*smr**2*xmn2/xmmq2 *( -sqrt(3.d0)*a12 &
            - a32 + sqrt(6.d0)*(smr**2-xmn2-q2)/sqrt(xmpq2*xmmq2)*s12 )
       
       call heliamp07('D13(1520)','n',tem,a12,a32,s12) 
       fvem(2,3)=rk/sqrt(xmpq2)*xmn*smr/xmmq2*(sqrt(3.d0)*a12-a32)
       fvem(2,4)=2.d0*rk/xmpq2/sqrt(xmpq2)*xmn2/xmmq2 *( -sqrt(3.d0)*smr*xmn*a12 &
            +(smr**2+smr*xmn+xmn2+q2)*a32 - sqrt(6.d0)*smr**2*(smr**2-xmn2+q2)/sqrt(xmpq2*xmmq2)*s12 )
       fvem(2,5)=2.d0*rk/xmpq2/sqrt(xmpq2)*smr**2*xmn2/xmmq2 *( -sqrt(3.d0)*a12 &
            - a32 + sqrt(6.d0)*(smr**2-xmn2-q2)/sqrt(xmpq2*xmmq2)*s12 )
       
       call ffem_nc(3,fvem,fvnc)     
       fanc(1,3)=0.d0
       fanc(1,5)=fca50*(1.d0+q2/x2ma2)**(-2)*0.5    
       fanc(1,6)=xmn2/(q2+xmpi**2)*fanc(1,5)
       fanc(1,4)=-0.25*fanc(1,5)

       do i=3,6
         fanc(2,i)=-fanc(1,i)
       enddo 
         ! Luis choice
          fanc(1,4)=0.d0
          fanc(2,4)=0.d0 
       END SUBROUTINE
!===================================================================================
!----------------------Form Factor For P33(1232)------------------------------------   
!===================================================================================

       SUBROUTINE ffP33(t1,fem,fvem,fvnc,fanc)
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)
       real*8,dimension(2,6) :: fem,fvem,fvnc,fanc

       fca50=fca5p33  !1.2
       smr=xmp33 !1.232
       !xma=0.93
       sw=1.-2.*sin_theta_w2
     !  print*,fvnc

       rk=sqrt(xmn*(smr**2-xmn2)/pi/alpha)
       xmp2=(smr+xmn)**2
       xmm2=(smr-xmn)**2 
       Q2=-t1
       xmpq2=(smr+xmn)**2 + Q2
       xmmq2=(smr-xmn)**2 + Q2


      ! EM form factor 
       tem=0.d0
       call heliamp07('P33(1232)','p',tem,a12,a32,s12)
         
       fem(1,3)=rk/sqrt(xmm2)*xmn*smr/xmp2*(sqrt(3.d0)*a12+a32)
       fem(1,4)=2.d0*rk/xmm2/sqrt(xmm2)*xmn2/xmp2 *( sqrt(3.d0)*smr*xmn*a12 &
            -(smr**2-smr*xmn+xmn2)*a32 + sqrt(6.d0)*smr**2*(smr**2-xmn2)/sqrt(xmp2*xmm2)*s12 )
       fem(1,5)=2.d0*rk/xmm2/sqrt(xmm2)*smr**2*xmn2/xmp2 *( -sqrt(3.d0)*a12 &
            + a32 - sqrt(6.d0)*(smr**2-xmn2)/sqrt(xmp2*xmm2)*s12 )
       
      ! NC form factor
       tem=t1
       call heliamp07('P33(1232)','p',tem,a12,a32,s12) 
       
       fvem(1,3)=rk/sqrt(xmmq2)*xmn*smr/xmpq2*(sqrt(3.d0)*a12+a32)
       fvem(1,4)=2.d0*rk/xmmq2/sqrt(xmmq2)*xmn2/xmpq2 *( sqrt(3.d0)*smr*xmn*a12 &
            -(smr**2-smr*xmn+xmn2+q2)*a32 + sqrt(6.d0)*smr**2*(smr**2-xmn2+q2)/sqrt(xmpq2*xmmq2)*s12 )
       fvem(1,5)=2.d0*rk/xmmq2/sqrt(xmmq2)*smr**2*xmn2/xmpq2 *( -sqrt(3.d0)*a12 &
            + a32 - sqrt(6.d0)*(smr**2-xmn2-q2)/sqrt(xmpq2*xmmq2)*s12 )
       
       !print *,fem(1,3),tem,a12,a32,s12
       !stop

         do j=3,5
             fvnc(1,j)=(1.d0-2.d0*sin_theta_w2)*fvem(1,j)
         enddo

       fanc(1,3)=0.d0
       fanc(1,5)=fca50*(1.+q2/xma2)**(-2)  * (-1.d0)    !  Fnc = -Fcc
       fanc(1,6)=xmn2/(q2+xmpi**2)*fanc(1,5)
       fanc(1,4)=-0.25*fanc(1,5)

       
       END SUBROUTINE
!*****************************************************************************
!    Helicity amplitude
!
!****************************************************************************
!     Helicity amplitides for nucleon-resonance em transitions according to MAID2007
      subroutine heliamp07(res,nucl,q2,A12,A32,S12)
        implicit none
        
        character(len=9),intent(in)::res  ! P33(1232),P11(1440), etc
        character(len=1),intent(in)::nucl ! p or n 
        real*8,intent(in)::q2  ! in GeV2
        real*8,intent(out)::A12,A32,S12  ! helicity ampls. in units of 10-3 GeV^(-1/2)

        real*8,parameter::pi=3.1416
        real*8,parameter::mn=0.9382723

        real*8 ::Q2G   ! -q2
        real*8 ::mres   ! resonance Breit-Wigner mass
        real*8 ::kgcm0  ! energy of an equivalent real photon in cm
        real*8 ::egcm,qcm    ! photon energy and momentum in cm
        
        real*8 ::AE,AM,AS,Fq,tau

        Q2G=-q2

        select case(res)
        case ('P33(1232)')

           mres=1.232
           kgcm0 = (mres**2-mn**2)/2./mres
           egcm = (mres**2-Q2G-mn**2)/2./mres
           qcm=sqrt(egcm**2+Q2G)
           tau=Q2G/4./mn**2

           Fq=1./(1+Q2G/0.71)**2*(qcm/kgcm0)

           AM=300*(1.+ 0.01*Q2G)*exp(-0.23*Q2G)*Fq
           AE=-6.37*(1.- 0.021*Q2G)*exp(-0.16*Q2G)*Fq

           AS=-12.40*(1.+0.12*Q2G)/(1.+4.9*tau)*(qcm/kgcm0)*exp(-0.23*Q2G)*Fq

           A12=-.5*(3*AE+AM)
           A32=sqrt(3.)/2.*(AE-AM)
           S12=-sqrt(2.)*AS

        case ('P11(1440)')

           if (nucl.eq.'p') then
              A12= -61.4*(1.- 1.22*Q2G -0.55*Q2G**4)*exp(-1.51*Q2G)
              S12= 4.2*(1.+ 40.*Q2G+1.5*Q2G**4)*exp(-1.75*Q2G)
              A32= 0.
           else
              A12=54.1*(1.+0.95*Q2G)*exp(-1.77*Q2G)
              S12=-41.5*(1.+ 2.98*Q2G)*exp(-1.55*Q2G)
              A32=0.
           end if

        case ('D13(1520)')

           if (nucl.eq.'p') then
              
              A12=-27.*(1.+7.77*Q2G)*exp(-1.09*Q2G)
              A32= 161.*(1.+0.69*Q2G)*exp(-2.1*Q2G)
              S12= -63.6*(1.+4.19*Q2G)*exp(-3.40*Q2G)

           else
                            
              A12= -77.*(1.-0.53*Q2G)*exp(-1.55*Q2G)
              A32= -154.*(1.+0.58*Q2G)*exp(-1.75*Q2G)
              S12= 13.6*(1.+15.7*Q2G)*exp(-1.57*Q2G)

           end if

        case ('S11(1535)')

            if (nucl.eq.'p') then
               A12=66.*(1.+1.61*Q2G)*exp(-0.70*Q2G)
               S12=-2.*(1.+23.9*Q2G)*exp(-0.81*Q2G)
               A32=0.
            else
               A12=-51.*(1.+4.75*Q2G)*exp(-1.69*Q2G)
               S12=28.5*(1.+0.36*Q2G)*exp(-1.55*Q2G)
               A32=0.
            end if
               
         end select
       a12=a12*1.e-3
       a32=a32*1.e-3
       s12=s12*1.e-3  
       
       end subroutine heliamp07





!***************************************************************************
!
!
!
!***************************************************************************
       SUBROUTINE diracmatrix()       
       use parameter
       implicit real*8 (a,b,d-h,o-z)
       implicit complex*16 (c)

!       The following is the definition of diracmatici
      do n=0,3
        do i=1,4
           do j=1,4
           c_ga(n,i,j)=0.d0               
          enddo
        enddo
      enddo

      do i=1,4
        do j=1,4
           c_ga5(i,j)=0.d0
           unm0(i,j)=0.d0
           if(i==j) then
             unm(i,j)=1.d0
            else
             unm(i,j)=0.d0
           endif 
        enddo 
      enddo

      do i=0,3
        do j=0,3

           if(i==j) then
          ! xunm(i,j)=1.d0
             if(i==0) then                               
               g(i,j)=1.d0
             else
               g(i,j)=-1.d0
             endif
           else
           !  xunm(i,j)=0.d0
             g(i,j)=0.d0
           endif
        enddo
      enddo 
   ! define the matrix of Gamma 0
      c_ga(0,1,1)=1.d0
      c_ga(0,2,2)=1.d0
      c_ga(0,3,3)=-1.d0
      c_ga(0,4,4)=-1.d0

   !define the matrix of Gamma 1 
      c_ga(1,1,4)=1.d0
      c_ga(1,2,3)=1.d0
      c_ga(1,3,2)=-1.d0
      c_ga(1,4,1)=-1.d0  

   !define the matrix of Gamma 2 
      c_ga(2,1,4)=-cu
      c_ga(2,2,3)=cu
      c_ga(2,3,2)=cu
      c_ga(2,4,1)=-cu   
 
   !define the matrix of Gamma 3
      c_ga(3,1,3)=1.d0
      c_ga(3,2,4)=-1.d0
      c_ga(3,3,1)=-1.d0
      c_ga(3,4,2)=1.d0   

   !define the matrix of Gamma 5 
      c_ga5(1,3)=1.d0
      c_ga5(2,4)=1.d0
      c_ga5(3,1)=1.d0
      c_ga5(4,2)=1.d0  


    !define the matrix dirac tensor sigma(a,b)=i/2(a*b-b*a)
      
     do n1=0,3
       do n2=0,3
         do i=1,4
           do j=1,4

              c_12=0.d0
              c_21=0.d0
             do k=1,4
              c_12=c_12+c_ga(n1,i,k)*c_ga(n2,k,j)
              c_21=c_21+c_ga(n2,i,k)*c_ga(n1,k,j)
             enddo

            c_gsigma(n1,n2,i,j)=cu/2.d0*(c_12-c_21)
 
           enddo
         enddo
       enddo
     enddo  
!   The above is the definition of diracmatric

!   the following is the lvtensor
      do i=0,3
        do j=0,3
          do k=0,3
            do l=0,3
         !     Lvtensor(i,j,k,l)=Lvten(i,j,k,l)
            enddo
          enddo
         enddo
      enddo



        
       END SUBROUTINE































      

       











           


