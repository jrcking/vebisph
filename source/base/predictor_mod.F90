module predictor_mod
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  use calculation_gradient
  use nnf_models
  use mirror_boundaries
  use omp_lib

  implicit none

  real(rkind),dimension(:,:),allocatable :: rhs
  real(rkind),dimension(dims) :: uij,rhs_tmp
  real(rkind),dimension(dims,dims) :: idm
  real(rkind) :: L0,eps_PTT,alpha_Giesekus
  real(rkind),dimension(:,:,:),allocatable :: phi_evss,grad_vel

  private
  public :: prediction_step

contains
!! ------------------------------------------------------------------------------------------------
  subroutine prediction_step
    integer(ikind) :: i,j,k
    real(rkind), dimension(dims,dims) :: tauij
    real(rkind) :: aatau
    real(rkind) :: coeff_solvent,coeff_polymeric,coeff_evss

    !! Calculate the velocity gradient 
    allocate(grad_vel(np,dims,dims))
    call grad_operator(grad_vel(:,1,:),u(:,1),1)
    call grad_operator(grad_vel(:,2,:),u(:,2),1)

    !! Evolution of the polymeric stress, if viscoelastic model
#if const_mod!=1
    !! Find the new polymeric stress - at t^n+1/2
    call constitutive_eqn_evolve
   
    !! EVSS
    coeff_evss = alpha_evss*El*sqrt(Ra/Pr)/(1.0d0-beta)
    allocate(phi_evss(np,dims,dims))
    !$OMP PARALLEL DO
    do i=1,np
       phi_evss(i,:,:) = tau_p(i,:,:) - coeff_evss*(grad_vel(i,:,:) + transpose(grad_vel(i,:,:)))
    end do
    !$OMP END PARALLEL DO
#endif

    !! Calculate contribution to RHS due to solvent and polymeric stresses
    allocate(rhs(npfb,dims))
    rhs(:,:)=0_rkind
    
    coeff_solvent = (beta+alpha_evss)*sqrt(Pr/Ra)
    coeff_polymeric = (1.0d0-beta)*Pr/(Ra*El)

#if const_mod==1
    !$OMP PARALLEL DO PRIVATE(k,j,uij,rhs_tmp)
    do i=1,npfb                                       !! Newtonian option
       rhs_tmp(:) = 0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k)
          uij(:) = u(i,:) - u(j,:)         
          rhs_tmp(:) = rhs_tmp(:) + coeff_solvent*uij(:)*ij_w_L(i,k)
       end do
       rhs(i,:) = rhs_tmp(:)
    end do
    !$OMP END PARALLEL DO
#else
    !$OMP PARALLEL DO PRIVATE(k,j,uij,rhs_tmp,tauij)
    do i=1,npfb                                      !! Viscoelastic option
       rhs_tmp(:) = 0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k)
          uij(:) = u(i,:) - u(j,:)         
#ifdef use_labfm
          tauij(:,:) = phi_evss(j,:,:) - phi_evss(i,:,:)    
          rhs_tmp(:) = rhs_tmp(:) + coeff_polymeric*matmul(tauij,ij_w_G(i,k,:)) & 
                                   + coeff_solvent*uij(:)*ij_w_L(i,k)
#else          
          tauij(:,:) = phi_evss(j,:,:) - phi_evss(i,:,:)    
          rhs_tmp(:) = rhs_tmp(:) + coeff_polymeric*matmul(tauij,matmul(kgcm(i,:,:),ij_w_G(i,k,:))) & 
                                   + coeff_solvent*uij(:)*ij_w_L(i,k)
#endif                                   
       end do
       rhs(i,:) = rhs_tmp(:)
    end do
    !$OMP END PARALLEL DO
#endif
    
    !! Calculate contribution to RHS of some external forcing if required
    if(external_forcing)then
       do i=1,npfb
          rhs(i,1) = rhs(i,1) + sin(4.0d0*r(i,2))   !! Kolmogorov forcing
       end do
    end if
    
    !! Calculate contribution to RHS of advection terms
    !$OMP PARALLEL DO
    do i=1,npfb
       rhs(i,1) = rhs(i,1) + dot_product(ushift(i,:),grad_vel(i,1,:))
       rhs(i,2) = rhs(i,2) + dot_product(ushift(i,:),grad_vel(i,2,:))
    end do
    !$OMP END PARALLEL DO

    !! Find U* and store in u-array
    !$omp parallel do
    do i=1+nbS,npfb
       u(i,:) = u(i,:) + rhs(i,:)*dt
    end do
    !$omp end parallel do

    deallocate(rhs)
    deallocate(grad_vel)
#if const_mod!=1
    deallocate(phi_evss)
#endif    
  end subroutine prediction_step
!! ------------------------------------------------------------------------------------------------
  subroutine constitutive_eqn_evolve

!    integer(ikind),intent(in) :: pc_flag
    real(rkind), dimension(dims,dims) :: ddt,AmI
    integer(ikind) :: i,j,k
    real(rkind), dimension(:,:),allocatable :: grad_xx,grad_xy,grad_yy
    real(rkind),dimension(:,:,:),allocatable :: adv_term

    real(rkind),dimension(:,:,:),allocatable :: EigenVecs,Psi,Psi0
    real(rkind),dimension(:,:),allocatable :: EigenVals
    real(rkind),dimension(dims,dims) :: A,tmp_mat,B,Omega
    real(rkind),dimension(dims) :: tmp_vec
    real(rkind) :: tmp1,tmp2
    real(rkind) :: alpha_R,beta_R
    

    !! Useful to have an identity matrix
    idm = 0.0d0;idm(1,1) = 1.0d0;idm(2,2) = 1.0d0  

    !! Constitutive eqn parameters (should come from param file in due course)
    L0 = sqrt(ve_nonlinearity)!5.0
    eps_PTT = ve_nonlinearity
    alpha_Giesekus = ve_nonlinearity
       

    allocate(EigenVecs(npfb,dims,dims),EigenVals(npfb,dims));EigenVecs = 0.0d0;EigenVals = 0.0d0
    !$OMP PARALLEL DO PRIVATE(A)
    do i=1,npfb
       !! Step 1: Find the conformation tensor at time n-1/2 (i.e. * in the prev time-step)
       A = f_strain_backward(tau_p(i,:,:),trace_conf(i))

       !! Step 2: Find Eigen Vecs/vals of conformation tensor - A=R(Lambda)R^T  (all at n-1/2)
       call calc_eigens(A(:,:),EigenVals(i,:),EigenVecs(i,:,:))
    end do
    !$OMP END PARALLEL DO


    !! Step 3: Find the log-conformation tensor Psi(n-1/2)=Rlog(Lambda)R^T  (all at n-1/2)
    allocate(Psi(np,dims,dims),Psi0(npfb,dims,dims));Psi = 0.0d0  ! Need to include memory for mirrors, so goes to np
    !$OMP PARALLEL DO PRIVATE(tmp_mat)
    do i=1,npfb
       tmp_mat(1,1) = log(EigenVals(i,1));tmp_mat(2,2) = log(EigenVals(i,2));tmp_mat(1,2)=0.0d0;tmp_mat(2,1)=0.0d0
       Psi(i,:,:) = matmul(EigenVecs(i,:,:),matmul(tmp_mat,transpose(EigenVecs(i,:,:))))
    end do  
    !$OMP END PARALLEL DO        
    Psi0(1:npfb,:,:) = Psi(1:npfb,:,:) !! 

    !! Step 4: Advect Psi using shifting velocity Psi* = Psi(n-1/2) + dt*u_shift.grad(Psi(n-1/2))
    !! Not necessary whilst LAGRANGIAN
!       if(lagrangian.eqv..false.)then
    if(.true.)then   !! Do this for Lagrangian SPH too, (correcting for shifting)
       do j=npfb+1,np
          i = irelation(j) ! Parent of j
          Psi(j,:,:) = Psi(i,:,:)    !! Psi mirror is same as parent, so dPsi/dn=0 (homogeneous Neumann)
       end do
       allocate(grad_xx(npfb,dims),grad_xy(npfb,dims),grad_yy(npfb,dims))
       call grad_operator(grad_xx(:,:),Psi(:,1,1),0)  !! Find gradients of Psi (flag 0 means don't find grad in mirrors)
       call grad_operator(grad_xy(:,:),Psi(:,1,2),0)  !! but do include mirror contributions for gradients
       call grad_operator(grad_yy(:,:),Psi(:,2,2),0)
       do i=nbS+1,npfb
          Psi(i,1,1) = Psi(i,1,1) + dt*dot_product(ushift(i,:),grad_xx(i,:))
          Psi(i,1,2) = Psi(i,1,2) + dt*dot_product(ushift(i,:),grad_xy(i,:))
          Psi(i,2,2) = Psi(i,2,2) + dt*dot_product(ushift(i,:),grad_yy(i,:))
          Psi(i,2,1) = Psi(i,1,2)   ! Psi is symmetric :)
       end do
       deallocate(grad_xx,grad_xy,grad_yy)
    end if

    !! Next few steps in a single do loop to save memory
    call random_seed()
    !$OMP PARALLEL DO PRIVATE(tmp_mat,tmp1,B,Omega,tmp_vec,A)
    do i=1,npfb
       !! Step 5: Find the elements B,Omega (at n) of the decomposed velocity gradient (at n), using Lambda, R at n-1/2
       tmp_mat = matmul(transpose(EigenVecs(i,:,:)),matmul(grad_vel(i,:,:),EigenVecs(i,:,:)))
       tmp1 = (EigenVals(i,2)*tmp_mat(1,2) + EigenVals(i,1)*tmp_mat(2,1))/(EigenVals(i,2)-EigenVals(i,1)+1d-10)
       tmp_mat(1,2) = 0.0d0;tmp_mat(2,1) = 0.0d0
       B = matmul(EigenVecs(i,:,:),matmul(tmp_mat,transpose(EigenVecs(i,:,:))))
       tmp_mat = 0.0d0;tmp_mat(1,2) = tmp1;tmp_mat(2,1) = -tmp1
       Omega = matmul(EigenVecs(i,:,:),matmul(tmp_mat,transpose(EigenVecs(i,:,:))))  

       !! Step 6: Integrate Psi* to Psi**: Psi** = Psi* + dt(2B + OmegaPsi* - Psi*Omega) (or should it be Psi(n-1/2)??)
       if(itime.eq.2.and.lagrangian.eqv..false.)then  !! We need to kick it if Eulerian...           
          call random_number(tmp1);B = B + 1.0d13*epsilon(tmp1)*(tmp1-0.5d0)
       end if
       tmp_mat = (2.0*B + matmul(Omega,Psi0(i,:,:)) - matmul(Psi0(i,:,:),Omega))
       Psi(i,:,:) = Psi(i,:,:) + dt*tmp_mat(:,:) 
             
       !! Step 7: Decompose Psi** to obtain Lambda, R at **, Psi** = Rlog(Lambda)R^T (all at **), then A** = RLambdaR^T
       call calc_eigens(Psi(i,:,:),EigenVals(i,:),EigenVecs(i,:,:)) 
       tmp_vec(:) = EigenVals(i,:)
       EigenVals(i,:) = exp(EigenVals(i,:))  !! remember to inverse the logarithm...  

       !! A** = RLambdaR^T
       tmp_mat=0.0d0;tmp_mat(1,1) = EigenVals(i,1);tmp_mat(2,2) = EigenVals(i,2)
       A(:,:) = matmul(EigenVecs(i,:,:),matmul(tmp_mat,transpose(EigenVecs(i,:,:))))

       !! Step 8: Analytic integration of relaxation function: dA/dt=-f_{R}(A**)/Wi 
       !! to obtain A at n+1/2  i.e. *-time in the current time-step.
       tmp_mat = A(:,:)
       A(:,:) = f_relaxation_integrate(tmp_mat,trace_conf(i),El*sqrt(Ra/Pr))  !! final arg is Weissenberg number
  
       !! Step 9: Determine tau_p*=f_{S}(A*)
       trace_conf(i) = A(1,1) + A(2,2)        !! Store the trace of the conformation tensor (for next step)
       tau_p(i,:,:) = f_strain_forward(A,trace_conf(i))           
a_out(i) = trace_conf(i)
    end do
    !$OMP END PARALLEL DO
    deallocate(EigenVecs,EigenVals,Psi,Psi0)    

    !! Step 10: Determine tau_p in mirror particles, for use in div.tau_p in main visc subroutine.
    call mirror_polymeric_stress

  end subroutine constitutive_eqn_evolve
!! ------------------------------------------------------------------------------------------------
  subroutine calc_eigens(Amat,Lvec,Rmat)
    implicit none
    real(rkind),dimension(:,:),intent(in) :: Amat
    real(rkind),dimension(:),intent(out) :: Lvec          
    real(rkind),dimension(:,:),intent(out) :: Rmat

    double complex,dimension(2) :: Lvec_tmp
    integer(ikind) :: l
    real(rkind) :: tr,det,nrm,tmp

    tr = Amat(1,1) + Amat(2,2)                        !! Trace and determinant
    det = Amat(1,1)*Amat(2,2) - Amat(2,1)*Amat(1,2)

    tmp = 0.25d0*tr*tr - det

    if(tmp.le.epsilon(tmp)) then    !! If machine precision will push EigenValues complex, stop it from doing so...
       tmp = 0.0d0       
    end if
    Lvec(1) = 0.5d0*tr + sqrt(tmp)    !! EigenValues     
    Lvec(2) = 0.5d0*tr - sqrt(tmp)  
!! This bit... 
!    if(Lvec(1).le.0.0d0) Lvec(1)=epsilon(Lvec(1))  !! Avoiding negative/zero eigens...
!    if(Lvec(2).le.0.0d0) Lvec(2)=epsilon(Lvec(2))
!if(Lvec(1).le.0.0d0.or.Lvec(2).le.0.0d0) write(6,*) Lvec(1),Lvec(2)
!!N.B above fix necessary to remove noise (I think). SOmething odd going on...


    if(Amat(2,1).ne.0.0d0)then     !! If A is not diagonal (remember, it is always symmetric)
       Rmat(1,1) = Lvec(1)-Amat(2,2);Rmat(2,1)=Amat(2,1)   ! Calculate EigenVec #1
       nrm=sqrt(Rmat(1,1)**2 + Rmat(2,1)**2)
       Rmat(:,1) = Rmat(:,1)/nrm                            ! normalise
       Rmat(1,2) = -Rmat(2,1); Rmat(2,2) = Rmat(1,1)     ! EigenVec #2 is orthogonal to #1
    else                          !! A is diagonal
       Lvec(1)=Amat(1,1);Lvec(2)=Amat(2,2)   
       Rmat(1,1) = 1.0;Rmat(2,2) = 1.0;Rmat(1,2) = 0.0;Rmat(2,1)=0.0
    end if

  end subroutine calc_eigens
!! ------------------------------------------------------------------------------------------------  
  function f_strain_backward(fstmp,trA) result (Atmp)     !! Subroutine to calculate the inverse of fs
     real(rkind),dimension(:,:) :: fstmp
     real(rkind) :: trA
     real(rkind),dimension(dims,dims) :: Atmp
     real(rkind) :: beta_S,alpha_S
#if const_mod==2    
     beta_S=1.0d0;alpha_S=1.0d0  !Oldroyd B
#elif const_mod==3        
     beta_S = L0*L0/(L0*L0 - trA);alpha_S=1.0d0 !Fene-P
#elif const_mod==4
     beta_S = 1.0d0;alpha_S= L0*L0/(L0*L0 - trA) !Fene-CR
#elif const_mod==5
     beta_S = 1.0d0;alpha_S = 1.0d0   !linear PTT
#elif const_mod==6
     beta_S = 1.0d0;alpha_S = 1.0d0  !Exponential PTT
#else
     beta_S=1.0d0;alpha_S=1.0d0     !Giesekus
#endif
   
     Atmp = (fstmp/alpha_S + idm)/beta_S 
  end function f_strain_backward
!! ------------------------------------------------------------------------------------------------  
  function f_strain_forward(Atmp,trA) result (fstmp)     !! Calculate tau from A
     real(rkind),dimension(dims,dims) :: fstmp
     real(rkind) :: trA
     real(rkind),dimension(:,:) :: Atmp
     real(rkind) :: beta_S,alpha_S
#if const_mod==2    
     beta_S=1.0d0;alpha_S=1.0d0  !Oldroyd B
#elif const_mod==3        
     beta_S = L0*L0/(L0*L0 - trA);alpha_S=1.0d0 !Fene-P
#elif const_mod==4
     beta_S = 1.0d0;alpha_S= L0*L0/(L0*L0 - trA) !Fene-CR
#elif const_mod==5
     beta_S = 1.0d0;alpha_S = 1.0d0   !linear PTT
#elif const_mod==6
     beta_S = 1.0d0;alpha_S = 1.0d0  !Exponential PTT
#else
     beta_S=1.0d0;alpha_S=1.0d0     !Giesekus
#endif

     fstmp = alpha_S*(beta_S*Atmp - idm)
  end function f_strain_forward
!! ------------------------------------------------------------------------------------------------  
  function f_relaxation_integrate(Atmp,trA,l) result(Ares)   !! Integration of RHS of conformation evol eqn
     real(rkind),dimension(:,:) :: Atmp
     real(rkind) :: trA,l,beta_R,alpha_R,tmp1,tmp2
     real(rkind),dimension(dims,dims) :: Ares
#if const_mod==2    
     beta_R=1.0d0;alpha_R=1.0d0
#elif const_mod==3        
     beta_R = L0*L0/(L0*L0 - trA);alpha_R=1.0d0
#elif const_mod==4
     beta_R = 1.0d0;alpha_R= L0*L0/(L0*L0 - trA)
#elif const_mod==5
     beta_R = 1.0d0;alpha_R = 1.0d0 + eps_PTT*(trA-2.0d0)
#elif const_mod==6
     beta_R = 1.0d0;alpha_R = exp(eps_PTT*(trA-2.0d0))
#endif

#if const_mod<=6
     !! Integrate analytically
     tmp1 = exp(-alpha_R*beta_R*dt/l);tmp2 = (1.0 - tmp1)/beta_R          
     Ares = Atmp*tmp1 + tmp2*idm
#else        
     tmp1 = 1.0d0 - 2.0d0*alpha_Giesekus;tmp2 = 1.0 - alpha_Giesekus
       
     !! 1st order Euler integration...
     Ares = Atmp - (dt/l)*(alpha_Giesekus*matmul(Atmp,Atmp) + tmp1*Atmp - tmp2*idm)   

     !! RK2 integration     
!     Ares = Atmp - 0.5d0*(dt/l)*(alpha_Giesekus*matmul(Atmp,Atmp) + tmp1*Atmp - tmp2*idm)        
!     Ares = Atmp - (dt/l)*(alpha_Giesekus*matmul(Ares,Ares) + tmp1*Ares - tmp2*idm)
#endif

  end function f_relaxation_integrate
!! ------------------------------------------------------------------------------------------------  
end module predictor_mod
