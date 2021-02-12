module nnf_models
  use kind_parameters
  use omp_lib
  implicit none

contains

  subroutine nnf_newtonian(mu,n,mu_eff)
    real(rkind),dimension(:),intent(out) :: mu_eff
    real(rkind),dimension(:),intent(in) :: mu
    integer(ikind),intent(in) :: n

    integer(ikind) :: i    

    !$OMP PARALLEL DO
    do i=1,n
       mu_eff(i) = mu(i)
    end do
    !$OMP END PARALLEL DO
    return
  end subroutine nnf_newtonian

  subroutine nnf_power_law(SR,mu,n,mu_eff)
    real(rkind),dimension(:),intent(in)::SR
    real(rkind),dimension(:),intent(out) :: mu_eff
    real(rkind),intent(in) :: mu
    integer(ikind),intent(in) :: n

    integer(ikind) :: i
    real(rkind) :: nn

    ! power law exponent: <1=shear thinning; >1=shear thickening
    nn = 1.5_rkind !0.15
    !$OMP PARALLEL DO
    do i=1,n
       mu_eff(i) = mu*(SR(i)+1.0_rkind)**(nn - 1.0)  ! power law
    end do
    !$OMP END PARALLEL DO
  end subroutine nnf_power_law
  
  subroutine nnf_bingham(SR,mu,n,mu_eff)
    ! Bi-linear model for a Bingham fluid
    real(rkind),dimension(:),intent(in)::SR
    real(rkind),dimension(:),intent(out) :: mu_eff
    real(rkind),intent(in) :: mu
    integer(ikind),intent(in) :: n

    integer(ikind) :: i
    real(rkind) :: tau_yield,alpha
    
    ! parameters of bingham fluid
    tau_yield = 1.0e-2_rkind
    alpha = 1d2
    !$OMP PARALLEL DO
    do i=1,n
       if (SR(i).le.tau_yield/(mu*alpha)) then
          mu_eff(i) = mu*alpha
       else
          mu_eff(i) = mu + tau_yield/SR(i)
       end if
    end do
    !$OMP END PARALLEL DO
  end subroutine nnf_bingham

  subroutine nnf_cross(SR,mu,n,mu_eff)
    ! cross model for shear thinning
    real(rkind),dimension(:),intent(in)::SR
    real(rkind),dimension(:),intent(out) :: mu_eff
    real(rkind),intent(in) :: mu
    integer(ikind),intent(in) :: n

    integer(ikind) :: i
    real(rkind) :: mu_inf,mu_0,alpha,kk,mm
    
    ! parameters of cross model
    alpha = 0.05_rkind
    mu_0 = mu
    mu_inf = alpha*mu_0
    kk = 1.0e-3_rkind
    mm = 3.0_rkind    
    !$OMP PARALLEL DO
    do i=1,n
       mu_eff(i) = mu_inf + (mu_0-mu_inf)/(1.0_rkind + (kk*SR(i))**mm)
    end do
    !$OMP END PARALLEL DO
  end subroutine nnf_cross

  subroutine nnf_carreau(SR,mu,n,mu_eff)
    ! carreau model for shear thinning
    real(rkind),dimension(:),intent(in)::SR
    real(rkind),dimension(:),intent(out) :: mu_eff
    real(rkind),intent(in) :: mu
    integer(ikind),intent(in) :: n

    integer(ikind) :: i
    real(rkind) :: mu_inf,mu_0,kk,mm,aa
  
    ! parameters of carreau model

    mu_0 = mu
    mu_inf = 1.0e-3_rkind
    kk = 0.14_rkind
    mm = 0.1_rkind
    aa = 1.5_rkind
    !$OMP PARALLEL DO   
    do i=1,n
       mu_eff(i) = mu_inf + (mu_0-mu_inf)*(1.0_rkind + (kk*SR(i))**aa)**((mm-1.0)/aa)
    end do
    !$OMP END PARALLEL DO
  end subroutine nnf_carreau

end module nnf_models


