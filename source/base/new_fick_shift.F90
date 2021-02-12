module part_shift
  use kind_parameters
  use common_parameter
  use common_2d
  use calculation_gradient
  use sphtools
  implicit none

  real(rkind), dimension(:,:), allocatable :: grad_C

  private 
  public :: new_fick_shift,value_correction
contains

  subroutine new_fick_shift
    !!locals
    integer(ikind) :: i
    real(rkind) :: temp, dCds,maxU,absU,shift_mag
    real(rkind), dimension(dims) :: surf_t,surf_n
    real(rkind), dimension(:), allocatable :: Diff_coeff

    !! STEP 1: Allocation and initial values, calculate conc gradient
    allocate(Diff_coeff(npfb))
    call shifting_forces

    !! Initialise shifting vector to zero
    dr_fs(:,:)=0.0d0  
    av_conc = sum(conc(nbS+1:npfb))/dble(npfb-nbS)

    ! Find the maximum particle speed
    maxU = 1.0d-10
    !$omp parallel do private(temp) reduction(max:maxU)
    do i=1,npfb
       temp = sqrt(dot_product(u(i,:),u(i,:)))
       maxU = temp
    end do
    !$omp end parallel do

    !! STEP 2: Set the shifting coefficient
    !$OMP PARALLEL DO PRIVATE(absU,temp)
    do i=1,npfb
       ! Shifting coefficient - SHOULD HAVE UNITS OF L^2
       if(i_shift_coeff.eq.0)then
          Diff_coeff(i) = 0.25d0*h*h ! If running dam break, dam break wet bed, TG vortices       0.25*h*h
       else if(i_shift_coeff.eq.1)then
          !! 4 seems to work with Newtonian viscous jet.
          !Diff_coeff(i) = 1.0*h*dt*sqrt(dot_product(u(i,:),u(i,:)))! If Poiseuille flow,jet
          absU = sqrt(dot_product(u(i,:),u(i,:)))
          temp = 0.2*exp(-absU/(0.2*maxU))  !! 0.2 for jet!!!
          Diff_coeff(i) = 4.0*h*dt*(absU + temp*maxU)  !4.0
       end if
    end do
    !$OMP END PARALLEL DO

    !! STEP 3: Calculate the shifting vector dr_fs
    shift_mag = 0.0d0
    !$OMP PARALLEL DO PRIVATE(temp,surf_n,surf_t,dCds) reduction(max:shift_mag)
    do i=1+nbS,npfb
       ! OPTION 1: A free surface particle - remove shifting normal to surface
       if(n_surf(i) .eq. 1.and.i_open_domain.eq.1) then   !! If using n_surf in closed domain, ignore for shifting...

          ! Unit surface normal and tangent
          temp=sqrt(dot_product(surf_norm(i,:),surf_norm(i,:)))
          surf_n(:)=surf_norm(i,:)/max(temp,1e-8) ! stabilise for div by 0
          surf_t(1) = -surf_n(2);surf_t(2) = surf_n(1)

          ! dC/ds concentration gradient along surface 
          dCds = dot_product(surf_t(:),grad_C(i,:))

          ! Shifting at surface, cf. Lind et al 2012, equation 27.
          dr_fs(i,:) = -Diff_coeff(i)*dCds*surf_t(:)
          !! extra surface tension??
          !dr_fs(i,:) = dr_fs(i,:) - Diff_coeff(i)*max(dot_product(surf_n(:),grad_C(i,:)),0.0d0)*surf_n(:) 

       ! OPTION 2: An ordinary internal particle, shift according to Fick's law
       else
          dr_fs(i,:) = -Diff_coeff(i)*grad_C(i,:)
       endif

       ! Set limit |dr_fs| < 0.1dx (Xu recommends, esp. for violent flows)
       temp=sqrt(dot_product(dr_fs(i,:),dr_fs(i,:)))
       shift_mag = max(shift_mag,temp)
       if(temp.gt.0.1d0*h) then
             dr_fs(i,:) = 0.1d0*h*dr_fs(i,:)/temp
       end if
    enddo
    !$OMP END PARALLEL DO

    !! STEP 4: Decallocation etc
    deallocate(Diff_coeff)
    deallocate(grad_C)
  end subroutine new_fick_shift

  subroutine value_correction   !! Not necessary if correction is done in prediction step before PPE
    ! Correct the velocities for Fickian shifting
    integer i
    real(rkind),allocatable,dimension(:,:) :: gradu,gradv

    !! STEP 1: Calculate the velocity gradients
    allocate(gradu(npfb,dims),gradv(npfb,dims))
    call grad_operator(gradu,u(:,1),0)
    call grad_operator(gradv,u(:,2),0)

    !! STEP 2: update the velocities : u' = u + grad_u.dr
    do i=1+nbS,npfb
       if(inbin(i)) cycle
       ! correct the hudrodynamic variables
       u(i,1)=u(i,1) + dot_product(gradu(i,:),dr_fs(i,:))
       u(i,2)=u(i,2) + dot_product(gradv(i,:),dr_fs(i,:))
    enddo

    !! STEP 3: Deallocation
    deallocate(gradu,gradv)

  end subroutine value_correction

  subroutine shifting_forces
    ! Subroutine calculates concentrations, concentration gradients and components of conc grads
    integer(ikind) :: i,j,k
    real(rkind) :: temp
    real(rkind) :: rad,rad_sq, qq,Wab1, fab,nt,Rt
    real(rkind), dimension(dims) :: rij
    real(rkind), allocatable, dimension(:,:) :: sum_f,sum_ft,sum_coh
    real(rkind),dimension(dims) :: sum_f_tmp,sum_ft_tmp,sum_coh_tmp,grad_C_tmp

    !! surface cohesion coefficients
    real(rkind) :: coh_f1,coh_f2

    !! Step 1: Tensile correction parameter values c.f. Eqn (23) in Lind et al. JCP (2011)
    nt=4.0d0
    Rt=0.25d0
    qq=dx/h  ! JRCK put this in (so h/dx can vary) (previously qq=1.3 was hard-coded
    Wab1=Wab(qq)

    !! Step 2: Allocation and initial zero-ing 
    allocate(grad_C(npfb,dims))
    allocate(sum_f(npfb,dims),sum_ft(npfb,dims),sum_coh(npfb,dims))

    !! Step 3: Loop over particles and then neighbours to calculate forces, concentrations
! NB. tensile correction bit is commented out now - Rt*(Wab(qq)/Wab1)**nt is the most expensive part...
    !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,fab,sum_f_tmp,sum_ft_tmp)
    do i=1,npfb
       sum_f_tmp=0.0d0;sum_ft_tmp=0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k) 
          rij(:) = r(i,:) - r(j,:)
          rad=sqrt(dot_product(rij,rij))
          qq  = rad/h
          fab=Rt*(Wab(qq)/Wab1)**nt
          ! sum_f: Sum of kernel gradients for use concentration gradient calculation
          ! see Eqns (22,24) in Lind et al, JCP 2011. 
          sum_f_tmp(:)=sum_f_tmp(:)+ij_w_G(i,k,:)
          ! Tensile correction for conc. gradient calc: see Eqn (23) in Lind et al., JCP 2011.
          sum_ft_tmp(:)=sum_ft_tmp(:)+ij_w_G(i,k,:)*fab
       end do
       sum_f(i,:) = sum_f_tmp(:);sum_ft(i,:)=sum_ft_tmp(:)
    end do
    !$OMP END PARALLEL DO

    !! Step 4: Cohesion forces as in Akinci
!    if(akinci.and.i_open_domain.eq.1) then
    if(akinci) then
       coh_f1 = 32.0d0/(pi*h)
       coh_f2 = 1.0d0 !1
       !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,temp,sum_coh_tmp)
       do i=1,npfb
          sum_coh_tmp=0.0d0
          do k=1,ij_count(i)
             j=ij_link(i,k) 
             rij(:) = r(i,:) - r(j,:)
             rad = sqrt(dot_product(rij,rij))
             qq  = rad/(2.0*dx)

             !! cohesion forces
             temp = coh_f1*coh_f2*shifting_cohesion(qq)!*1.0d0/(sum_W(i)+sum_W(j))
             sum_coh_tmp(:) = sum_coh_tmp(:) - temp*rij(:)/rad
          end do
          sum_coh(i,:) = sum_coh_tmp(:)
       end do
       !$OMP END PARALLEL DO
    end if

    ! Step 5: Final step for grad_C and contribs
    if(i_shift_coeff.eq.0)then
       !$omp parallel do
       do i=1,npfb
          grad_C(i,:) = sum_f(i,:) + sum_ft(i,:)
       end do
       !$omp end parallel do 
    else
       !$OMP PARALLEL DO
       do i=1,npfb
          ! The concentration gradient, from Lind et al
          grad_C(i,:) =sum_f(i,:)+sum_ft(i,:)!/av_conc!See Eqn 24, Lind et al, JCP,2011

          if(akinci.and.i_open_domain.eq.1)then
             if(n_surf(i).eq.0)then
                grad_C(i,:) = grad_C(i,:) + sum_coh(i,:)  !! add cohesion to interior particles
             end if
          end if
          if(akinci.and.i_open_domain.eq.0)then
             grad_C(i,:) = grad_C(i,:) + sum_coh(i,:) !! Add cohesion to all particles
          end if
       end do
       !$OMP END PARALLEL DO  
    end if
    
    deallocate(sum_f,sum_ft,sum_coh)
    
    return
  end subroutine shifting_forces
end module part_shift
