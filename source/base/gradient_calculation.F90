module calculation_gradient
  use kind_parameters
  use common_parameter
  use common_2d
  use mirror_boundaries
  use sphtools
  use omp_lib
  implicit none

  real(rkind) :: rad,qq
  real(rkind),dimension(dims) :: rij
  private
  public :: calc_pressure_gradient,grad_operator
contains
!! ------------------------------------------------------------------------------------------------
  subroutine calc_pressure_gradient (grad_p)
  !! Specific routine for pressure gradient
    integer(ikind) :: i,j,k
    real(rkind) :: temp
    real(rkind),dimension(dims) :: grad_p_tmp
    real(rkind), dimension(:,:), intent(inout) :: grad_p

    ! set the pressure in mirror particles using dP_mp
    !$OMP PARALLEL DO PRIVATE(j)
    do i=npfb+1,np
       j = irelation(i)
       P(i) = P(j) + dP_mp(i)    
    end do
    !$OMP END PARALLEL DO

    grad_p=0_rkind

    ! The form of the gradient approximation: Eqn 5, Lind et al.,JCP,2011
    !$OMP PARALLEL DO PRIVATE(k,j,temp,grad_p_tmp)
    do i=1,npfb
       grad_p_tmp(:) = 0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k)
          temp = (P(j)-P(i))
          grad_p_tmp(:) = grad_p_tmp(:) + temp*ij_w_G(i,k,:)
       end do
     
       grad_p(i,:) = matmul(kgcm(i,:,:),grad_p_tmp(:))

    end do
    !$OMP END PARALLEL DO

  end subroutine calc_pressure_gradient
!! ------------------------------------------------------------------------------------------------
  subroutine grad_operator(grad,pha,inc_mirrors)
  !! General gradient operator (corrected)
    integer(ikind) :: i, j, k, ib, ibp,ibm1
    real(rkind) :: temp
    integer(ikind), intent(in) :: inc_mirrors
    real(rkind), dimension(:), intent(in) :: pha
    real(rkind), dimension(:,:), intent(inout) ::grad
    real(rkind), dimension(dims) :: grad_bfr,grad_tmp
    real(rkind), dimension(dims,dims) :: Qij,QijT

    grad(:,:)=0_rkind
    !$OMP PARALLEL DO PRIVATE(k,j,temp,grad_tmp)
    do i=1,npfb
       grad_tmp(:)=0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k)
          temp = pha(j) - pha(i)
          grad_tmp(:) = grad_tmp(:) + temp*ij_w_G(i,k,:)
       end do

       grad(i,:) = matmul(kgcm(i,:,:),grad_tmp(:))   !! Apply kernel gradient correction
    end do
    !$OMP END PARALLEL DO
    

    ! Calculate the gradient for mirror particles: rotate/reflect, multiply by -1 as required
    if(inc_mirrors.ne.0)then
       !$OMP PARALLEL DO PRIVATE(i,ib,ibp,Qij,rij,grad_bfr,QijT)
       do j=npfb+1,npfb+nmirror
          i=irelation(j)
          !! First look at mirrors which are not corner mirrors
          if(vrelation(j).eq.999) then !! Mirrored in a circle
             rij(:) = r(i,:) - r(j,:) ! Calculate rotation matrix
             Qij(:,1) = rij(:)/sqrt(dot_product(rij,rij))
             Qij(1,2) = -Qij(2,1)
             Qij(2,2) = Qij(1,1)
             grad_bfr(:) = matmul(Qij(:,:),grad(i,:)) !gradu,gradv in FoR aligned with boundary
             grad_bfr(2) = -grad_bfr(2) ! Reverse tangential components ofr gradu,gradv
             QijT = transpose(Qij)                
             grad(j,:) = matmul(QijT,grad_bfr) !! Rotate gradients back to x,y frame of reference
          end if
          if(vrelation(j).gt.0.and.vrelation(j).lt.999) then 
             if(b_type(vrelation(j)).eq.2)then ! Periodic boundary, rotate gradu,gradv as required.
                ib = vrelation(j)
                ibp = b_periodic_parent(ib)
                Qij(:,:) = angle_between_patches(ib,ibp)
                grad(j,:) = matmul(Qij(:,:),grad(i,:))
             else if(b_type(vrelation(j)).eq.3) then ! Inflow boundary, copied from i (because it works?)
                grad(j,:) = grad(i,:)!0.0d0
             else if(b_type(vrelation(j)).eq.1) then ! Wall boundary, grad. tang.to wall multply by -1
                rij(:) = r(i,:) - r(j,:) ! Calculate rotation matrix
                Qij(:,1) = rij(:)/sqrt(dot_product(rij,rij))
                Qij(1,2) = -Qij(2,1)
                Qij(2,2) = Qij(1,1)
                grad_bfr(:) = matmul(Qij(:,:),grad(i,:)) !gradu,gradv in FoR aligned with boundary
                grad_bfr(2) = -grad_bfr(2) ! Reverse tangential components ofr gradu,gradv
                QijT = transpose(Qij)                
                grad(j,:) = matmul(QijT,grad_bfr) !! Rotate gradients back to x,y frame of reference
             else if(b_type(vrelation(j)).eq.4) then ! Outflow boundary. normal gradient 0, transverse copied from i
                rij(:) = r(i,:) - r(j,:) ! Calculate rotation matrix
                Qij(:,1) = rij(:)/sqrt(dot_product(rij,rij))
                Qij(1,2) = -Qij(2,1)
                Qij(2,2) = Qij(1,1)
                grad_bfr(:) = matmul(Qij(:,:),grad(i,:)) !gradu,gradv in FoR aligned with boundary
                grad_bfr(1) = 0.0d0 ! normal component is 0...
                QijT = transpose(Qij)                
                grad(j,:) = matmul(QijT,grad_bfr) !! Rotate gradients back to x,y frame of reference
             end if
          else  
             !! Next look at corner mirrors
             if(vrelation(j).lt.0.and.vrelation(j).gt.-npfb)then ! wall-wall or periodic-periodic...
                ib = -vrelation(j)
                ibm1 = mod(ib+nb_patches-2,nb_patches)+1      
                if(b_type(ib).eq.1.and.b_type(ibm1).eq.1) then   ! wall-wall corner 
                    grad(j,:) = - grad(i,:)  !! Gradients reversed
                else if(b_type(ib).eq.2.and.b_type(ibm1).eq.2) then  ! periodic,periodic corner
                    grad(j,:) = grad(i,:)    !! Gradients unchanged
                else if(b_type(ib).eq.1.and.b_type(ibm1).eq.2.or.b_type(ib).eq.2.and.b_type(ibm1).eq.1)then !wall-periodic
                   if(b_type(ib).eq.1)then
                      rij(1) = -1.0d0*b_edge(ib,2);rij(2)=b_edge(ib,1)
                   else if(b_type(ibm1).eq.1)then
                      rij(1) = -1.0d0*b_edge(ibm1,2);rij(2)=b_edge(ibm1,1)
                   end if                   
                   Qij(:,1) = rij(:)/sqrt(dot_product(rij,rij))
                   Qij(1,2) = -Qij(2,1)
                   Qij(2,2) = Qij(1,1)
                   grad_bfr(:) = matmul(Qij(:,:),grad(i,:)) !gradu,gradv in FoR aligned with boundary
                   grad_bfr(2) = -grad_bfr(2) ! Reverse tangential components ofr gradu,gradv
                   QijT = transpose(Qij)                
                   grad(j,:) = matmul(QijT,grad_bfr) !! Rotate gradients back to x,y frame of reference       
                end if
             end if
          end if
       end do
       !$OMP END PARALLEL DO
    end if
    
  end subroutine grad_operator

end module calculation_gradient
