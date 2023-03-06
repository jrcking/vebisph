!------------------------------------------------
subroutine kernel_calculation
  !------------------------------------------------
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  implicit none

  integer i,j,k

  real(rkind) rad,qq,temp,temp2
  real(rkind) rr2_tmp,eta_small,oneoverh
  real(rkind),dimension(dims) :: rij,gradw
  real(rkind),dimension(dims,dims) :: kgcf_tmp
  
#define KGC 1
  !! Calculate the kernel correction - Bonet and Lok tensor..
  !! Calculate the weights for gradients and laplacians..
  
  !! Allocation of memory is expensive, so only re-allocate if npfb or maxneighbours has increased
  maxneighbours = maxval(ij_count(1:npfb)) 

  if(npfb.gt.npfb_old)then
     deallocate(kgcm)
     allocate(kgcm(npfb,dims,dims));kgcm = 0.0d0 !! kernel gradient correction matrix
  end if  
  if(maxneighbours.gt.maxneighbours_old.or.npfb.gt.npfb_old) then
     deallocate(ij_w_L,ij_w_G)
     allocate(ij_w_L(npfb,maxneighbours));ij_w_L=0.0d0
     allocate(ij_w_G(npfb,maxneighbours,dims));ij_w_G=0.0d0
     maxneighbours_old = maxneighbours
     npfb_old = npfb
  end if
  
  
  eta_small = 0.01*h*h
  oneoverh = 1.0d0/h

  !$OMP PARALLEL DO PRIVATE(k,j,rij,rr2_tmp,rad,qq,gradw,temp,kgcf_tmp,temp2) &
  !$OMP SHARED(ij_w_G,ij_w_L)
  do i=1,npfb !loop over fluid particles
     kgcf_tmp(:,:)=0.0d0
     do k=1,ij_count(i)
        j=ij_link(i,k)  ! Particles j within interaction distance of i
        rij(:) = r(i,:) - r(j,:)       
        rr2_tmp = rij(1)*rij(1) + rij(2)*rij(2)  ! Faster than dot_product...
       
        rad=sqrt(rr2_tmp);
        qq = rad*oneoverh
        gradw(:) = rij(:)*fac(qq)/max(rad,eta2)   !uncorrected kernel gradient
        ij_w_G(i,k,:) = gradw(:)*dv
        ij_w_L(i,k) = 2.0*dot_product(rij,ij_w_G(i,k,:))/(rr2_tmp+eta_small)                
#if KGC==1
        kgcf_tmp(1,:) = kgcf_tmp(1,:) - dv*gradw(:)*rij(1)
        kgcf_tmp(2,:) = kgcf_tmp(2,:) - dv*gradw(:)*rij(2)  ! preferable to running over l=1,dims
#endif        
     end do
#if KGC==1
     !! And the inverse...
     temp = 1.0d0/((kgcf_tmp(1,1)*kgcf_tmp(2,2)-kgcf_tmp(1,2)*kgcf_tmp(2,1))+1.0d-15)
     temp2 = kgcf_tmp(1,1)
     kgcm(i,1,1) = temp*kgcf_tmp(2,2)
     kgcm(i,1,2) = -temp*kgcf_tmp(1,2)
     kgcm(i,2,1) = -temp*kgcf_tmp(2,1)
     kgcm(i,2,2) = temp*temp2
#else
     kgcm(i,1,1)=1.0d0;kgcm(i,2,2)=1.0d0;kgcm(i,1,2)=0.0d0;kgcm(i,2,1)=0.0d0
#endif
  end do
  !$OMP END PARALLEL DO
  

!----------------------------------------------------------------------
end subroutine kernel_calculation
!------------------------------------------------
