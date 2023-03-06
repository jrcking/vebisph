module freesurf
  !! Module contains routines to locate the free surface and various other bits (e.g. concentration)
  use kind_parameters
  use common_2d
  use common_parameter
  use sphtools
  implicit none

  real(rkind),dimension(:),allocatable :: beta_div_r
  real(rkind),dimension(dims) :: rij
  private
  public :: locate_free_surf_part


contains
!! ------------------------------------------------------------------------------------------------
  subroutine locate_free_surf_part 
  ! subroutine to identify the free surface
    integer(ikind) :: i,ii

    !! find div.r, conc and surf_norm and beta_div_r
    if(i_open_domain.eq.1)then
       call calc_part_dist_vars
    end if

    !! Indentify free surface particles for open and closed domains
    if(i_open_domain.eq.1)then
       !$OMP PARALLEL DO 
       do i=1,npfb
          ! Free surface identity criterion
          beta_div_r(i)= 1.4d0 + 0.5*max(conc(i)-1.0d0,-0.2d0) 
          if(div_r(i).gt.beta_div_r(i)) then 
             n_surf(i)=0_ikind
          else
             n_surf(i)=1_ikind 
          end if
       end do
       !$OMP END PARALLEL DO
    else
       n_surf(:)=0
    end if
    !$OMP PARALLEL DO
    do i=npfb+1,np  !! set mirrors or surface to surface
       n_surf(i) = n_surf(irelation(i))
    end do
    !$OMP END PARALLEL DO

    !! smooth the surface normals
    if(i_open_domain.eq.1)then
       call smooth_normals
    end if

    !! A variable to smooth the PPE
    if(i_open_domain.eq.1)then
       if(allocated(ppe_smooth).eqv..false.) allocate(ppe_smooth(npfb))
       ppe_smooth(:) = 1.0d0
       do i =1,npfb
          if(div_r(i).gt.1.2d0*beta_div_r(i)) cycle
          if(n_surf(i).eq.1) then
             ppe_smooth(i) = 0.0d0
          else
             ppe_smooth(i) = 0.5*(1.0d0 - cos((div_r(i)-beta_div_r(i))*pi/(0.2d0*beta_div_r(i))))
          end if
       end do  
       deallocate(beta_div_r)       
    end if  

  end subroutine locate_free_surf_part
!! ------------------------------------------------------------------------------------------------
  subroutine calc_part_dist_vars
  !! subroutine to calculate variables related to particle 
  !! distribution: div.r, surf_norm, conc, beta_div_r
    integer(ikind) :: i,j,k
    real(rkind) :: av_conc,qq,rad,temp,wtmp,Aij
    real(rkind) :: conc_tmp,div_r_tmp
    real(rkind),dimension(dims) :: surf_norm_tmp
   
    !! initial allocation
    if(allocated(conc).eqv..false.) allocate(conc(npfb))
    if(allocated(beta_div_r).eqv..false.) allocate(beta_div_r(npfb))
    if(allocated(surf_norm).eqv..false.) allocate(surf_norm(npfb,dims))

    !! calculate div.r, normals and concentration 
    !$OMP PARALLEL DO PRIVATE(conc_tmp,div_r_tmp,surf_norm_tmp,k,j,rij,rad,qq,temp,wtmp)
    do i=1,npfb
       conc_tmp =0.0d0;div_r_tmp=0.0d0;surf_norm_tmp=0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k) 
          rij(:) = r(i,:) - r(j,:);rad = sqrt(dot_product(rij,rij));qq = rad/h
          wtmp = Wab(qq)
          conc_tmp = conc_tmp + wtmp
          div_r_tmp=div_r_tmp - dot_product(rij,ij_w_G(i,k,:))          ! Divergence of position vector r
          surf_norm_tmp(:)=surf_norm_tmp(:) - ij_w_G(i,k,:)         ! Surface-normal vector = sum(V*gradW)
       end do
       conc(i) = conc_tmp*dv;div_r(i)=div_r_tmp;surf_norm(i,:) = surf_norm_tmp(:)
    end do
    !$OMP END PARALLEL DO
    av_conc = sum(conc(nbS+1:npfb))/dble(npfb-nbS)
    wtmp = Wab(0.0d0)
    conc(1:npfb) = conc(1:npfb) + dv*wtmp
    surf_norm(:,:) = surf_norm(:,:)*h  ! normalise so |n|=0.504989 for planar surface

  end subroutine calc_part_dist_vars
!! ------------------------------------------------------------------------------------------------
  subroutine smooth_normals
  !! Subroutine to smooth normal vectors (via Shepard filtering)
    integer(ikind) :: i,j,k
    real(rkind) :: rad,qq,wtmp
    real(rkind),dimension(dims) :: sn_temp_tmp
    real(rkind),dimension(:,:),allocatable :: sn_temp
    allocate(sn_temp(np,dims))

    !! Calculate filtered result
    wtmp = Wab(0.0d0)    
    !$OMP PARALLEL DO PRIVATE(sn_temp_tmp,k,j,rij,rad,qq)
    do i=1,npfb
       sn_temp_tmp=0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k)
          if(j.gt.npfb) cycle ! there are no surf_norm for j in mirror
          rij(:) = r(i,:) - r(j,:)
          rad = sqrt(dot_product(rij,rij))
          qq = rad/h
          sn_temp_tmp(:) = sn_temp_tmp(:) + Wab(qq)*surf_norm(j,:)
       end do
       sn_temp(i,:) = sn_temp_tmp(:)*dv + surf_norm(i,:)*wtmp*dv
    end do
    !$OMP END PARALLEL DO
    
    !! Pass filtered result back to original array
    !$OMP PARALLEL DO
    do i=1,npfb
       surf_norm(i,:) = sn_temp(i,:)     
    end do
    !$OMP END PARALLEL DO
    deallocate(sn_temp)
  end subroutine smooth_normals 
!! ------------------------------------------------------------------------------------------------
end module freesurf
