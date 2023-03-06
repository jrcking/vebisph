module input
  use common_parameter
  use common_2d
  use cell_griddin
  implicit none

  private 
  public :: getdata

contains
!! ------------------------------------------------------------------------------------------------
  subroutine getdata()
    real(rkind) :: re_time,p_dist,error_mag
    real(rkind),dimension(dims) :: tmp
    integer(ikind) :: i,j,ib

    !! read in data on run parameters
    open(11,file='INDAT')
    read(11,*) i_RESTART
    read(11,*) i_open_domain
    read(11,*) tmax
    read(11,*) i_shift_coeff
    read(11,*) np              
    read(11,*) npfb                !! Number of particles (fluid + boundaries)
    read(11,*) nbS                 !! Number of boundary particles
    read(11,*) dt_out
    read(11,*) h                   !! Smoothing length 
    read(11,*) dx
    read(11,*) dt
    read(11,*) body_force(1),body_force(2)
    read(11,*) Ra
    read(11,*) Pr
    read(11,*) El
    read(11,*) beta
    read(11,*) dt_coef_cfl,dt_coef_visc
    read(11,*) ve_nonlinearity
    read(11,*) alpha_evss    
    close(11)

    itime = 0
    time = 0.0d0

    !! estimate the number of mirror particles, set some initial (small) sizes and allocate ij weights
    nmirror_esti=int(0.5*npfb)
    if(npfb.le.500) nmirror_esti = 10*npfb
    npfb_old = 1
    maxneighbours_old = 1
    allocate(ij_w_G(1,1,1),ij_w_L(1,1),kgcm(1,1,1))

    !! read in data on boundaries JRCK
    open(10,file='IBOUND')
    read(10,*) nb_patches
    allocate(b_node(nb_patches,2),b_edge(nb_patches,2),b_vel(nb_patches),b_thermal(nb_patches))
    allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
    do i=1,nb_patches
       read(10,*) b_node(i,:)
       read(10,*) b_edge(i,:)
       read(10,*) b_type(i)
       read(10,*) b_vel(i)   !! transverse velocity on (wall) patch, in direction of b_edge...
       read(10,*) b_thermal(i)
       if(b_type(i).eq.2)then !period boundary patch
          read(10,*) b_periodic_parent(i)
       end if
       if(b_type(i).eq.3)then ! inflow patch
          read(10,*) U_inflow
          U_inflow0=U_inflow
       end if      
    end do
    !! Circles !!
    read(10,*) nb_circles
    if(nb_circles.ne.0)then 
       allocate(c_centre(4*nb_circles,2),c_radius(4*nb_circles),c_omega(4*nb_circles),c_vel(4*nb_circles,2),c_thermal(4*nb_circles))
       do i=1,nb_circles
          read(10,*) c_centre(i,:)
          read(10,*) c_radius(i)
          read(10,*) c_omega(i)
          read(10,*) c_vel(i,:)
          read(10,*) c_thermal(i)
       end do
    end if
    close(10)
    write(6,*) "Boundaries read in. There are ",nb_patches," boundary patches"

    
    !! check boundaries JRCK
    call check_boundaries

    !! create a box around the closed boundaries JRCK
    xb_min = minval(b_node(:,1))
    xb_max = maxval(b_node(:,1))
    yb_min = minval(b_node(:,2))
    yb_max = maxval(b_node(:,2))
    write(6,*) "xb_min, xb_max=",xb_min,xb_max
    write(6,*) "yb_min, yb_max=",yb_min,yb_max


    print*,'xb_min,xb_max',xb_min, xb_max
    print*,'yb_min,yb_max',yb_min, yb_max
    print*,'np   ',np
    print*,'itime',itime

    time=dble(itime)*dt  ! need check if it is useful here

    !! Read in particle data
    allocate(r(npar,dims))
    allocate(u(npar,dims))
    allocate(p(npar),conc(npar),a_out(npar),theta(npar))
    allocate(vort(npar))
#if const_mod!=1
       allocate(tau_p(npar,dims,dims)) !only allocate if VE fluid
       allocate(trace_conf(npar));trace_conf = 2.0
#endif
    allocate(n_surf(npar),div_r(npar))
    allocate(inbin(npar))
    open(13,file='IPART')
    open(14,file='TEMP_REC.dat')

    if(i_RESTART) then
       read(14,*)itime,re_time,n_dump
       read(14,*)np,npfb,nbS,np_inflow
       read(14,*)xb_min,xb_max,yb_min,yb_max, dx, av_conc
       time=re_time
       print*,'CASE RESTART FROM',itime,'TH ITERATION'
       print*,'RESTART TIME',re_time
       do i=1,npfb
          read(14,*)r(i,1), r(i,2),u(i,1),u(i,2), p(i), &
               dv, n_surf(i)
#if const_mod!=1
             read(14,*) tau_p(i,1,1),tau_p(i,1,2),tau_p(i,2,1),tau_p(i,2,2)
#endif
       enddo
       if(np_inflow.ne.0)then
          allocate(p_inflow(np_inflow+2))
          do i=1,np_inflow
             read(14,*) p_inflow(i)
          end do
       end if
    else
       do i=1,npfb
          read(13,*)r(i,1), r(i,2), u(i,1), u(i,2), theta(i), dv
#if const_mod!=1
          tau_p(i,:,:)=0.0d0   !! Initialise polymeric stress to zero
#endif
       enddo
    end if

    !! Initialise all fluid particles to free surface if open domain
    n_surf(:) = 0_ikind
    if(i_open_domain.eq.0)then
       n_surf(1:npfb) = 0_ikind
    else
       n_surf(1:npfb) = 1_ikind
    end if
    close(13)
   

    !! Setup which particles are inflow particles...
    if(i_RESTART.eqv..false.)then
       call set_inflow_particles
       n_free = 0
       allocate(p_free(npar))
    end if

    h2 = h*h  
    sup_size     = 3.0d0*h
    sup_size_2 = 9.0d0*h2!sup_size*sup_size
    uno_sup_size = 1.0d0/(3.0d0*h)
    ad_7=7.0d0/478.0d0/pi/h/h
    ad_7h=ad_7/h

    !! eta2
    eta2  = h*h*1.d-10 ! Eta2 appears in the denominator of many SPH approximations

    !! initially, all particles are in the domain. If they leave for
    !! some reason, they will be placed 'in the bin'
    inbin(:) = .false.

    close(14)
    return
  end subroutine getdata
!! ------------------------------------------------------------------------------------------------
  subroutine check_boundaries
    !! This subroutine checks the boundary conditions read in from
    !! IBOUND for Jack's boundary condition framework. JRCK
    real(rkind),dimension(dims) :: tmp
    real(rkind) :: error_mag,tmp1,tmp2,p_dist
    integer(ikind) :: i,j,n_fail,ibm1,ibp1
    logical :: b_check_pass

    ! check 1: node(i) + edge(i) = node(i+1) ??
    b_check_pass = .true.
    n_fail = 0
    do i=1,nb_patches
       if(i.eq.nb_patches) then
          j=1
       else
          j=i+1
       end if
       tmp(:) = b_node(i,:) + b_edge(i,:) - b_node(j,:)
       error_mag = sqrt(dot_product(tmp,tmp))
       if(error_mag.gt.1.0d-10)then
          n_fail = n_fail + 1
          write(6,*) "Check 1 failed: b_nodes",i,"and ",j,"not connected"
          b_check_pass = .false.
       end if
    end do
    if(b_check_pass)then
       write(6,*) "Check 1 passed: all nodes and edges connected within tolerance"
    end if

    !check 2: sum(edges) = 0 ??
    tmp(:) = 0.0
    do i=1,nb_patches
       tmp(:) = tmp(:) + b_edge(i,:)
       write(6,*) tmp(:)
    end do
    error_mag = sqrt(dot_product(tmp,tmp))
    if(error_mag.gt.1.0d-10)then
       n_fail = n_fail + 1
       write(6,*) "Check 2 failed: edges sum to",tmp
    else
       write(6,*) "Check 2 passed: edges sum to zero"
    end if

    ! checks 3: periodic patches matching ??
    b_check_pass = .true.
    do i=1,nb_patches
       if(b_type(i).eq.2) then !periodic
          tmp = b_edge(i,:) + b_edge(b_periodic_parent(i),:)
          error_mag = sqrt(dot_product(tmp,tmp))
          if(abs(error_mag).gt.1.0d-10)then
             n_fail = n_fail + 1
             write(6,*) "Check 3 failed: periodic patches ",i,b_periodic_parent(i),"different length/direction"
             write(6,*) error_mag
             b_check_pass = .false.
          end if
       end if
    end do
    if(b_check_pass)then
       write(6,*) "Check 3 passed: all (if any) periodic patches well matched"
    end if
    if(n_fail.ne.0)then
       stop
    end if

    ! check 4: calculate angles at corners...
    allocate(b_corner_angle(nb_patches))
    error_mag = 0.0
    write(6,*) "For information, the corner angles are:"
    do i=1,nb_patches  ! for each corner
       ibm1 = mod(i+nb_patches-2,nb_patches)+1
       ibp1 = mod(i,nb_patches)+1
       ! find the perpendicular distance (signed) between edge ibm1 and node ibp1
       p_dist = b_edge(ibm1,2)*b_node(ibp1,1) - b_edge(ibm1,1)*b_node(ibp1,2) + &
         b_node(i,1)*b_node(ibm1,2) - b_node(i,2)*b_node(ibm1,1)
       p_dist = -1.0d0*p_dist/sqrt(dot_product(b_edge(ibm1,:),b_edge(ibm1,:)))
       ! find the magnitude of edge i
       tmp1 = sqrt(dot_product(b_edge(i,:),b_edge(i,:)))
       ! find the angle
       if(p_dist.le.0.0d0.and.abs(p_dist).gt.abs(tmp1))then
          p_dist = -tmp1
       else if(p_dist.gt.0.0d0.and.abs(p_dist).gt.abs(tmp1))then
          p_dist = tmp1
       end if
!       p_dist = min(p_dist,tmp1) ! if p_dist=1.0000000002(eg) tmp2=NaN, hence this line
       tmp2 = asin(p_dist/tmp1)
       b_corner_angle(i) = tmp2
       write(6,*) ibm1,i,tmp2
       error_mag = error_mag + tmp2
    end do
    write(6,*) "Check 4: Angles sum to : ",error_mag
    return
  end subroutine check_boundaries
!! ------------------------------------------------------------------------------------------------
  subroutine set_inflow_particles

    integer(ikind) :: i,j,ib,ibp1,np_inflow_esti
    real(rkind) :: p_dist,tmp,alph_inflow
    real(rkind),dimension(2) :: nrm 
    integer(ikind),dimension(:),allocatable :: p_inflow_tmp
    real(rkind),dimension(:),allocatable :: p_inflow_dist,p_inflow_dist_tmp
    logical,dimension(:),allocatable :: p_inflow_mask

    ! I assume there is only 1 inflow patch - if >1, need to make changes
    ! I also assume that the order I find particles will be in the order they appear
    ! along the boundary patch... (have to be careful here)
    np_inflow = 0
    alph_inflow = 1.0d0*dx  ! this is for cartesian particle spacing
    do ib = 1,nb_patches
       if(b_type(ib).ne.3) cycle ! Only for the inflow patches
       
       write(6,*) "WARNING. INFLOW NOT CODED IN THIS VERSION. STOPPING"
       stop
       
       ibp1 = mod(ib,nb_patches)+1
       ! calculate the boundary normal
       tmp = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
       nrm(1) = -1.0*b_edge(ib,2)/tmp
       nrm(2) = b_edge(ib,1)/tmp

       ! allocate p_inflow
       np_inflow_esti = 2*int(tmp/dx) + 2! *2 for safety!!
       allocate(p_inflow_tmp(np_inflow_esti));p_inflow_tmp(:)=0
       allocate(p_inflow_dist_tmp(np_inflow_esti));p_inflow_dist_tmp(:)=0.0d0
       
       do i=1,npfb ! look through all particles
          ! find the perpendicular distance (signed) between particle i
          ! and boundary patch ib.
          p_dist = b_edge(ib,2)*r(i,1) - b_edge(ib,1)*r(i,2) + &
               b_node(ibp1,1)*b_node(ib,2) - b_node(ibp1,2)*b_node(ib,1)
          p_dist = -1.0d0*p_dist/sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
          if(p_dist.lt.0.0) cycle ! discard if i is the wrong side of ib
          if(p_dist.gt.alph_inflow) cycle ! discard if it is too far from ib
          ! find position along the boundary patch which is closest to particle i
          if(abs(b_edge(ib,1)).gt.abs(b_edge(ib,2)))then
             tmp = (r(i,1) - b_node(ib,1) - p_dist*nrm(1))/b_edge(ib,1)
          else
             tmp = (r(i,2) - b_node(ib,2) - p_dist*nrm(2))/b_edge(ib,2)
          end if
          if(tmp.lt.0.0d0.or.tmp.gt.1.0d0) cycle  ! discard if it is not between ends of ib
          if(i.le.nbS.and.p_dist.gt.1.0d-10) cycle  ! if boundary, and not on the corner/ b_node, cycle

          ! if still here, it should be labelled an inflow particle!!
          np_inflow = np_inflow + 1
          p_inflow_tmp(np_inflow) = i
          p_inflow_dist_tmp(np_inflow) = tmp
       end do

       !! Look through inflow particles, and order along patch...
       allocate(p_inflow(np_inflow),p_inflow_dist(np_inflow));p_inflow_dist(1:np_inflow)=p_inflow_dist_tmp(1:np_inflow)
       allocate(p_inflow_mask(np_inflow));p_inflow_mask(:)=.true.
       do i=1,np_inflow
          j=minloc(p_inflow_dist,1,p_inflow_mask)  !! index of *next* nearest particle
          p_inflow_mask(j)=.false.  !! set j's mask to false, to exclude from next search
          p_inflow(i)=p_inflow_tmp(j)   
       end do       
       deallocate(p_inflow_tmp,p_inflow_dist)

    end do
    return
  end subroutine set_inflow_particles
!! ------------------------------------------------------------------------------------------------
end module input
