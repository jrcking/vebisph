module mirror_boundaries
  use kind_parameters
  use common_parameter
  use common_2d

  implicit none
  real(rkind),dimension(dims) :: norm_be
  real(rkind) :: perp_dist
  integer(ikind),allocatable,dimension(:) :: n_near_patch,n_near_corner
  integer(ikind),allocatable,dimension(:,:) :: b_neighbour,bc_neighbour
  real(rkind),allocatable,dimension(:,:,:) :: b_neigh_pos
  logical :: LDC

  public :: create_mirror_particles,angle_between_patches,mirror_velocities, &
            mirror_polymeric_stress
contains
!! ------------------------------------------------------------------------------------------------
  subroutine create_mirror_particles
!! -- NB:

!!  vrelation +ve and <nb_patches  : edge
!!  vrelation +ve and >=999        : circle
!!  vrelation -ve and >-nb_patches : wall-wall or periodic-periodic corner
!!  vrelation -ve and <-nb_patches : wall-periodic
!!

!! If j is a wall mirror, vrelation(j) = ib (index of wall patch j is near)
!! If j is a periodic mirror, vrelation(j) = ib (the index of the periodic patch j is near)
!! If j is an inflow mirror, vrelation(j) = ib (index of inflow patch j is near)
!! If j is a double wall corner, vrelation(j) = -ib (the index of the boundary node (corner) it is near)
!! If j is a double periodic corner, vrelation(j) = -ib (the index of the boundary node (corner) it is near)
!! If j is a wall-periodic corner, vrelation(j) = -i  (index of mirror parent)
!! If j is a circle particle, vrelation(j) = 999
    real(rkind),dimension(2) :: trans_bp,bp,ric,tang_be
    real(rkind) tmp,tmp_angle,m_be,rad,c_dir,U_outflow
    integer(ikind) :: i,j,ib,imp,imp_temp,k,l
    integer(ikind) :: ibp,ibm1,ibb,ibp1
    real(rkind),dimension(dims,dims) :: Qij,Qijt

    ! STEP 1: Replace any particles which have escaped
    call replace_escapees

    ! STEP 2: Insert or remove particles for inflow/outflow
    !! call insert_or_remove_particles ! COMMENTED OUT AS NOT REQUIRED  
    
    ! STEP 3: Create boxes of particles close to each boundary patch (except non-boundaries)
    allocate(n_near_patch(nb_patches),n_near_corner(nb_patches))
    allocate(b_neighbour(nb_patches,npfb),bc_neighbour(nb_patches,npfb))
    allocate(b_neigh_pos(nb_patches,npfb,2))
    call create_boundary_boxes
    
    ! STEP 4: Allocate arrays for mirror-parent relations
    nmirror_esti = 2*npfb  ! adjust estimate for max number of mirrors (as npfb increases)
    if(npfb.le.500) nmirror_esti=10*npfb
    allocate(irelation(npfb+1:npfb+nmirror_esti))
    allocate(vrelation(npfb+1:npfb+nmirror_esti))
    allocate(dP_mp(npfb+1:npfb+nmirror_esti))
    imp = 0 ! Initially there are no mirror particles

    ! STEP 5a: Mirrors for the circles...
    if(nb_circles.ne.0)then
       do i=nbS+1,npfb      !! For each particle
          do ib=1,nb_circles   !! For each circle
             c_dir = c_radius(ib)/abs(c_radius(ib)) ! +1 for inner, -1 for outer
             ric(:) = r(i,:) - c_centre(ib,:)
             rad = sqrt(dot_product(ric,ric))
             if(c_dir.gt.0.0d0)then  !! regular (empty circle)
                if(rad.le.c_radius(ib)+sup_size) then
                   ric = ric/rad !! Normalise
                   imp = imp + 1
                   k = npfb + imp
                   irelation(k) = i
                   vrelation(k) = 999+ib
                   r(k,:) = r(i,:) - ric(:)*2.0*(rad-c_radius(ib))
                   tang_be(1) = ric(2);tang_be(2)=-ric(1)
                   u(k,:) = 2.0d0*c_vel(ib,:) + 2.0d0*tang_be(:)*c_radius(ib)*c_omega(ib)-1.0d0*u(i,:)   
                   dP_mp(k) = dot_product(body_force,r(k,:)-r(i,:)) 
                end if            
             else if(c_dir.lt.0.0d0)then  !! outer circle (i.e. fluid inside, boundary outside)
                if(rad.ge.abs(c_radius(ib))-sup_size) then
                   ric = ric/rad !! Normalise
                   imp = imp + 1
                   k = npfb + imp
                   irelation(k) = i
                   vrelation(k) = 999+ib
                   r(k,:) = r(i,:) + ric(:)*2.0*(abs(c_radius(ib))-rad)
                   tang_be(1) = ric(2);tang_be(2)=-ric(1)
                   u(k,:) = 2.0d0*c_vel(ib,:) + 2.0d0*tang_be(:)*abs(c_radius(ib))*c_omega(ib)-1.0d0*u(i,:)   
                   dP_mp(k) = dot_product(body_force,r(k,:)-r(i,:)) 
                end if
             end if
          end do   
       end do
    end if

    ! STEP 5: Loop over each patch, and create mirrors as required
    do ib = 1,nb_patches
       ibp1 = mod(ib,nb_patches)+1

       ! OPTION 0: A non-boundary, do nothing
       if(b_type(ib).eq.0) cycle 

       ! OPTION 1: A wall boundary
       if(b_type(ib).eq.1) then
          norm_be(:) = bound_norm(ib)
          do j=1,n_near_patch(ib) ! loop over particles near patch ib
             i=b_neighbour(ib,j)
             bp = b_neigh_pos(ib,j,:) ! the position of particle relative to boundary
             if(bp(2).le.1.0d-10) cycle  ! if wrong side of bound

             ! if the corner angle is less than 0 and a mirror of i inteferes with corner ib
             ! allow 0.1dx so mirror particles don't come too close (only for one bit of corner)
             if(b_corner_angle(ib).lt.0.0)then
                tmp_angle = 0.5*(pi + b_corner_angle(ib))
                m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
                if(bp(2)+0.1d0*dx.ge.bp(1)*tan(tmp_angle)*m_be) cycle
             end if
             ! if the corner angle is less than 0 and a mirror of i inteferes with corner ibp1
             if(b_corner_angle(ibp1).lt.0.0)then
                tmp_angle = 0.5*(pi + b_corner_angle(ibp1))
                m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
                if(bp(2).ge.(1.0d0-bp(1))*tan(tmp_angle)*m_be) cycle
             end if

             ! if still in loop, i needs a mirror k
             imp = imp + 1
             k = npfb + imp   ! index of the mirror particle
             irelation(k) = i ! parent of the mirror
             vrelation(k) = ib !mirror-parent velocity relationship

             ! position of the mirror
             r(k,:) = r(i,:) - 2.0d0*bp(2)*norm_be(:)

             ! velocity of the mirror is *-1.0
             tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
             u(k,:) = 2.0*b_vel(ib)*tang_be(:) -u(i,:)
             if(LDC)then
                tmp = b_vel(ib)*16.0d0*bp(1)*bp(1)*((1.0d0-bp(1))**2.0)*(0.5+0.5*tanh(80.0*(time-0.05)))
                u(k,:) = 2.0*tmp*tang_be(:) -u(i,:)
             end if

             ! set the pressure difference between mirror and parent
             dP_mp(k) = dot_product(body_force,r(k,:)-r(i,:)) 
          end do
       end if

       ! OPTION 2: A periodic boundary
       if(b_type(ib).eq.2)then
          ! we are looking for particles near boundary patch b_periodic_parent(ib)
          ibp = b_periodic_parent(ib)  !identify parent patch for mirror particles
          m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))    
          norm_be(:) = bound_norm(ib)
          ! loop over every particle near patch ibp
          do j=1,n_near_patch(ibp)
             i=b_neighbour(ibp,j)
             bp(:) = b_neigh_pos(ibp,j,:)
             if(bp(2).le.0.0d0) cycle ! if perp_dist negative (or 0)

             ! if still in loop, i needs a mirror k
             imp = imp + 1
             k = npfb + imp   ! index of the mirror particle
             irelation(k) = i ! parent of the mirror
             vrelation(k) = ib !mirror-parent velocity relationship
             
             ! position the mirror (relative to patches now, to account for non-parallel patches)
             trans_bp(:) = -1.0d0*bp(1)*b_edge(ib,:) - bp(2)*norm_be(:)
             r(k,:) = b_node(ib,:) + b_edge(ib,:) + trans_bp(:)

             ! velocity of the mirror is unchanged (patches must be parallel)
             u(k,:) = u(i,:)

             ! set the pressure difference between mirror and parent
             dP_mp(k) = 0.0d0
          end do
       end if
     end do       

     !! Make sure mirrors of solid particles have the same velocity as their parents...
     do j=npfb+1,npfb+imp
        i=irelation(j)
        tmp = dot_product(u(i,:),u(i,:))
        if(i.le.nbS.and.tmp.ne.0.0d0) then  !! if mirroring a moving solid particle
           u(j,:)=u(i,:)
        end if
     end do

    ! STEP 6: Create mirrors in the corners...
    do ib = 1,nb_patches !loop over the corners       
       ibm1 = mod(ib+nb_patches-2,nb_patches)+1
       
       ! OPTION 1: One of the patches is a non-boundary, do nothing
       if(b_type(ib).eq.0.or.b_type(ibm1).eq.0) cycle

       ! OPTION 2: The two patches are co-linear (cross product =0), do nothing
       tmp = b_edge(ibm1,1)*b_edge(ib,2) - b_edge(ibm1,2)*b_edge(ib,1)
       if(abs(tmp).lt.1d-10) cycle

       ! OPTION 3: a double wall corner
       if(b_type(ib).eq.1.and.b_type(ibm1).eq.1)then   ! a wall corner
          ! Check to see if the corner needs doing...
          if(b_corner_angle(ib).le.0.0) cycle

          do j=1,n_near_corner(ib) ! loop over paticles near corner
             i = bc_neighbour(ib,j)
             
             trans_bp(:) = r(i,:)-b_node(ib,:)
             perp_dist = sqrt(dot_product(trans_bp,trans_bp))

             ! create a mirror particle
             imp = imp + 1
             k = npfb + imp
             irelation(k) = i ! parent of the mirror
             vrelation(k) = -ib !mirror-parent velocity relationship

             ! position the mirror
             r(k,:) = r(i,:) - 2.0d0*trans_bp(:)

!             u(k,:) = u(i,:)              ! don't reverse the velocity !! Change with b_vel ??
             norm_be=bound_norm(ib);tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
             u(k,:) = 2.0*tmp*tang_be(:) -u(i,:)    !! Contributions from both moving bounds...
             norm_be=bound_norm(ibm1);tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1) 
             u(k,:) = u(k,:) + 2.0*tmp*tang_be(:)

             ! set the pressure difference between mirror and parent
             dP_mp(k) = dot_product(body_force,r(k,:)-r(i,:))
          end do
       end if

       ! OPTION 4: A double periodic corner 
       if(b_type(ib).eq.2.and.b_type(ibm1).eq.2)then
          ibp = b_periodic_parent(ib)  !identify parent patch for mirror particles
          do j=1,n_near_corner(ibp) ! look for particles near b_node(ibp)
             i=bc_neighbour(ibp,j)

             ! create a mirror particle
             imp = imp + 1
             k = npfb + imp
             irelation(k) = i ! parent of the mirror
             vrelation(k) = -ib !mirror-parent velocity relationship

             ! position the mirror
             trans_bp(:) = b_node(ib,:) - b_node(ibp,:)
             r(k,:) = r(i,:) + trans_bp(:)

             ! don't change the velocity
             u(k,:) = u(i,:)

             ! set the pressure difference between mirror and parent
             dP_mp(k) = 0.0d0

          end do
       end if

       ! OPTION 5: The two patches are of different type
       if(b_type(ib).ne.b_type(ibm1))then
          ! OPTION 5a. One of the two is a wall, we will look near this first...
          if(b_type(ib).eq.1.or.b_type(ibm1).eq.1)then 
             if(b_type(ibm1).eq.1)then
                ibb = ibm1   !! Ibb is the wall patch
             else
                ibb = ib
             end if
             imp_temp = imp
             do j=npfb+1,npfb + imp_temp  ! look through existing mirror particles near wall
                perp_dist = bound_dist(j,ibb)
                if(perp_dist.gt.sup_size+1.0d-10) cycle  ! cycle if it's not near the patch
                if(perp_dist.lt.1.0d-10) cycle ! cycle if it's actually aligned with the patch
                norm_be(:) = bound_norm(ibb)  ! find the normal...
                tmp = position_along_patch(j,ibb)            
                ! is it off the end (the correct end) of the patch?
                if(tmp.lt.0.0d0.and.ibb.eq.ib.or.tmp.gt.1.0d0.and.ibb.eq.ibm1) then
                   imp = imp + 1
                   k = npfb + imp
                   ! parent of the mirror is the parent of the mirror in which it is mirrored!
                   i = irelation(j)
                   irelation(k) = i
                   vrelation(k) = -ib !j !mirror-parent velocity relationship

                   ! position the mirror
                   r(k,:) = r(j,:) - 2.0d0*perp_dist*norm_be(:)

                   tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
                   u(k,:) = 2.0*b_vel(ibb)*tang_be(:)-u(j,:)          ! velocity of the mirror is *-1.0 + 2b_velocity
                   
                   ! set the pressure difference between mirror and parent
                   dP_mp(k) = dP_mp(j) ! diff between parent-mirror and parent-mirror-parent
                   dP_mp(k) = dP_mp(k) + dot_product(body_force,r(k,:)-r(j,:))
                end if
             end do             
          else
             ! OPTION 5b: Neither of the two is a wall...(TBCompleted)  
          end if
       end if
    end do
    ! STEP 7: Update np and nmirror
    nmirror = imp
    np = npfb + nmirror

    deallocate(n_near_patch,n_near_corner)
    deallocate(b_neighbour,bc_neighbour)
    deallocate(b_neigh_pos)

    return
  end subroutine create_mirror_particles
!! ------------------------------------------------------------------------------------------------
  subroutine replace_escapees
  !! Replaces any particles which have "accidentally" passed through solid boundaries
    real(rkind),dimension(2) :: trans_bp,ric
    real(rkind) tmp,m_be,rad
    integer(ikind) :: i,ib,ibp,n_escapees

    !! Circles first:
    n_escapees = 0
    if(nb_circles.ne.0)then
       do ib=1,nb_circles   !! For each circle
          do i=nbS+1,npfb      !! For each particle
             ric(:) = r(i,:) - c_centre(ib,:)
             rad = sqrt(dot_product(ric,ric)) 
             if(rad.lt.abs(c_radius(ib)).and.c_radius(ib).gt.0.0d0) then   !! if it's inside, put it back outside...
                r(i,:) = r(i,:) + (ric(:)/rad)*2.0*(c_radius(ib)-rad)
                write(6,*) "Replaced particle from inside circle",ib,i
                write(6,*) r(i,:),ric
             end if
             if(rad.gt.abs(c_radius(ib)).and.c_radius(ib).lt.0.0d0)then   !! if it's outside, put it back inside...
                r(i,:) = r(i,:) - (ric(:)/rad)*2.0*(rad + c_radius(ib))
                write(6,*) "Replaced particle from outside circle"
             end if
          end do   
       end do
    end if

    !! Then boundary patches.
    do ib = 1,nb_patches
       if(b_type(ib).eq.0) cycle  ! if it's a non-boundary, do nothing
       if(b_type(ib).eq.1.or.b_type(ib).eq.2) then   ! it's a wall or periodic
          ! loop over every particle
          norm_be(:) = bound_norm(ib)
          ibp = b_periodic_parent(ib)
          m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
          do i=nbS+1,npfb
             if(inbin(i)) cycle ! don't replace particles in bin
             ! find the perpendicular distance between i and the boundary edge
             perp_dist = bound_dist(i,ib) 
             ! if perp_dist is positive, cycle
             if(perp_dist.gt.0.0d0) cycle
             if(abs(perp_dist).gt.sup_size+1.0d-10) cycle
             tmp = position_along_patch(i,ib)
             ! if i is close to the line but beyond the ends, cycle
             ! it needs replacing, are we at a wall or a periodic?
             if(b_type(ib).eq.1) then ! WALL
                if(tmp.lt.0.0d0.or.tmp.gt.1.0d0) cycle
                ! put particle back in (reflect in boundary)
                r(i,:) = r(i,:) - 2.0d0*perp_dist*norm_be(:)
                write(6,*) "replaced particle which had passed through wall",ib,r(i,:)
                n_escapees = n_escapees + 1
             else ! PERIODIC
                if(tmp.lt.0.0d0.or.tmp.gt.1.0d0) then  !! if unexpectedly about to leave via a corner
                   trans_bp(:) = -1.0*tmp*b_edge(ibp,:) - perp_dist*bound_norm(ibp)
                   r(i,:) = b_node(ibp,:) + b_edge(ibp,:) + trans_bp(:) 
                   n_escapees = n_escapees + 1               
                else
                   ! reposition the particle relative to the boundary patches, which
                   ! allows for non-parallel patch-parents (eg. bent Poiseuille flow)
                   trans_bp(:) = -1.0*tmp*b_edge(ibp,:) - perp_dist*bound_norm(ibp)
                   r(i,:) = b_node(ibp,:) + b_edge(ibp,:) + trans_bp(:)
                   ! velocity of the mirror is unchanged
                   ! alter the pressure (only for pretty outputs)
                   n_escapees = n_escapees + 1
                end if
             end if
             ! set the o-arrays, as we're doing inflow/outflow after o-array allocation
             r0(i,:) = r(i,:) - u(i,:)*dt
             u0(i,:) = u(i,:)
          end do
       end if
    end do
  end subroutine replace_escapees
!! ------------------------------------------------------------------------------------------------
  subroutine create_boundary_boxes
    ! This subroutine finds all particles which are close to each patch and corner
    ! and makes linked lists for use in mirror particle creation
    integer(ikind) :: i,ib
    real(rkind) :: tmp

    n_near_patch(:) = 0
    n_near_corner(:) = 0
    do ib = 1,nb_patches
       if(b_type(ib).eq.0) cycle ! if it's a non-boundary, i don't care
       norm_be(:) = bound_norm(ib)
       do i=1,npfb
          ! find the perpendicular distance between i and the boundary edge
          perp_dist = bound_dist(i,ib)
          ! if perp_dist too big, cycle
          if(abs(perp_dist).gt.sup_size+1.0d-10) cycle
          ! find the normal to the boundary (pointing inwards (fluid-wards))
          tmp = position_along_patch(i,ib)
          ! if i is close to the line but beyond the ends, cycle
          if(tmp.lt.0.0d0.or.tmp.gt.1.0d0) cycle
          n_near_patch(ib) = n_near_patch(ib) + 1 ! how many in box?
          b_neighbour(ib,n_near_patch(ib)) = i   ! index of this within box
          b_neigh_pos(ib,n_near_patch(ib),1) = tmp  ! position relative to bnode
          b_neigh_pos(ib,n_near_patch(ib),2) = perp_dist ! position relative to bnode
       end do

       do i=nbS+1,npfb
          perp_dist = sqrt((r(i,1) - b_node(ib,1))**2 + (r(i,2)-b_node(ib,2))**2)
          if(perp_dist.gt.sup_size+1.0d-10) cycle
          ! it is near a corner!!
          n_near_corner(ib) = n_near_corner(ib) + 1
          bc_neighbour(ib,n_near_corner(ib)) = i
       end do
    end do

    return
  end subroutine create_boundary_boxes
!! ------------------------------------------------------------------------------------------------
  function bound_dist(i,ib) result(p_dist)
    !! Returns the perpendicular (signed) distance between particle i and patch ib
    integer(ikind),intent(in) :: i,ib
    real(rkind) :: p_dist
    integer(ikind) :: ibp1

    ! the ib+1'th patch is 1 if ib = nb_patches
    ibp1 = mod(ib,nb_patches)+1

    ! find the perpendicular distance (signed) between particle i
    ! and boundary patch ib.
    p_dist = b_edge(ib,2)*r(i,1) - b_edge(ib,1)*r(i,2) + &
         b_node(ibp1,1)*b_node(ib,2) - b_node(ibp1,2)*b_node(ib,1)
    p_dist = -1.0d0*p_dist/sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))

  end function bound_dist
!! ------------------------------------------------------------------------------------------------
  function bound_norm(ib) result(nrm)
    ! returns the unit normal vector to the boundary edge
    integer(ikind),intent(in) :: ib
    real(rkind) :: tmp1
    real(rkind),dimension(2) :: nrm

    tmp1 = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))

    nrm(1) = -1.0*b_edge(ib,2)/tmp1
    nrm(2) = b_edge(ib,1)/tmp1
  end function bound_norm
!! ------------------------------------------------------------------------------------------------
  function position_along_patch(i,ib) result(tmp)
    !! Calculates the (relative) position along patch ib or particle i
    integer(ikind),intent(in) :: i,ib
    real(rkind) :: tmp

    ! assume (it should be so) that perp_dist and norm_be are already in
    ! the correct variables
    
    ! find position along the boundary patch which is closest to particle i
    if(abs(b_edge(ib,1)).gt.abs(b_edge(ib,2)))then
       tmp = (r(i,1) - b_node(ib,1) - perp_dist*norm_be(1))/b_edge(ib,1)
    else
       tmp = (r(i,2) - b_node(ib,2) - perp_dist*norm_be(2))/b_edge(ib,2)
    end if
             
  end function position_along_patch
!! ------------------------------------------------------------------------------------------------
  function angle_between_patches(a,b) result(Rmat)
     !! Calculates the angle between two boundary patches, returns a rotation matrix
     integer(ikind),intent(in) :: a,b
     real(rkind), dimension(dims,dims) :: Rmat
     real(rkind) tmp_angle

     tmp_angle = dot_product(b_edge(a,:),b_edge(b,:))/ &
        sqrt(dot_product(b_edge(a,:),b_edge(a,:))*  &
        dot_product(b_edge(b,:),b_edge(b,:)))
     tmp_angle = acos(tmp_angle) - pi
     Rmat(1,1) =  cos(tmp_angle)
     Rmat(1,2) = -sin(tmp_angle)
     Rmat(2,1) = sin(tmp_angle)
     Rmat(2,2) = cos(tmp_angle)
   end function angle_between_patches
!! ------------------------------------------------------------------------------------------------
   subroutine mirror_velocities
   ! set mirror velocities (ie. after viscous, before ppe setup)
     integer(ikind) i,j,k,ib,ibm1
     real(rkind),dimension(dims,dims) :: Qij,Qijt
     real(rkind),dimension(dims) :: rij,tang_be
     real(rkind) :: tmp
!     !$OMP PARALLEL DO PRIVATE(i,norm_be,Qij,QijT,k,rij,ib,tang_be,tmp,ibm1)
     do j=npfb + 1,np
        i=irelation(j)
        if(vrelation(j).gt.999) then !! The mirror is across a circle
           ib = vrelation(j)-999
           u(j,:) = -u(i,:)
           rij = 0.5d0*(r(i,:)-r(j,:));tmp=sqrt(dot_product(rij,rij))
           tang_be(1)=rij(2)/tmp;tang_be(2)=-rij(1)/tmp
           u(j,:) = 2.0d0*c_vel(ib,:) + 2.0d0*tang_be(:)*c_radius(ib)*c_omega(ib)-1.0d0*u(i,:)   
        end if
        if(vrelation(j).gt.0.and.vrelation(j).lt.999) then  ! not a corner
           if(abs(b_type(vrelation(j))).eq.1)then !wall
              rij = 0.5d0*(r(i,:)-r(j,:));tmp=sqrt(dot_product(rij,rij))
              tang_be(1)=rij(2)/tmp;tang_be(2)=-rij(1)/tmp
              u(j,:) = 2.0d0*b_vel(vrelation(j))*tang_be(:) - u(i,:)
              if(LDC.and.vrelation(j).eq.3)then
                 tmp = b_vel(3)*16.0d0*r(i,1)*r(i,1)*((1.0d0-r(i,1))**2.0)*(0.5+0.5*tanh(80.0*(time-0.05)))
                 u(j,:) = 2.0*tmp*tang_be(:) -u(i,:)
              end if
           else if(b_type(vrelation(j)).eq.2) then ! periodic
              u(j,:) = u(i,:)  ! no change in velocity (assuming parallel periodic patches)
           else if(b_type(vrelation(j)).eq.3) then ! inflow, u* is U_inflow

           else if(b_type(vrelation(j)).eq.4)then ! outflow
              !! Constant velocity across the patch...
           end if
        end if
        if(vrelation(j).lt.0.and.vrelation(j).gt.-npfb)then ! wall-wall or periodic-periodic or wall-periodic
           ib = -vrelation(j)
           ibm1 = mod(ib+nb_patches-2,nb_patches)+1      
           if(b_type(ib).eq.1.and.b_type(ibm1).eq.1) then   ! wall-wall corner
              u(j,:) = u(i,:)  ! don't reverse velocity
              if(b_vel(ib).ne.0.0d0)then 
                 norm_be=bound_norm(ib);tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
                 tmp = b_vel(ib)
                 u(j,:) = 2.0*tmp*tang_be(:) -u(i,:)    !! Contributions from both moving bounds...
              end if
              if (b_vel(ibm1).ne.0.0d0)then                
                 norm_be=bound_norm(ibm1);tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1) 
                 tmp = b_vel(ibm1)
                 u(j,:) = u(j,:) + 2.0*tmp*tang_be(:)
              end if
           else if(b_type(ib).eq.2.and.b_type(ibm1).eq.2) then  ! periodic,periodic corner
              u(j,:) = u(i,:)   ! assuming parallel patches
           else if(b_type(ib).eq.1.and.b_type(ibm1).eq.2.or.b_type(ib).eq.2.and.b_type(ibm1).eq.1)then ! wall-periodic corner
              if(b_vel(ib).ne.0.0d0)then
                 norm_be=bound_norm(ib);tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
                 tmp = b_vel(ib)                 
              else if(b_vel(ibm1).ne.0.0d0)then
                 norm_be=bound_norm(ibm1);tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1) 
                 tmp = b_vel(ibm1)             
              else
                 tmp=0.0d0;tang_be=0.0d0
              end if
              u(j,:) = 2.0d0*tmp*tang_be(:) - u(i,:)
           end if
        end if
     end do
!     !$OMP END PARALLEL DO

     !! Make sure mirrors of solid particles have the same velocity as their parents...
     do j=npfb+1,np
        i=irelation(j)
        tmp = dot_product(u(i,:),u(i,:))
        if(i.le.nbS.and.tmp.ne.0.0d0) then  !! if mirroring a moving solid particle
           u(j,:)=u(i,:)
        end if
     end do
   return
   end subroutine mirror_velocities
!! ------------------------------------------------------------------------------------------------
   subroutine mirror_polymeric_stress
     integer(ikind) i,j,k,ib,ibm1
     real(rkind),dimension(dims) :: n_n,rij,xv,yv,n_t
     real(rkind),dimension(dims,dims) :: Qij,QijT,tmp
     real(rkind) :: T11,T12,T22,tmp0,m_be

     xv = (/1,0/);yv = (/0,1/)
     !$OMP PARALLEL DO PRIVATE(i,rij,n_n,n_t,Qij,QijT,tmp,tmp0,m_be,T11,T12,T22,k)
     do j=npfb + 1,np ! look through mirrors
        i=irelation(j)
        if(vrelation(j).gt.999) then !! Mirrorred in a circle
           ib = vrelation(j)-999
           rij(:) = r(i,:) - r(j,:)
           n_n = rij(:)/sqrt(dot_product(rij,rij))   ! wall normal
           n_t(1) = -n_n(2);n_t(2) = n_n(1)
           Qij(1,1) = dot_product(xv,n_n);Qij(1,2) = dot_product(xv,n_t)
           Qij(2,1) = dot_product(yv,n_n);Qij(2,2) = dot_product(yv,n_t)
           QijT = transpose(Qij)
           tmp = matmul(QijT,matmul(tau_p(i,:,:),Qij))    ! rotate into frame of ref w/ wall
           tmp(1,:) = -tmp(1,:)                       ! tau_p.n=0 on the wall (tau_p.n in mirror = -tau_p.n in fluid)
           tau_p(j,:,:) = matmul(Qij,matmul(tmp,QijT))
        end if

        if(vrelation(j).gt.0.and.vrelation(j).lt.999) then  ! not a corner
           if(b_type(vrelation(j)).eq.1)then !wall
              rij(:) = r(i,:) - r(j,:)
              n_n = rij(:)/sqrt(dot_product(rij,rij))   ! wall normal
              n_t(1) = -n_n(2);n_t(2) = n_n(1)
              Qij(1,1) = dot_product(xv,n_n);Qij(1,2) = dot_product(xv,n_t)
              Qij(2,1) = dot_product(yv,n_n);Qij(2,2) = dot_product(yv,n_t)
              QijT = transpose(Qij)
              tmp = matmul(QijT,matmul(tau_p(i,:,:),Qij))    ! rotate into frame of ref w/ wall
              tmp(1,:) = -tmp(1,:)                       ! tau_p.n=0 on the wall (tau_p.n in mirror = -tau_p.n in fluid)
              tau_p(j,:,:) = matmul(Qij,matmul(tmp,QijT))
           else if(b_type(vrelation(j)).eq.2) then ! periodic
              tau_p(j,:,:) = tau_p(i,:,:)
           else if(b_type(vrelation(j)).eq.3) then ! 
              tau_p(j,:,:) = 0.0d0
            else if(b_type(vrelation(j)).eq.4) then ! outflow - don't change tau_p
              tau_p(j,:,:) = tau_p(i,:,:)
           end if      
        end if
        if(vrelation(j).lt.0.and.vrelation(j).gt.-npfb)then ! wall-wall or periodic-periodic...
           ib = -vrelation(j)
           ibm1 = mod(ib+nb_patches-2,nb_patches)+1    
           if(b_type(ib).eq.1.and.b_type(ibm1).eq.1) then   ! wall-wall corner
              tau_p(j,:,:) = -tau_p(i,:,:) ! 0.0d0  ! assuming right-angle corner (close enough...)
           else if(b_type(ib).eq.2.and.b_type(ibm1).eq.2) then  ! periodic,periodic corner
              tau_p(j,:,:) = tau_p(i,:,:)
           else if(b_type(ib).eq.1.and.b_type(ibm1).eq.2.or.b_type(ib).eq.2.and.b_type(ibm1).eq.1)then !wall-periodic
              if(b_type(ib).eq.1) then       !! Find bound-norm of wall
                 norm_be(:) = bound_norm(ib)
              else if(b_type(ibm1).eq.1)then
                 norm_be(:) = bound_norm(ibm1)
              end if
              n_n=norm_be(:);n_t(1)=-n_n(2);n_t(2)=n_n(1)  !! Rotation matrix
              Qij(1,1) = dot_product(xv,n_n);Qij(1,2) = dot_product(xv,n_t)
              Qij(2,1) = dot_product(yv,n_n);Qij(2,2) = dot_product(yv,n_t)
              QijT = transpose(Qij)              
              tmp = matmul(QijT,matmul(tau_p(i,:,:),Qij))
              tmp(1,:) = -tmp(1,:)
              tau_p(j,:,:) = matmul(Qij,matmul(tmp,QijT))              
           end if
        end if
     end do
     !$OMP END PARALLEL DO
   return
   end subroutine mirror_polymeric_stress
!! ------------------------------------------------------------------------------------------------
end module mirror_boundaries
