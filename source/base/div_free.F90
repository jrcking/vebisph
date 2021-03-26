module divfree
  !! This module is the heart of the code, and the main routine (div_free) performs a single complete time-step
  use kind_parameters
  use common_parameter
  use common_2d
  use calculation_gradient
  use neighbour_finding
  use mirror_boundaries
  use freesurf
  use part_shift
  use predictor_mod
  use velocity_correction
  use sphtools
  use omp_lib
#ifdef use_labfm  
  use labfm
#endif  
  implicit none

  private
  public :: div_free

contains
!! ------------------------------------------------------------------------------------------------
  subroutine div_free(x_mint,y_mint,x_maxt,y_maxt)

    !! args
    real(rkind) :: x_mint,x_maxt,y_mint,y_maxt
    !! locals
    integer(ikind) :: i
    real(rkind), allocatable, dimension(:,:) :: grad_p

    !! These flags control the ISPH algorithm
    lagrangian = .true.          !! Lagrangian or Eulerian
    external_forcing = .true.    !! Switch on external (maybe divergent) forcing term in predictor module
    fick_shifting = .true.       !! Use Fickian shifting as in Lind et al. 2012 (JCP paper)
    akinci = .false.             !! Use Akinci shifting (as in King 2020 ArXiv preprint). Only for free-surface flows
    output_everytime = .false.   !! Output every time-step
    output_mirrors = .false.     !! Include mirror/ghost particles in output

    !! STEP 1: Pass position & velocity into o variables for storage
    allocate(r0(2*npfb,dims),u0(2*npfb,dims))
    u0(1:npfb,:) = u(1:npfb,:)
    r0(1:npfb,:) = r(1:npfb,:)

    !! STEP 2: calculate the value of the time-step
    call set_time_step

    !! STEP 3: initial advection:     r* = r + dt*u
    if(lagrangian)then
       !$OMP PARALLEL DO
       do i=nbS+1, npfb   !! Don't move boundary particles...
          r(i,:) = r0(i,:) + u0(i,:)*dt
       end do
       !$OMP END PARALLEL DO
    end if

    !! STEP 4: Housekeeping, build boundary conditions, find neighbours, find interparticle weights
    !! and set free surface if required.
    call put_in_bin(x_mint,x_maxt,y_mint,y_maxt) 
    call create_mirror_particles     
    call particle_allocation(x_mint,x_maxt,y_mint,y_maxt,0) ! neighbour finding routine
    call kernel_calculation  !! Find interparticle weights
    call locate_free_surf_part !! Free surface identification and various bits


    !! STEP 5: Calculate the shifting velocity
    allocate(ushift(npar,dims))
    if(lagrangian)then
       if(fick_shifting)then
          allocate(dr_fs(npfb,dims))
          call new_fick_shift
          ushift(1:npfb,:) = dr_fs(:,:)/dt
          deallocate(dr_fs)      
       else
          ushift(1:npfb,:) = 0.0d0    
       end if
       ushift(1:nbS,:) = -u(1:nbS,:) !! Solid boundary nodes are always Eulerian (fixed in space)
    else  !! Eulerian
       ushift(1:npfb,:) = -u(1:npfb,:)  
    end if
 
    !! STEP 6: Calculate high-order weights using LABFM if desired (if unsure, ignore this)
#ifdef use_labfm
    call labfm_weights
#endif    
    
    !! STEP 7: Prediction step - Calculate viscous, polymeric and advective terms, and update velocity to u*
    call prediction_step       

    !! STEP 8: Build and solve the Pressure (correction) Poisson equation    
    call mirror_velocities    !! Calculate u* for mirror particles
    call ppe_solve_in_house_simple_dirichlet     !! Setup and solve the Poisson equation

    !! STEP 9: Calculate pressure gradient and final velocity
    allocate(grad_p(npfb,dims))
    call calc_pressure_gradient(grad_p) 

    !$OMP PARALLEL DO
    do i=nbS+1,npfb
       u(i,:) = u(i,:) + dt*(body_force(:) - grad_p(i,:))
    enddo
    !$OMP END PARALLEL DO
    deallocate(grad_p)

    !! STEP 10: Correct velocity - only required for high Reynolds flows
!    call correct_velocity ! Remove component of velocity towards the wall. 
 
    !! STEP 11: Calculate secondary properties, such as vorticity, and perform any analysis for outputting
    call secondary_properties
    
    !! STEP 12: Particle advection - move particles to final positions
    if(lagrangian)then
       !$OMP PARALLEL DO
       do i=nbS+1,npfb    !! Don't move boundary particles
          r(i,:)=r0(i,:)+0.5d0*(u(i,:)+u0(i,:))*dt + ushift(i,:)*dt 
       end do
       !$OMP END PARALLEL DO
    end if

    !! STEP 13: Housekeeping: deallocate linking lists and other things
    deallocate(irelation,dP_mp,vrelation,ij_count,ij_link,ushift)
    if(allocated(surf_norm))deallocate(surf_norm)
    deallocate(r0,u0)

  end subroutine div_free
!! ------------------------------------------------------------------------------------------------
  subroutine set_time_step
    !! subroutine to calculate the time-step as in Lee et al, JCP 227(2008)
    real(rkind) :: u_abs
    real(rkind) :: dt_visc, dt_CFL
    integer(ikind) :: i

    ! Constraint 1: CFL condition - coeff*h/max(|u|)
    umax = 1.0d-10
    do i=nbS+1,npfb
       u_abs = sqrt(dot_product(u(i,:),u(i,:)))
       if(u_abs.gt.umax)then
          umax = u_abs
       end if
    end do
    dt_CFL = dt_coef_cfl*h/umax
    ! Constraint 2: Viscous condition - coeff*min(ro*h**2/mu)
    dt_visc = dt_coef_visc*h*h*sqrt(Ra/Pr)      

    dt = min(dt_CFL,dt_visc)
    return
  end subroutine set_time_step
!! ------------------------------------------------------------------------------------------------
  subroutine secondary_properties
     integer(ikind) :: i,j,k,nf
     real(rkind) :: tmp,vort_tmp
     real(rkind) :: tmpa,tmpb

     !! Everything within this loop is done only when required for output
     if(time.gt.next_dump2.or.output_everytime) then
        next_dump2 = next_dump2 + dt_out

        !! Vorticity calculation for output...
        !$OMP PARALLEL DO PRIVATE(k,j,tmp,vort_tmp)
        do i=1,npfb
           vort_tmp=0.0d0
           do k=1,ij_count(i)
              j=ij_link(i,k)
              tmp = (u(i,2)-u(j,2))*ij_w_G(i,k,1) - (u(i,1)-u(j,1))*ij_w_G(i,k,2)
              vort_tmp = vort_tmp - tmp
           end do
           vort(i) = vort_tmp
        end do
        !$OMP END PARALLEL DO

        if(.true.)then
        !! Centre-line velocity for poiseuille flow
        k=0;tmpa=0.0d0
        !$OMP PARALLEL DO REDUCTION(+:k,tmpa)
        do i=nbS+1,npfb
           if(abs(r(i,2)-0.5d0).le.dx)then
              k=k+1
              tmpa=tmpa+u(i,1)
           end if
        end do
        !$OMP END PARALLEL DO
        tmpb=tmpa/dble(k)
        write(192,*) time,tmpb
        flush(192)
        end if
     
     end if
   
     !! Put here any calculations we want to do every time step
   
              
     return
  end subroutine secondary_properties
!! ------------------------------------------------------------------------------------------------
end module divfree
