program datgen
  use kind_parameters
  use common_parameter
  use global_variables 
  implicit none

  real(rkind) :: x,y,rho1,v_mu1
  
  !! Dimensionless groups:
  real(rkind) :: beta,Pr,Ra,El   !! Primary
  real(rkind) :: Re,Wi           !! Secondary
  real(rkind) :: ve_nonlinearity,alpha_evss    
  
  real(rkind),dimension(2) :: body_force

  logical :: i_RESTART
  real(rkind) :: dt_coef_cfl,dt_coef_acc,dt_coef_visc

  integer ipart,itest,shift_coeff
  integer i,j,icha,nn
  double precision h0,tmax

  double precision yh,lambda
  double precision Cr,grv,u,v,pp,yl
  double precision thta


  body_force = 0.0d0   
  nb_circles=0


  write(*,*) 'Cases: '
  write(*,*) '  case  1:  2D dam break'
  write(*,*) '  case  2:  Taylor Green vortices'
  write(*,*) '  case  3:  Poiseuille flow' 
  write(*,*) '  case  4:  Rayleigh-Benard flow'
  write(*,*) '  case  5:  Grilli flow'
  write(*,*) '  case  6:  Free-surface Rayleigh-Benard'
  write(*,*) '  case  7:  Plane Couette flow'

  write(*,*) '  '
  write(*,*) 'Input test case number: '
  read(*,*) itest

  select case (itest) 
     !___________________________________________________________________________
     !   
  case(1)

     !  2 d dam break
     !  ooooooooooooooooooooooooooooooo
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  oxxxxxxxxxxxx                 o
     !  oxxxxxxxxxxxx                 o
     !  oxxxLIQUIDxxx                 o
     !  oxxxxxxxxxxxx                 o
     !  oxxxxxxxxxxxx                 o
     !  ooooooooooooooooooooooooooooooo 

     !! Domain size:
     xl=1.0d0;yl=1.0d0   !! Should be of order unity, as simulation is dimensionless
     h0=0.25d0

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 0 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=1
     dt_out = 0.01;tmax = 100.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/0.0d0,-1.0d0/)    
     ve_nonlinearity=100d0;alpha_evss=0.0d0
     beta = 1.0d0;Ra=1.0d4;Pr=1.0d0;El=0.0d0

     !! Resolution:
     dx=h0/40.0;h=1.3d0*dx      !! N.B. you should only need to change dx

     !   Build domain
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches),b_thermal(nb_patches))
     b_type(:) = (/ 1, 1, 1, 1 /) 
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_node(1,:) = (/ 0.0d0, 0.0d0 /) 
     b_node(2,:) = (/ xl, 0.0d0 /)
     b_node(3,:) = (/ xl, yl /)
     b_node(4,:) = (/ 0.0d0, yl /)
     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     !! outsource making wall particles... 
     call make_boundary_particles
     ipart = nbS

      ! -- fluid particles-- 
     x = xb_min + dx*0.5
     do while(x.lt.xb_max-0.5*dx)
        y = yb_min + dx*0.5
        do while(y.lt.yb_max - 0.5*dx)
           if(x.le.h0.and.y.le.2.0*h0)then
              ipart=ipart+1
              xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0
           end if
           y = y + dx
        enddo
        x = x + dx
     enddo
     
     npfb = ipart
     np=ipart

  !! ----------------------------------------------------------------------------------------------
  case(2)
     !      Taylor-Green Vortices
     !
     !     ---------------------------
     !    |   |    ^        |    ^   |
     !    |   v    |        v    |   |
     !    |<--      <------       <--| 
     !    |                          |
     !    |-->       ------>      -->|
     !    |   |    ^        |    ^   |
     !    |   |    |        |    |   |
     !    |   |    |        |    |   |
     !    |   v    |        v    |   |
     !    |<--      <------       <--|
     !    |                          |
     !    |-->       ------>      -->|
     !    |   |    ^        |    ^   |
     !    |   v    |        v    |   |
     !     --------------------------
     
     !! Domain size:
     xl=1.0d0;yl=1.0d0   !! Should be of order unity, as simulation is dimensionless

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 0 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=0
     dt_out = 0.01;tmax = 100.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/0.0d0,0.0d0/)    
     ve_nonlinearity=100d0;alpha_evss=0.0d0
     beta = 1.0d0;Ra=1.0d4;Pr=1.0d0;El=0.0d0

     !! Resolution:
     dx=yl/100.0;h=1.3d0*dx      !! N.B. you should only need to change dx


     ! Build domain
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2),b_vel(nb_patches),b_thermal(nb_patches))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))     
     b_type(:) = (/ 2, 2, 2, 2 /)  ! all periodic boundaries
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_vel(:) = (/ -0.0d0,0.0d0,0.0d0,0.0d0/)
     b_thermal(:) = (/ 1.0d0,0.0d0,-1.0d0,0.0d0/)
     b_node(1,:) = (/ 0d0, 0.0d0 /)
     b_node(2,:) = (/ xl, 0.0d0 /)
     b_node(3,:) = (/ xl, yl /)
     b_node(4,:) = (/ 0.0d0, yl /)
     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1))
     xb_max = maxval(b_node(:,1))
     yb_min = minval(b_node(:,2))
     yb_max = maxval(b_node(:,2))

     !! outsource making wall particles... JRCK
     call make_boundary_particles
     ipart = nbS

     ! -- fluid particles--    
     y=yb_min+dx/2.0d0

     do while (y .lt. yb_max)
        x=xb_min+dx/2.0d0
        do while(x .lt. xb_max)
           ipart=ipart+1
           xp(ipart)=x;yp(ipart)=y;theta(ipart)=0.0d0;
           x=x+dx
        enddo
        y=y+dx
     enddo
     npfb=ipart
     np = ipart  

     do i=1,npfb
        up(i)=-cos(2.0d0*xp(i)*pi)*sin(2.0d0*yp(i)*pi)
        vp(i)=sin(2.0d0*xp(i)*pi)*cos(2.0d0*yp(i)*pi)
     enddo


  !! ----------------------------------------------------------------------------------------------
  case(3)

     !! POISEUILLE FLOW

     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! |                                        |
     ! |                                        |
     ! |          ----> pressure gradient       |
     ! |                                        |
     ! |                                        |
     ! |                periodic bounds in x    |
     ! |                                        |  
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo

     !! Domain size:
     xl=1.0d0;yl=1.0d0   !! Should be of order unity, as simulation is dimensionless

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 0 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=0
     dt_out = 0.1;tmax = 100.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/8.0d0,0.0d0/)    
     ve_nonlinearity=0.0d0;alpha_evss=0.0d0
     beta = 0.1d0;Ra=1.0d0;Pr=1.0d0;El=1.0d0

     !! Resolution:
     dx=yl/30.0;h=1.3d0*dx      !! N.B. you should only need to change dx
    
     !! Build domain
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2),b_vel(nb_patches),b_thermal(nb_patches))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 1, 2/)
     b_vel(:) = (/ -0.0d0,0.0d0,0.0d0,0.0d0/)
     b_thermal(:) = (/ 0.0d0,0.0d0,-0.0d0,0.0d0/)
     b_periodic_parent(:) = (/ 3, 4, 1, 2/)
     b_node(1,:) = (/ 0.0d0, 0.0d0 /)
     b_node(2,:) = (/ xl, 0.0d0 /)
     b_node(3,:) = (/ xl, yl /)
     b_node(4,:) = (/ 0.0d0, yl /)
     nb_circles = 0
   
     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     call make_boundary_particles
     ipart = nbS

      ! -- fluid particles-- 
     x = xb_min + 0.5*dx
     do while(x.lt.xb_max)
        y = yb_min + dx
        do while(y.lt.yb_max - 0.5*dx)
           ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0;theta(ipart)=0.0d0
           if(y.le.0.5d0*yl + 0.25d0*yl*sin(2.0*pi*x))then
              theta(ipart) = 0.0d0
           else
              theta(ipart)=-0.0d0
           end if
           y = y + dx
        enddo
        x = x + dx
     enddo
    
     npfb = ipart
     np=ipart

!! ------------------------------------------------------------------------------------------------
  case(4)

     !! Rayleigh-Benard flow

     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! |                                        |
     ! |                                        |
     ! |                                        |
     ! |                                        |
     ! |                                        |
     ! |                periodic bounds in x    |
     ! |                                        |  
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo

     !! Domain size:
     xl=1.0d0;yl=1.0d0   !! Should be of order unity, as simulation is dimensionless

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 0 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=0
     dt_out = 0.1;tmax = 100.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/0.0d0,0.0d0/)    
     ve_nonlinearity=0.1d0;alpha_evss=0.0d0
     beta = 0.9d0;Ra=1.0d4;Pr=1.0d0;El=1.0d-4

     !! Resolution:
     dx=yl/100.0;h=1.3d0*dx      !! N.B. you should only need to change dx
    
     !! Build domain
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2),b_vel(nb_patches),b_thermal(nb_patches))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 1, 2/)
     b_vel(:) = (/ -0.0d0,0.0d0,0.0d0,0.0d0/)
     b_thermal(:) = (/ 1.0d0,0.0d0,-1.0d0,0.0d0/)
     b_periodic_parent(:) = (/ 3, 4, 1, 2/)
     b_node(1,:) = (/ 0.0d0, 0.0d0 /)
     b_node(2,:) = (/ xl, 0.0d0 /)
     b_node(3,:) = (/ xl, yl /)
     b_node(4,:) = (/ 0.0d0, yl /)
     nb_circles = 0
   
     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     call make_boundary_particles
     ipart = nbS

      ! -- fluid particles-- 
     x = xb_min + 0.5*dx
     do while(x.lt.xb_max)
        y = yb_min + dx
        do while(y.lt.yb_max - 0.5*dx)
           ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0;theta(ipart)=0.0d0
           if(y.le.0.5d0*yl + 0.25d0*yl*sin(2.0*pi*x))then
              theta(ipart) = 1.0d0
           else
              theta(ipart)=-1.0d0
           end if
           y = y + dx
        enddo
        x = x + dx
     enddo
    
     npfb = ipart
     np=ipart


!! ------------------------------------------------------------------------------------------------
  case(5)

     !! POISEUILLE FLOW with circular obstacle...
     !! From Grilli et al (2013) PRL 

     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! |                                        |
     ! |          ----> pressure gradient       |
     ! |                                        |
     ! |                ooo                     |
     ! |               ooooo                    |
     ! |               ooooo                    |
     ! |                ooo                     |
     ! |                                        |
     ! |                periodic bounds in x    |
     ! |                                        |  
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo

     !! Domain size:
     xl=3.0d0;yl=2.0d0   !! Should be of order unity, as simulation is dimensionless
     h0=0.5d0

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 0 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=0
     dt_out = 0.01;tmax = 100.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/00.0d0,0.0d0/)    
     ve_nonlinearity=100d0;alpha_evss=0.0d0
     beta = 0.1d0;Ra=10000.0d0;Pr=1.0d0;El=20.0d0

     !! Resolution:
     dx=h0/20.0;h=1.3d0*dx      !! N.B. you should only need to change dx

     ! Build domain
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2),b_thermal(nb_patches))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 1, 2/)  
     b_periodic_parent(:) = (/ 3, 4, 1, 2/)
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)
      b_thermal(:) = (/ 0.0d0,0.0d0, 0.0d0,0.0d0/)
     nb_circles = 1
     allocate(c_centre(nb_circles,2),c_radius(nb_circles),c_omega(nb_circles),c_vel(nb_circles,2),c_thermal(nb_circles))
     c_centre(1,:) = (/ 0.0d0,0.0d0 /);c_radius(1)=0.25*yl;c_omega(1)=0.0d0;c_vel(1,:)=0.0d0
     c_thermal(:) = 1.0d0     

     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     call make_boundary_particles
     ipart = nbS     

      ! -- fluid particles-- 
     x = xb_min + 0.5*dx
     do while(x.lt.xb_max)
        y = yb_min + dx
        do while(y.lt.yb_max - 0.5*dx)
           if((x-c_centre(1,1))**2+(y-c_centre(1,2))**2.gt.(c_radius(1)+0.25*dx)**2) then  !! Don't fill the circle...
              ipart=ipart+1         
              xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0;theta(ipart)=0.0d0
           end if
           y = y + dx
        enddo
        x = x + dx
     enddo   
     npfb = ipart
     np=ipart
!! ------------------------------------------------------------------------------------------------
  case(6)

     !! Free-surface Rayleigh-Benard flow
     
     !! Domain size:
     xl=1.0d0;yl=xl   !! Should be of order unity, as simulation is dimensionless
     h0=0.5d0         !! Fluid depth

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 1 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=1
     dt_out = 0.1;tmax = 1000.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/0.0d0,-1.0d0/)    
     ve_nonlinearity=0.1d0;alpha_evss=0.0d0
     beta = 0.5d0;Ra=2.0d4;Pr=1.0d0;El=0.01d0

     !! Resolution:
     dx=yl/150.0;h=1.3d0*dx      !! N.B. you should only need to change dx

     !   JRCK boundary conditions
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches),b_vel(nb_patches),b_thermal(nb_patches))
     b_type(:) = (/ 1, 2, 0, 2/)  
     b_periodic_parent(:) = (/ 3, 4, 1, 2/)
     b_vel(:) = (/ 0.0d0,0.0d0,0.0d0,0.0d0/)
     b_thermal(:) = (/ 1.0d0,0.0d0,-1.0d0,0.0d0/)     
     b_node(1,:) = (/ 0d0, 0.0d0 /)
     b_node(2,:) = (/ xl, 0.0d0 /)
     b_node(3,:) = (/ xl, yl /)
     b_node(4,:) = (/ 0.0d0, yl /)

     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     !! outsource making wall particles... JRCK
     call make_boundary_particles
     ipart = nbS     

     ! -- fluid particles-- 
     x = xb_min + 0.5*dx
     do while(x.lt.xb_max)
        y = yb_min + 0.5*dx
        do while(y.lt.h0)
           ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;p(ipart) =0.0d0;up(ipart)=0.0d0;vp(ipart)=0.0d0;theta(ipart)=0.0d0

           if(y.le.0.2*h0 + 0.1*h0*sin(2.0*pi*x/xl)) then
              theta(ipart)=1.0d0
           else
              theta(ipart)=0.0d0
           end if
           y = y + dx
        enddo
        x = x + dx
     enddo
    
     npfb = ipart
     np=ipart
!! ------------------------------------------------------------------------------------------------     
  case(7)

     ! Planar Couette flow
     !  ooooooooooooooooooooooooooooo     
     !          --------->
     ! 
     ! 
     !      <---------
     !  oooooooooooooooooooooooooooooo


     !! Domain size:
     xl=1.0d0;yl=1.0d0   !! Should be of order unity, as simulation is dimensionless
     h0=0.5d0

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 0 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=0
     dt_out = 0.01;tmax = 100.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/0.0d0,0.0d0/)    
     ve_nonlinearity=100d0;alpha_evss=0.0d0
     beta = 0.1d0;Ra=1.0d0;Pr=1.0d0;El=1.0d0

     !! Resolution:
     dx=yl/30.0;h=1.3d0*dx      !! N.B. you should only need to change dx

     !! Build domain
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2),b_vel(nb_patches),b_thermal(nb_patches))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))     
     b_type(:) = (/ 1,2,1,2/)
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_vel(:) = (/ -1.0d0,0.0d0,-1.0d0,0.0d0/)
     b_thermal(:) = (/ 1.0d0,0.0d0,-1.0d0,0.0d0/)
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl,-0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl,0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)

     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     !! outsource making wall particles... JRCK
     call make_boundary_particles
     ipart = nbS

      ! -- fluid particles-- 
     x = xb_min + 0.5*dx
     do while(x.lt.xb_max)
        y = yb_max - dx
        do while(y.gt.yb_min  + 0.5*dx)
              ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0;theta(ipart)=0.0d0
           y = y - dx
        enddo
        x = x + dx
     enddo
    
     npfb = ipart
     np=ipart
!! ------------------------------------------------------------------------------------------------
  case(8)

     !! FREE SURFACE PERIODIC

     ! |----------------------------------------|
     ! |                                        |
     ! |          ----> pressure gradient       |
     ! |                                        |
     ! |                                        |
     ! |                periodic bounds in x    |
     ! |                                        |  
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo

     !! Domain size:
     xl=2.0d0;yl=2.0d0   !! Should be of order unity, as simulation is dimensionless

     !! Simulation controls:
     i_RESTART=.false. !! Option to restart from TEMP_REC file
     shift_coeff = 1 ! choose shift coeff, 0=.5h*h,1=A*dt*h*|u|
     dt_coef_cfl = 0.4d0;dt_coef_visc = 0.2d0
     i_open_domain=1
     dt_out = 0.1;tmax = 100.0d0    !! Output interval and final time
     
     !! Physical parameters
     body_force = (/0.0d0,0.0d0/)    
     ve_nonlinearity=400.0d0;alpha_evss=0.0d0
     beta = 0.1d0;Ra=1000.0d0;Pr=10.0d0;El=1.0d0

     !! Resolution:
     dx=yl/75.0;h=1.3d0*dx      !! N.B. you should only need to change dx
    
     !! Build domain
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2),b_vel(nb_patches),b_thermal(nb_patches))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 0, 2/)
     b_vel(:) = (/ -0.0d0,0.0d0,0.0d0,0.0d0/)
     b_thermal(:) = (/ 1.0d0,0.0d0,-1.0d0,0.0d0/)
     b_periodic_parent(:) = (/ 3, 4, 1, 2/)
     b_node(1,:) = (/ 0.0d0, 0.0d0 /)
     b_node(2,:) = (/ xl, 0.0d0 /)
     b_node(3,:) = (/ xl, yl /)
     b_node(4,:) = (/ 0.0d0, yl /)
     nb_circles = 0
   
     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     call make_boundary_particles
     ipart = nbS

      ! -- fluid particles-- 
     x = xb_min + 0.5*dx
     do while(x.lt.xb_max)
        y = yb_min + dx
        do while(y.lt.0.5d0*yl)
           ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0;theta(ipart)=0.0d0
           if(y.le.0.125d0*yl + 0.0625d0*yl*sin(2.0*pi*x))then
              theta(ipart) = 1.0d0
           else
              theta(ipart)=-1.0d0
           end if
           y = y + dx
        enddo
        x = x + dx
     enddo
    
     npfb = ipart
     np=ipart


  end select
!! ===============================================================================================  
!! ------------------------------------------------------------------------------------------------
  !! Write info to screen...
  write(6,*) "beta = ",beta,"El=",El
  write(6,*) "Ra=",Ra,"Pr=",Pr
  write(6,*) "Reynolds number =",sqrt(Ra/Pr)
  write(6,*) "Weissenberg number = ",El*sqrt(Ra/Pr)

!! ------------------------------------------------------------------------------------------------

  !! write data to files


  !! Particle data  
  open(13,file='./IPART')
  do i=1,np
     write(13,*) xp(i), yp(i), up(i), vp(i), theta(i), dx*dx
     
  end do

  close(13)
  deallocate(xp, yp,up,vp,p,vol,theta)
  write(*,*) 'np = ',np,"nbS = ",nbS
  print*,'dx',dx,'h',h


  !!  output run parameters
  open(20,file='./INDAT')
  write(20,*)i_RESTART,"        :i_RESTART"
  write(20,'(I5,A)') i_open_domain,"        :open domain (0-closed,1-open)"
  write(20,'(F18.10,A)') tmax,"        :tmax - max time of simulation"
  write(20,'(I5,A)')shift_coeff,"        :shift_coeff - how much to shift"
  write(20,'(I7,A)') np,"        :np"
  write(20,'(I7,A)') npfb,"        :npfb"
  write(20,'(I7,A)') nbS,"        :nbS"
  write(20,'(F18.10,A)') dt_out," :dt_out - time interval between outputs"
  write(20,'(F18.10,A)') h," :h"
  write(20,'(F18.10,A)') dx," :dx"
  write(20,'(F18.10,A)') dt," :dt"
  write(20,'(F18.10,F18.6,A)') body_force(1),body_force(2)," :body_force"
  write(20,'(F18.10,A)') Ra," Rayleigh Number"
  write(20,'(F18.10,A)') Pr," Prandtl Number"
  write(20,'(F18.10,A)') El," Prandtl Number"
  write(20,'(F18.10,A)') beta," Viscosity ratio"   
  write(20,'(F18.10,F18.10,A)') dt_coef_cfl,dt_coef_visc," :dt_coefs cfl,visc"
  write(20,'(F18.10,A)') ve_nonlinearity,"  :ve_nonlinearity parameter"
  write(20,'(F18.10,A)') alpha_evss,"  :EVSS coefficient"     
  close(20)


  !!  output boundary data
  open(21,file='./IBOUND')
  write(21,*) nb_patches,"                                           :nb_patches"              
  do i=1,nb_patches
     write(21,*) b_node(i,:),"   :b_node",i
     write(21,*) b_edge(i,:),"   :b_edge",i
     write(21,*) b_type(i),"                                           :b_type",i
     if(allocated(b_vel))then
        write(21,*) b_vel(i),"   :b_vel",i
     else
        write(21,*) 0.0d0
     end if
     if(allocated(b_thermal))then
        write(21,*) b_thermal(i),"   :b_thermal",i
     else
        write(21,*) 0.0d0
     end if
     if(b_type(i).eq.2)then
        write(21,*) b_periodic_parent(i), &
             "                                           :b_periodic_parent"
     end if
     if(b_type(i).eq.3)then
        write(21,*) U_bulk,"                     :U_inflow"
     end if
  end do
  write(21,*) nb_circles
  if(nb_circles.ne.0)then
     do i=1,nb_circles
        write(21,*) c_centre(i,:),"    :c_centre",i
        write(21,*) c_radius(i),"      :c_radius",i
        write(21,*) c_omega(i),"       :c_omega",i
        write(21,*) c_vel(i,:),"       :c_vel",i
        write(21,*) c_thermal(i),"   :c_thermal",i
     end do
  end if
  close(21)
  deallocate(b_node,b_edge)
  deallocate(b_type,b_periodic_parent)

  !! Everything done, finish!!

  write(*,*) 'END of DATCLASS'
  stop
end program datgen
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_edge_vectors
     !! Take boundary nodes, and construct boundary edge vectors
     use kind_parameters
     use common_parameter
     use global_variables
     implicit none
     integer(ikind) ib,ibp1

     do ib = 1,nb_patches ! loop over all boundary patches
        ibp1 = mod(ib,nb_patches) + 1   
        b_edge(ib,:) = b_node(ibp1,:) - b_node(ib,:)  ! calculate b_edge
     end do

     return 
   end subroutine make_boundary_edge_vectors
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_particles
     !! Run along boundary edge vectors and place nodes (spaced dx) along all walls...
     use kind_parameters
     use common_parameter
     use global_variables 

     implicit none
     integer(ikind) :: ipart,ib,ibm1
     real(rkind) :: x,y,m_be,tmp,tmp2
     integer(ikind) :: nround,iround

     if(allocated(b_vel).eqv..false.)then
        allocate(b_vel(nb_patches));b_vel(:)=0.0d0
     end if
     
     ipart=0
     allocate(xp(npar), yp(npar),up(npar), vp(npar),p(npar),theta(npar), vol(npar))

     ! circle particles
     do ib=1,nb_circles
        tmp = 2.0*abs(c_radius(ib))*pi;nround=ceiling(tmp/dx)  !! circumference and number round circle
        do iround =1,nround
           tmp2 = iround*2.0*pi/nround
           x = c_centre(ib,1) + abs(c_radius(ib))*cos(tmp2)
           y = c_centre(ib,2) + abs(c_radius(ib))*sin(tmp2)
           ipart = ipart + 1;xp(ipart)=x;yp(ipart)=y
           up(ipart) = c_vel(ib,1)+ c_omega(ib)*abs(c_radius(ib))*sin(tmp2);
           vp(ipart)=c_vel(ib,2)-c_omega(ib)*abs(c_radius(ib))*cos(tmp2) !! Velocity...
           p(ipart) = 0.0d0
           theta(ipart)= c_thermal(ib)
        end do
     end do
     !! Wall particles
     do ib=1,nb_patches  ! loop over all boundary patches
        ibm1 = mod(ib+nb_patches-2,nb_patches)+1
        if(abs(b_type(ib)).eq.1)then  ! if it is a wall boundary patch
           m_be = dsqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
           x = b_node(ib,1);y = b_node(ib,2)
           tmp = 0.0
           do while(tmp.lt.1.0-1.0d-10)   ! move along the patch in increments of dx
              ipart = ipart + 1
              xp(ipart) = x;yp(ipart) = y
              up(ipart) = b_vel(ib)*b_edge(ib,1)/m_be!0.0d0
              vp(ipart) = b_vel(ib)*b_edge(ib,2)/m_be!0.0d0
              if(abs(b_type(ibm1)).eq.1.and.b_vel(ibm1).ne.0.0.and.tmp.eq.0.0) then !! If prev patch is MOVING wall and this is particle 1 on current patch...
                 up(ipart) = b_vel(ibm1)*b_edge(ibm1,1)/dsqrt(dot_product(b_edge(ibm1,:),b_edge(ibm1,:)))
                 vp(ipart) = b_vel(ibm1)*b_edge(ibm1,2)/dsqrt(dot_product(b_edge(ibm1,:),b_edge(ibm1,:)))
              end if
              p(ipart) = 0.0d0
              theta(ipart)= b_thermal(ib)
              tmp = tmp + dx/m_be  
              x = b_node(ib,1) + tmp*b_edge(ib,1)
              y = b_node(ib,2) + tmp*b_edge(ib,2)
           end do
        else if(abs(b_type(ibm1)).eq.1)then  ! if the previous boundary patch was a wall (but ib is not a wall...)
           ipart = ipart + 1
           xp(ipart) = b_node(ib,1)  ! place a single particle at the node (to complete the previous wall)
           yp(ipart) = b_node(ib,2)
           m_be = dsqrt(dot_product(b_edge(ibm1,:),b_edge(ibm1,:)))
           up(ipart) = b_vel(ibm1)*b_edge(ibm1,1)/m_be!0.0d0
           vp(ipart) = b_vel(ibm1)*b_edge(ibm1,2)/m_be!0.0d0
           p(ipart) = 0.0d0
           theta(ipart)= b_thermal(ibm1)    !Write 
        end if
     end do
     nbS=ipart    

     write(*,*) 'no. of solid boundary particles: ',nbS
     return
   end subroutine make_boundary_particles
!! ------------------------------------------------------------------------------------------------
   subroutine is_inside_bounds(x,y,inflag)
     !! Only for simple-ish domains: check whether a point is inside or outside the boundaries, return a flag
     use kind_parameters
     use common_parameter
     use global_variables 
     implicit none
     real(rkind), intent(in) :: x,y
     logical, intent(out) :: inflag
     integer(ikind) :: ib,ibp1
     real(rkind) :: p_dist,tmp1,pap
     real(rkind),dimension(dims) :: nrm
     
     inflag = .true.
     do ib=1,nb_patches   ! loop over each patch
        ibp1 = mod(ib,nb_patches)+1
        !! Step 1: calculate the perpendicular distance...
        p_dist = b_edge(ib,2)*x - b_edge(ib,1)*y + &
            b_node(ibp1,1)*b_node(ib,2) - b_node(ibp1,2)*b_node(ib,1)
        p_dist = -1.0d0*p_dist/dsqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))

        !! Step 2: calculate the boundary normal...
        tmp1 = dsqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
        nrm(1) = -1.0*b_edge(ib,2)/tmp1
        nrm(2) = b_edge(ib,1)/tmp1

        !! Step 3: Calculate whether it is between the ends...
        if(abs(b_edge(ib,1)).gt.abs(b_edge(ib,2)))then
           pap = (x - b_node(ib,1) - p_dist*nrm(1))/b_edge(ib,1)
        else
           pap = (y - b_node(ib,2) - p_dist*nrm(2))/b_edge(ib,2)
        end if
 
        !! Step 4: is it between the ends and on the wrong side
        if(pap.ge.0.0d0.and.pap.le.1.0d0.and.p_dist.le.0.25d0*dx)then
           inflag = .false.
        end if

     end do

   end subroutine is_inside_bounds
!! ------------------------------------------------------------------------------------------------
