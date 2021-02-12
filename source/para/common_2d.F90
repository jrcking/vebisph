module common_2d

  use iso_c_binding
  use common_parameter
  use kind_parameters      

  implicit none

  !! Lagrangian or Eulerian?
  logical :: lagrangian

  !! Dimensionless groups
  real(rkind) :: beta,Pr,Ra,El,eta,lambda
  
  !! Forcing, etc
  real(rkind),dimension(2) :: body_force
  real(rkind) :: ve_nonlinearity,alpha_evss
  
  !! Primary fluid properties
  real(rkind),dimension(:,:),allocatable :: r, u   !! Position and velocity
  real(rkind),dimension(:,:,:),allocatable :: tau_p !! polymeric stress
  real(rkind),dimension(:),allocatable :: theta     !! Temperature (relative, dimensionless)
  
  !! Secondary fluid properties
  real(rkind),dimension(:),allocatable :: p,trace_conf,vort,a_out
  
  !! Tertiary fields: related to temporal and spatial discretisation
  real(rkind),dimension(:),allocatable :: div_r,conc  !! Distribution divergence and density
  real(rkind),dimension(:,:,:),allocatable :: kgcm   !! Kernel correction matrices
  integer(ikind),dimension(:),allocatable :: n_surf   !! Free-surface flag
  real(rkind), dimension(:,:), allocatable :: r0, u0         !! Position and velocity at start of each time step
  logical( c_bool ), dimension(:), allocatable, target :: inbin  !! Bin flag
  real(rkind),allocatable,dimension(:) :: ppe_smooth
  integer(ikind), dimension(:), allocatable, target :: irelation,vrelation
  real(rkind),dimension(:),allocatable, target :: dP_mp
  real(rkind),dimension(:,:),allocatable, target :: surf_norm
  real(rkind), dimension(:,:), allocatable,target :: ushift,dr_fs !! For Fickian shifting

  !! Spatial and temporal discretisation
  real(rkind) :: dx, dv, av_conc,dt,dt_coef_cfl,dt_coef_visc ! space and time step related
  real(rkind) :: umax,time,tmax, dt_out,next_dump,next_dump2
  integer(ikind) :: itime,n_dump,n_dump2
  integer(ikind) :: np,npfb,nbS,nmirror,n_count    !! Numbers of particles etc.
  integer(ikind) :: nmirror_esti,npfb_old,maxneighbours,maxneighbours_old,n_inbin
  real(rkind) :: xmin,xmax,ymin,ymax        !! Domain size
  real(rkind) :: h,h2,uno_sup_size,sup_size,sup_size_2,eta2 !! Smoothing length and support size
  real(rkind) :: xb_min,xb_max,yb_min,yb_max,zb_min,zb_max  !! domain size (again)
  real(rkind) :: ad_7,ad_7h                !! Kernel pre-multiplication factors

  !! Switches and flags
  logical :: external_forcing,akinci,output_everytime,fick_shifting
  logical :: use_labfm,output_mirrors,i_RESTART  
  integer(ikind) ::  i_open_domain,i_shift_coeff

  !! Boundary condition framework
  real(rkind),dimension(:,:),allocatable, target :: b_node,b_edge,c_centre,c_vel,c_vel0,c_centre0
  integer(ikind),dimension(:),allocatable,target :: b_type,b_periodic_parent
  real(rkind),dimension(:),allocatable,target :: b_corner_angle,c_radius,b_vel,c_omega
  integer(ikind),dimension(:),allocatable,target :: p_inflow
  real(rkind) :: U_inflow,U_inflow0
  integer(ikind) :: nb_patches,np_inflow,i_pset_part,ip_dirichlet,nb_circles
  integer(ikind),dimension(:),allocatable,target :: p_free
  integer(ikind) :: n_free

  !! linear solver and PPE stuff
  real(rkind), dimension(:), allocatable, target :: A,Factor_B
  integer(ikind), dimension(:), allocatable, target :: ija
  integer(ikind) :: LS_iters
  real(rkind) :: LS_residual

  !! numbers of cells, including mirrors
  integer(ikind) :: ncx,ncy,nct,nsheet_m
  real(rkind), dimension(:),allocatable :: xc_m,yc_m
  integer(ikind) ::  ncell,ncell_m  

  !! neighbour variables
  integer(ikind),dimension(:),allocatable,target :: ij_count
  integer(ikind),dimension(:,:),allocatable,target :: ij_link,ij_ind,ij_cell,ij_slip
  real(rkind),dimension(:,:),allocatable :: ij_w_L
  real(rkind),dimension(:,:,:),allocatable :: ij_w_G
end module common_2d
