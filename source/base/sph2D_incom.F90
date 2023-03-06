program isph
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!         Viscoelastic Incompressible Smoothed Particle Hydrodynamics                         !!
  !!         -----------------------------------------------------------                         !!
  !!                                                                                             !!
  !!         AUTHOR: Jack King                                                                   !!
  !!                                                                                             !!
  !!         Department of Mechanical, Aerospace and Civil Engineering                           !!
  !!         University of Manchester                                                            !!
  !!                                                                                             !!
  !!                                                                                             !!
  !!         RELEVANT REFERENCES:                                                                !!  
  !!         1) Lind et al. (2012) J. Comput. Phys. 231:1499-1523                                !!
  !!         2) King & Lind (2020) ArXiv: 2009.12245 (https://arxiv.org/abs/2009.12245)          !!
  !!                                                                                             !!
  !!                                                                                             !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use kind_parameters
  use common_parameter
  use common_2d
  use input
  use cell_griddin
  use omp_lib
  use divfree

  implicit none
  real(rkind) :: x_mint,y_mint,x_maxt,y_maxt
  real(rkind) :: ts_start,ts_end,t_run,t_per_dt,t_last_10
  integer(ikind) :: n_threads,narg
  character(len=200) :: argch

  !! Set the number of threads:
  narg = command_argument_count()
  if(narg.ne.0)then
     call get_command_argument(1,argch)
     argch = adjustl(argch)
     read(argch,938) n_threads
     938 format(i4)     
     call omp_set_num_threads(n_threads)
     write(6,*) "User specified ",n_threads,"threads"
  else
     !$omp parallel
     n_threads = omp_get_num_threads()
     !$omp end parallel
     write(6,*) "Automatically allocated",n_threads,"threads"
  end if

  ! ------------------------------------------
  ! ------------------------------------------
  ! Read in data and initialize the variables
  call getdata()
  !! open some files for outputting
  if(i_RESTART) then
     open(unit=92,file='./data_directory/TIME',status='unknown',access='append')
     open(unit=93,file='./data_directory/LS',status='unknown',access='append')
     open(unit=94,file='./data_directory/TIME_OUT',status='unknown',access='append')
  else
     open(unit=92,file='./data_directory/TIME')
     open(unit=93,file='./data_directory/LS')
     open(unit=94,file='./data_directory/TIME_OUT')
  end if

  ! Adjust the domain size to give room for mirror particles
  x_mint=xb_min - sup_size*1.5d0
  x_maxt=xb_max + sup_size*1.5d0
  y_mint=yb_min - sup_size*1.5d0
  y_maxt=yb_max + sup_size*1.5d0

  !! Split the domain into a grid of cells, which are used for efficient particle-neighbour finding later
  call cell_generator(x_mint,x_maxt,y_mint,y_maxt,uno_sup_size,ncx,ncy,nct)
  ncell_m=nct
  nsheet_m = nct
  allocate(xc_m(nct),yc_m(nct))   !! find the positions of the cells 
  call cell_position(x_mint,x_maxt,y_mint,y_maxt,sup_size,     &
       xc_m,yc_m,ncx,ncy,nct)   !! cell positions with mirrors 

  !! Prepare for prifiling and output counting
  t_run = 0.0d0
  t_last_10 = 0.0d0
  next_dump = time;next_dump2 = time
  if(i_RESTART.eqv..false.)then
     n_dump = 0
     n_dump2= 0
  end if
  
  ! Output data to file
  call output  

  !! MAIN LOOP ====================================================================================
  do while (time.le.tmax)
     itime = itime+1

     ! profiling
     ts_start=omp_get_wtime()
     
     !! A single time step
     call div_free(x_mint,y_mint,x_maxt,y_maxt)
     time = time + dt

     ! profiling
     ts_end=omp_get_wtime()
     t_run = t_run + ts_end - ts_start
     t_per_dt = t_run/dble(itime)
     t_last_10 = t_last_10 + ts_end - ts_start

     ! Output data to screen
     if(mod(itime,10).eq.0)then
        write(6,*)"itime,time,dt=", itime,time,dt
        write(6,*) "np,npfb,nmirror",np,npfb,nmirror
        write(6,*) "Linear solver:",LS_iters,"iterations,error:",LS_residual
        write(6,*) "There are ",n_inbin,"particles in the bin"
        !$OMP PARALLEL
        n_threads = omp_get_num_threads()
        !$OMP END PARALLEL
        write(6,*) "Number of threads=",n_threads,"Run time=",t_run
!        write(6,*) "Run time=",t_run,"run-time/step=",t_per_dt
        write(6,*) "run-time/dt=",t_per_dt,"Moving avg=",0.1d0*t_last_10
        t_last_10 = 0.0d0
        write(6,'(/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,/,A)') "  "
     end if

     ! Output data to file
     call output
     
  end do
  ! END OF MAIN LOOP ==============================================================================

  ! Output at final step
  call output

  close(92)  ! data_directory/TIME for time,dt,umax,vmax,pmax
  close(93)  ! data_directory/LS for iters and residual every dt
  close(94)  ! data_directory/TIME_OUT for time,npfb every dump

  deallocate(r,u)
  deallocate(p)
  deallocate(xc_m,yc_m)
  deallocate(inbin)

end program isph
