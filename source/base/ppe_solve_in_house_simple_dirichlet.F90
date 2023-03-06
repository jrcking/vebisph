subroutine ppe_solve_in_house_simple_dirichlet

  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  use omp_lib
  use linear_solver
  implicit none


  real(rkind) :: temp
  integer(ikind) :: itermax
  real(rkind) :: tol

  integer ind,i,j,k,num,ij,n_t,nnz,nnzold
  double precision unorobar
  real(rkind),dimension(dims) :: rij,uij
  real(rkind) :: rr2_tmp, rad, qq,Factor_B_tmp,Adiag
  integer(ikind), dimension(:),allocatable :: ij_num,ia

  allocate(Factor_B(npfb))
  ! Initialise arrays
  Factor_B=0.0d0 !RHS of PPE (div u)

  !! PART 1: Calculate the RHS of the PPE (div.u*)
  temp = 0.0d0
  !$OMP PARALLEL DO PRIVATE(k,j,uij,Factor_B_tmp) reduction(+:temp)
  do i=1,npfb
     Factor_B_tmp = 0.0d0
     do k=1,ij_count(i)
        j=ij_link(i,k)
    
        uij(:) = u(i,:) - u(j,:)

        !R.H.S. Poisson equation cf Eqn 15 in Lind et al. 2011

        Factor_B_tmp=Factor_B_tmp - dot_product(uij(:),matmul(kgcm(i,:,:),ij_w_G(i,k,:)))
     enddo
     Factor_B(i) = Factor_B_tmp/dt
     temp = temp + Factor_B(i)
  enddo
  !$OMP END PARALLEL DO

  !! Remove DC from Factor_B if closed domain
  if(i_open_domain.eq.0)then
     temp = temp/dble(npfb)
     !$OMP PARALLEL DO
     do i=1,npfb
        Factor_B(i)=Factor_B(i)-temp
     end do
     !$OMP END PARALLEL DO
  end if

  !! PART 2: Construct the LHS of the Pressure Poisson Eqn using Morris operator

  !! index of first element in each row... 
  !! This is the only loop which we have to do in serial
  !! Could be parallelized as in Sphysics-woof, but it's actually not very expensive this way
  allocate(ia(npfb+1))
  ia(1)=npfb+2
  num=npfb+1
  do i=1,npfb
     ia(i+1) = ia(i) + ij_count(i)
  end do

  !! Size of ppe and allocation...
  nnz = ia(npfb+1)
  allocate(A(nnz),ija(nnz));ija=0;A=0.0d0
  allocate(ij_num(npfb));ij_num=0

  !$omp parallel do private(num,k,j,ij_num,unorobar,temp,Adiag,ij,n_t)
  do i=1,npfb   ! for all particles except compressible particles
     ija(i) = ia(i)
     !! First embedded loop builds ija so we can track where things go...
     num = ija(i) - 1
     do k=1,ij_count(i)   ! for all neighbours
        j=ij_link(i,k)     ! now considering particles i and j
        if(j.le.npfb)then
           num=num+1          
           ij_num(j)=num
           ija(num)=j
        endif
        if(j.gt.npfb)then
           num=num+1
           ij = irelation(j)
           ij_num(ij)=num
           ija(num) = ij
        end if
     enddo

     !! Second embedded loop is used to populate A
     Adiag=0.0_rkind
     do k=1,ij_count(i)
        j=ij_link(i,k)
        !Morris operator for Laplacian (cf. Eqn 8, Lind et al. 2011)
        temp = -ij_w_L(i,k)
        if(i.ne.j)then
           Adiag=Adiag - temp   ! Contribution to diagonal of A
           if(j.gt.npfb)then    ! If j is mirror, then
              ij=irelation(j)    ! look for parent of j, ij
              if(i.ne.ij)then
                 n_t = ij_num(ij)
                 A(n_t)=A(n_t)+temp   ! Contribution to A(i,ij)
              else
                 Adiag=Adiag+temp       ! Contribution to A(i,ij) if ij=i 
              end if
              ! Building pressure boundary condition
              Factor_B(i)=Factor_B(i)-dP_mp(j)*temp
           elseif(j.le.npfb)then   ! if j is not a mirror
              n_t = ij_num(j)   
              A(n_t)=A(n_t)+temp    ! Contribution to A(i,j)
           end if
        end if
     end do
     if(n_surf(i).ne.0)then
        Adiag=1.0d30            !! Simple dirichlet...
     end if
     A(i)=A(i)+Adiag
  end do
  !$omp end parallel do
  ija(npfb+1)=ia(npfb+1)  !! Final element of ija
  
  
!   if(.true..and.twophase)then   !! PPE SMOOTHING for low viscosity flows with twophase code..
  if(.true..and.i_open_domain.eq.1)then
     !$omp parallel do private(j)
     do i=1,npfb
        if(n_surf(i).eq.1) cycle ! already dealt with by Dirichlet BC
        if(ppe_smooth(i).eq.1.0_rkind) cycle  !! not near the FS
        Factor_B(i) = Factor_B(i)*ppe_smooth(i)
        do j=ija(i),ija(i+1)-1
           A(j) = A(j)*ppe_smooth(i)
        end do
     end do
     !$omp end parallel do
     deallocate(ppe_smooth)
  end if

  deallocate(ij_num,ia)
  
  !! Solve the PPE....
  tol = 1.0d-8!1.d-10
  itermax = 1000
  p(1:npfb) = 0.0d0

   
  if(i_open_domain.eq.0)then
     call SPARSE_PBICGSTAB_nonull(A,ija,P,Factor_B,npfb,nnz,itermax,tol,LS_iters,LS_residual)
  else
     call SPARSE_PBICGSTAB(A,ija,P,Factor_B,npfb,nnz,itermax,tol,LS_iters,LS_residual)
  end if

  deallocate(Factor_B,A,ija)    ! can safely deallocate A,B,ija here
  
!  write(6,*) LS_iters,LS_residual

end subroutine ppe_solve_in_house_simple_dirichlet
