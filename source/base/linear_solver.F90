module linear_solver
 use kind_parameters
 
contains

subroutine norm(r,n,norm_r)
  implicit none

  integer i,n
  real(rkind) r(n)
  real(rkind) norm_r

  norm_r=0.0d0
  !$OMP PARALLEL DO REDUCTION(+:norm_r)
  do i=1,n
     norm_r=norm_r+r(i)*r(i)
  enddo
  !$OMP END PARALLEL DO

  return
end subroutine norm

!-------------------------------------------

subroutine vtime(R,X,RX,n)
  implicit none 

  integer n,i
  real(rkind) R(n),X(n)
  real(rkind) temp,RX

  temp = 0.0d0
  !$OMP PARALLEL DO REDUCTION(+:temp)
  do i=1,n
     temp=temp+R(i)*X(i)
  enddo
  !$OMP END PARALLEL DO
  RX=temp  

  return
end subroutine vtime

!-------------------------------------------
!-------------------------------------------
!-------------------------------------------

subroutine inverseM(M,X,inv_X,n)
  implicit none

  integer n,i
  real(rkind) M(n),X(n),inv_X(n)

  !$OMP PARALLEL DO
  do i=1,n
     inv_X(i)=X(i)/(M(i)+1.0d-15)
  enddo
  !$OMP END PARALLEL DO

end subroutine inverseM

!-----------------------------------------------------------------

subroutine sparse_make_preconditioner  &
     &              (M,A,n,nmax) ! Jacobi preconditioner is made here
  implicit none
  integer n,nmax,i
  real(rkind) M(n),A(nmax)
 
  !$OMP PARALLEL DO
  do i=1,n
     M(i)=A(i)
  enddo
  !$OMP END PARALLEL DO

end subroutine sparse_make_preconditioner

!-------------------------------------------
!-------------------------------------------
!               SPARSE_PBiCGSTAB
!-------------------------------------------

subroutine SPARSE_PBICGSTAB(A,ija,x,B,n,nmax,itermax,tol,iter,residial)
  implicit none

  integer ::  n,nmax,i,itermax, iter, ija(nmax)
  real(rkind) :: tol,residial
  real(rkind) ::  temp,rho,rho0,alpha,omega,beta,  &
       &                  norm2_t,norm2_r,norm2_r0, &
       &                  norm2_s,norm2_s0 
  real(rkind) :: r(n),r_bar(n),  &
       &                  B(n), A(nmax), &
       !     &                  B(n),A(nmax),ija(nmax), &
       &                  x(n),v(n),AX(n), &
       &                  p(n),p_inv(n),s(n), &
       &                  s_inv(n),M(n),t(n)

  temp=0.0d0
  rho0=1.0d0
  alpha=1.0d0
  omega=1.0d0
  rho=0.0d0

  do i=1,n
     r(i)=0.0d0
     r_bar(i)=0.0d0
     v(i)=0.0d0
     p(i)=0.0d0
     p_inv(i)=0.0d0
     s(i)=0.0d0
     s_inv(i)=0.0d0
     AX(i)=0.0d0
     M(i)=0.0d0
     t(i)=0.0d0
  enddo

  call sparse_make_preconditioner  &
       &      (M,A,n,nmax)

  call sparse_atime(A,ija,x,AX,n,nmax)

  !$OMP PARALLEL DO
  do i=1,n
     r(i)=B(i)-AX(i)    
  enddo
  !$OMP END PARALLEL DO

  r_bar=r               

  iter=1
  norm2_r0=1.0d0
  norm2_s0=1.0d0

  do while(iter .le. itermax)

     call vtime(r_bar,r,rho,n)      !! wiki PBiCBSTAB 1

     if(rho==0.0d0) then
        print*,'BICGSTAB Fails'
        goto 10
     endif

     if (iter==1)then
        do i=1,n
           p(i)=r(i)
        enddo
     else
        beta=rho/(rho0+1.0d-15)*  &
             &        alpha/(omega+1.0d-15)                 !! wiki PBiCBSTAB 2
        !$OMP PARALLEL DO
        do i=1,n
           p(i)=r(i)+beta*(p(i)-omega*v(i))                 !! wiki PBiCBSTAB 3   p is the search direction.. (search vector)
        enddo
        !$OMP END PARALLEL DO
     endif

     call inverseM(M,p,p_inv,n)                             !! wiki PBiCBSTAB 4
     call sparse_atime(A,ija,p_inv,v,n,nmax)               !! wiki PBiCBSTAB 5
     call vtime(r_bar,v,temp,n)                             !! wiki PBiCBSTAB 6a

     alpha=rho/(temp+1.0d-15)                                !! wiki PBiCBSTAB 6b

     !$OMP PARALLEL DO
     do i=1,n
        s(i)=r(i)-alpha*v(i)                               !! wiki PBiCBSTAB 9 or 7?
     enddo
     !$OMP END PARALLEL DO

     call norm(s,n,norm2_s)                      !! wiki PBiCBSTAB 8

     if (iter==1) then
        if (norm2_s .le. 10-6) then
           norm2_s0=1.0d0
        else
           norm2_s0=norm2_s
        end if
     endif

     if (norm2_s.le.tol*tol*norm2_s0) then
        do i=1,n 
           x(i)=x(i)+alpha*p_inv(i)
        enddo
        residial=sqrt(norm2_s/norm2_s0)
        !          print*,'iter,residual',iter,residial
        goto 10
     endif

     call inverseM(M,s,s_inv,n)             !! wiki PBiCBSTAB 10
     call sparse_atime(A,ija,s_inv,t,n,nmax)       !! wiki PBiCBSTAB 11
     call vtime(t,s,temp,n)                      !! wiki PBiCBSTAB 12a
     call norm(t,n,norm2_t)                      !! wiki PBiCBSTAB 12b

     omega=temp/(norm2_t+1.0d-15)                    !! wiki PBiCBSTAB 12c
     !$OMP PARALLEL DO
     do i=1,n
        x(i)=x(i)+alpha*p_inv(i)+omega*s_inv(i)       !! wiki PBiCBSTAB 13
     enddo
     !$OMP END PARALLEL DO

     !$OMP PARALLEL DO
     do i=1,n
        r(i)=s(i)-omega*t(i)                            !! wiki PBiCBSTAB 15
     enddo
     !$OMP END PARALLEL DO
     call norm(r,n,norm2_r)

     if (iter==1) then
        if (norm2_r.le.10-6) then
           norm2_r0=1.0d0
        else
           norm2_r0=norm2_r
        end if
     endif

     residial=sqrt(norm2_r/norm2_r0)
     !        print*,'iter,residual',iter,residial

     if (norm2_r.le.tol*tol*norm2_r0) goto 10                !! wiki PBiCBSTAB 14

     iter=iter+1

     rho0=rho

  enddo

10 continue

end subroutine SPARSE_PBICGSTAB

!-------------------------------------------

subroutine remove_dc(a,n)
  implicit none
  integer ::  n
  real(rkind) ::  a(n),tmp
  integer :: i,k

  tmp=0.0d0
  !$OMP PARALLEL DO REDUCTION(+:tmp)
  do i=1,n
     tmp = tmp + a(i)
  end do
  !$OMP END PARALLEL DO
  tmp = tmp/dble(n)
  
  !$OMP PARALLEL DO
  do i=1,n
     a(i)=a(i)-tmp
  end do
  !$OMP END PARALLEL DO
end subroutine remove_dc

!-------------------------------------------

subroutine sparse_atime(A,ija,x,b,n,nmax)
  implicit none
  integer ::  n, nmax, ija(nmax)
  real(rkind) ::  b(n),x(n),A(nmax),bb
  integer :: i,k

  if(ija(1).ne.n+2)then
     print*,'Sparse storage mode is not configured properly!'
  endif
  !        b=0.0d0
  !$OMP PARALLEL DO PRIVATE(k,bb)
  do i=1,n
     bb=A(i)*x(i)
     do k=ija(i),ija(i+1)-1
        bb=bb+A(k)*x(ija(k))
     enddo
     b(i) = bb
  enddo
  !$OMP END PARALLEL DO
end subroutine sparse_atime
! ------------------------------------------------------
subroutine SPARSE_PBICGSTAB_nonull(A,ija,x,B,n,nmax,itermax,tol,iter,residial)
  use omp_lib
  implicit none

  real(rkind),dimension(:),intent(in)    :: A,B
  integer(ikind),dimension(:),intent(in) :: ija
  integer(ikind),intent(in)              :: n,nmax,itermax
  real(rkind),intent(in)                 :: tol
  integer(ikind),intent(out)             :: iter
  real(rkind),intent(out)                :: residial
  real(rkind),dimension(:),intent(inout) :: x
  integer(ikind)                         :: i
  
  real(rkind) ::  temp,rho,rho0,alpha,omega,beta,  &
       &                  norm2_t,norm2_r,norm2_r0, &
       &                  norm2_s,norm2_s0 
  real(rkind),dimension(:),allocatable :: r,r_bar,v,AX,p,p_inv,s,M,s_inv,t
  real(rkind) ::  eps_small
  
  !! Set a small number for mollification of div-by-zero
  eps_small = 1.0d-15!epsilon(tol) 

  allocate(r(n),r_bar(n),v(n),AX(n),p(n),p_inv(n),s(n),s_inv(n),t(n))
  allocate(M(n))

  temp=0.0d0
  rho0=1.0d0
  alpha=1.0d0
  omega=1.0d0
  rho=0.0d0

  call sparse_atime(A,ija,x,AX,n,nmax)

  !$OMP PARALLEL DO
  do i=1,n
     r(i)=B(i)-AX(i)    
     r_bar(i)=r(i)
     v(i)=0.0d0
     p(i)=0.0d0
     p_inv(i)=0.0d0
     s(i)=0.0d0
     s_inv(i)=0.0d0
     AX(i)=0.0d0
     M(i)=1.0d0!A(i)   !!! For very low Re, closed domains (e.g. Poiseuille Re=0.01) preconditioner makes poorly conditioned!
     t(i)=0.0d0
  enddo
  !$OMP END PARALLEL DO

  iter=1
  norm2_r0=1.0d0
  norm2_s0=1.0d0
  
  do while(iter .le. itermax)


     call vtime(r_bar,r,rho,n)      !! wiki PBiCBSTAB 1

     if(rho==0.0d0) then
        print*,'BICGSTAB Fails'
        goto 10
     endif

     if (iter==1)then
        do i=1,n
           p(i)=r(i)
        enddo
     else
        beta=rho/(rho0+1.0d-15)*  &
             &        alpha/(omega+1.0d-15)                 !! wiki PBiCBSTAB 2
        !$OMP PARALLEL DO
        do i=1,n
           p(i)=r(i)+beta*(p(i)-omega*v(i))                 !! wiki PBiCBSTAB 3   p is the search direction.. (search vector)
        enddo
        !$OMP END PARALLEL DO
     endif

     call remove_dc(p,n)                            !! Remove nullspace from search direction
     call inverseM(M,p,p_inv,n)                             !! wiki PBiCBSTAB 4
     call sparse_atime(A,ija,p_inv,v,n,nmax)               !! wiki PBiCBSTAB 5
     call vtime(r_bar,v,temp,n)                             !! wiki PBiCBSTAB 6a

     alpha=rho/(temp+1.0d-15)                                !! wiki PBiCBSTAB 6b

     !$OMP PARALLEL DO
     do i=1,n
        s(i)=r(i)-alpha*v(i)                               !! wiki PBiCBSTAB 9 or 7?
     enddo
     !$OMP END PARALLEL DO

     call norm(s,n,norm2_s)                      !! wiki PBiCBSTAB 8

     if (iter==1) then
        if (norm2_s .le. 1.0d-6) then  !     
           norm2_s0=1.0d0
        else
           norm2_s0=norm2_s
        end if
     endif

     if (norm2_s.le.tol*tol*norm2_s0) then
        do i=1,n 
           x(i)=x(i)+alpha*p_inv(i)
        enddo
        residial=sqrt(norm2_s/norm2_s0)
        !          print*,'iter,residual',iter,residial
        goto 10
     endif

     call remove_dc(s,n)      !! Remove nullspace from search direction
     call inverseM(M,s,s_inv,n)             !! wiki PBiCBSTAB 10
     call sparse_atime(A,ija,s_inv,t,n,nmax)       !! wiki PBiCBSTAB 11
     call vtime(t,s,temp,n)                      !! wiki PBiCBSTAB 12a
     call norm(t,n,norm2_t)                      !! wiki PBiCBSTAB 12b

     omega=temp/(norm2_t+1.0d-15)                    !! wiki PBiCBSTAB 12c
     !$OMP PARALLEL DO
     do i=1,n
        x(i)=x(i)+alpha*p_inv(i)+omega*s_inv(i)       !! wiki PBiCBSTAB 13
     enddo
     !$OMP END PARALLEL DO

     !$OMP PARALLEL DO
     do i=1,n
        r(i)=s(i)-omega*t(i)                            !! wiki PBiCBSTAB 15
     enddo
     !$OMP END PARALLEL DO
     call norm(r,n,norm2_r)

     if (iter==1) then
        if (norm2_r.le.1.0d-6) then !        
           norm2_r0=1.0d0
        else
           norm2_r0=norm2_r
        end if
     endif

     residial=sqrt(norm2_r/norm2_r0)
     !        print*,'iter,residual',iter,residial

     if (norm2_r.le.tol*tol*norm2_r0) goto 10                !! wiki PBiCBSTAB 14

     iter=iter+1

     rho0=rho

  enddo
  

10 continue

end subroutine SPARSE_PBICGSTAB_nonull
end module linear_solver

