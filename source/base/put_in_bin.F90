subroutine put_in_bin(x_mint,x_maxt,y_mint,y_maxt)
  use kind_parameters
  use common_parameter
  use common_2d

  ! subroutine to put any particle out of the domain in the bin.
  ! rest of program should then skip it when updating everything.

  !arguments
  real(rkind) :: x_mint,x_maxt,y_mint,y_maxt

  ! vars
  integer(ikind) :: i
  ! loop over all fluid particles
  !$omp parallel do
  do i=nbS+1,npfb
     ! is the particle already in the bin?
     if(.not.inbin(i)) then
        ! is the particle out of the domain?
        if(r(i,1).lt.x_mint.or.r(i,1).gt.x_maxt.or.r(i,2).lt.y_mint.or.r(i,2).gt.y_maxt)then
           inbin(i) = .true.
        end if
        ! is the particle position NaN
        if(isnan(r(i,1)).or.isnan(r(i,2)))then
           inbin(i) = .true.
        end if
     end if
  end do
  !$omp end parallel do

  ! now loop to correct position of everything in the bin
  n_inbin = 0
  n_free = 0
  !$omp parallel do
  do i=nbS+1,npfb
     if(inbin(i))then
        ! re-locate the particle to somewhere out of the way...
        r(i,1) = 2.0*x_mint - x_maxt
        r(i,2) = 2.0*y_mint - y_maxt
        ! give it a zero velocity
        u(i,:) = 0.0_rkind
        n_inbin = n_inbin + 1
        !! Put the particle in the free particles array
        n_free = n_free + 1   ! add 1 to the number of free indices for gas particles
        p_free(n_free) = i    ! store the free index in the list of free indices
     end if
  end do
  !$omp end parallel do
  
  !! Check whether everything is in the bin...
  if(n_inbin.eq.npfb-nbS)then
     write(6,*) "All ",npfb - nbS," liquid particles are in the bin"
     stop
  end if
  return
end subroutine put_in_bin
