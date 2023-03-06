module velocity_correction
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  implicit none

  public :: correct_velocity

  !! NOT OMP YET !!

contains
  subroutine correct_velocity
    integer(ikind) :: i, j, k, num
    real(rkind) :: temp1,u_norm
    integer(ikind), dimension(:), allocatable :: imirror_close

    real(rkind) :: r_min2
    real(rkind),dimension(dims) :: uij,wall_norm

    allocate(imirror_close(npfb))
    imirror_close=0_ikind

    ! Find the particles close to wall mirrors
    do i=1+nbS, npfb
       r_min2=dx*dx ! cut-off distance is dx
       do k=1,ij_count(i)
          j=ij_link(i,k)

          !! For particles close to wall mirrors, mark them as such
          if(j.gt.npfb .and. dot_product(r(i,:)-r(j,:),r(i,:)-r(j,:)).lt. r_min2) then
             if(vrelation(j).gt.0.and.b_type(vrelation(j)).eq.1) then
                r_min2=dot_product(r(i,:)-r(j,:),r(i,:)-r(j,:))
                imirror_close(i)=j
             end if
             if(vrelation(j).eq.-1.or.vrelation(j).lt.-2)then  ! corners
                r_min2=dot_product(r(i,:)-r(j,:),r(i,:)-r(j,:))
                imirror_close(i)=j
             end if
          end if
       enddo
    enddo

    ! Correct the velocity
    do i=1+nbS,npfb    ! for each particle i
       if(imirror_close(i).ne.0)then   ! if it is close to a wall mirror
          j=imirror_close(i)

	  ! calculate the unit vector between i and the wall mirror
          wall_norm(:) = r(i,:)-r(j,:)
          temp1=sqrt(dot_product(wall_norm,wall_norm))
          wall_norm(:) = wall_norm(:)/temp1

          !! relative velocity
          uij(:) = u(i,:) - u(j,:)

	  ! calculate velocity component towards the mirror
          u_norm=dot_product(uij,wall_norm)  

          !! remove particle velocity component toward mirror
          if(u_norm.le.0.0d0)then
             u(i,:)=u(i,:) - u_norm*wall_norm(:)
          endif
       endif
    enddo

    deallocate(imirror_close)


  end subroutine correct_velocity

end module velocity_correction
