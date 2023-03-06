!! Gridding cells using preconditioned dynamic vector.
!! See J.M. Dominguez, 2010 
module cell_griddin
  use kind_parameters
  use common_parameter
  !  use common_2d
  implicit none

  public :: cell_generator, cell_position

contains

  subroutine cell_generator(xmin,xmax,ymin,ymax,uno_sup_size,ncx_t,ncy_t,nct_t)

    real(rkind), intent(in) ::  xmin,xmax,ymin,ymax,uno_sup_size
    integer(ikind), intent(out) :: ncx_t,ncy_t,nct_t

    !! determine number of cells in x, y directions and total
    ncx_t = int( (xmax-xmin)*uno_sup_size ) + 1 !rui 29/06
    ncy_t = int( (ymax-ymin)*uno_sup_size ) + 1

    nct_t=ncx_t*ncy_t ! Number of cells in a XY sheet

  end subroutine cell_generator

  subroutine cell_position(xmin,xmax,ymin,ymax,sup_size,xc,yc,ncx,ncy,nct)

    implicit none

    integer(ikind) :: i,j,k
    integer(ikind) :: ncx,ncy,nct
    real(rkind) :: xc(nct),yc(nct)
    real(rkind) :: xmin,xmax,ymin,ymax,sup_size
    real(rkind) :: temp_x,temp_y


    i=1
    temp_x=xmin+sup_size/2.0d0
    do while(i.le.ncx)

       j=1
       temp_y=ymin+sup_size/2.0d0
       do while(j.le.ncy)
          k=i+(j-1)*ncx
          xc(k)=temp_x
          yc(k)=temp_y
          j=j+1
          if(j.lt.ncy)then
             temp_y=temp_y+sup_size
          elseif(j.eq.ncy)then
             temp_y=ymax-ymin
             temp_y=temp_y-dble(ncy-1)*sup_size
             temp_y=ymax-temp_y/2.0d0
          endif
       enddo

       i=i+1
       if(i.lt.ncx)then
          temp_x=temp_x+sup_size
       elseif(i.eq.ncx)then
          temp_x=xmax-xmin
          temp_x=temp_x-dble(ncx-1)*sup_size
          temp_x=xmax-temp_x/2.0d0
       endif

    enddo

    return
  end subroutine cell_position
end module cell_griddin
