module neighbour_finding
  use common_parameter
  use common_2d
  use sphtools
  implicit none

  private
  public :: particle_allocation

  integer(ikind), dimension(:), allocatable ::ist, ip, nc
  integer(ikind), dimension(:), allocatable :: cellpart 

contains

  subroutine particle_allocation(x_mint,x_maxt,y_mint,y_maxt,phase_flag)

    integer ncx_temp,ncy_temp,nct_temp
    integer(ikind) :: phase_flag

    double precision x_mint,x_maxt,y_mint,y_maxt,   &
         &                  xb_temp_max,xb_temp_min,       &
         &                  yb_temp_max,yb_temp_min
     

    ! allocate linking lists
    allocate(ij_count(npfb))
    allocate(ij_link(npfb,nplink))
    n_count=npfb+nmirror_esti
    
    ! preamble
    ncx_temp=ncx
    ncy_temp=ncy
    nct_temp=nct+ncx+ncy+3
    xb_temp_max=x_maxt
    xb_temp_min=x_mint
    yb_temp_max=y_maxt
    yb_temp_min=y_mint
    allocate(ist(nct_temp+1))
    allocate(nc(nct_temp))
    allocate(ip(n_count))
    nc=0

    ! put particles in different cells, including mirror particles	
    call divide(1,np,x_mint,y_mint) 
!    call neighboring_list
    call neighboring_list_parallel

    deallocate(ist)
    deallocate(nc)
    deallocate(ip)
    deallocate(cellpart)

  end subroutine particle_allocation

! -----------------------------------------------------------------------
  subroutine divide(nini,nend,xmin_t,ymin_t) 
    !!args
    integer(ikind), intent(in) :: nini,nend 
    real(rkind), intent(in) :: xmin_t,ymin_t


    !!locals
    real(rkind) :: ddx,ddy
    integer(ikind) :: nc_ii_max, k, icell, jcell, ii

    nc=0_ikind  ! nc(i,mkind)=0

    nc_ii_max = 0
    allocate(cellpart(nend-nini+1))

!    !$omp parallel do private(ddx,ddy,icell,jcell,ii) reduction(+:nc)
    do k=nini, nend    ! loop over specified range of particles
       if(inbin(k)) cycle  ! don't try and find a cell for particles in the bin.
       ddx = r(k,1) - xmin_t    ! find x,y of particle if origin at BL corner of domain
       ddy = r(k,2) - ymin_t


       icell = int( ddx *uno_sup_size) + 1       !! Cell label in x
       jcell = int( ddy *uno_sup_size) + 1       !! Cell label in y
       !! Cell counter ii=[1,nct_t]        !! ii is the linear cell position in the matrix of cells
       ii = icell + (jcell - 1)*ncx        !! nc is the number of particles in cell ii
       if(ii .gt. nct .or. ii .lt. 1) then
          print *, "particles are out of domain"
          print *, "ii=", ii, "nct_t=", nct
          print *, "nct=", nct
          print *, "particle in trouble is ",k,ddx,ddy
          stop
       endif
       nc(ii) = nc(ii)+1     ! add 1 to the number of particles in this cell
       !! If the number of particles in a cell is greater than
       !! some upper limit then stop
       if(nc(ii).gt.nc_ii_max) nc_ii_max = nc(ii)  ! ensure nc_ii_max is max(nc(ii) for all ii
       if(nc_ii_max.gt.nplink)then       !! nplink is the limit for #particles in a cell
          print*,'itime',itime
          print*,'nc_ii_max',nc_ii_max
          print*,'nplink',nplink
          print*,'nc_ii_max.gt.nplink'
          print*,'k',k
          print*,'ii',ii
          print*,'icell',icell
          print*,'jcell',jcell
          stop
       end if
       cellpart(k)=ii             !! the index of the cell in which particle i is
       !!     ibox(ii,nc(ii))=k
       ! ibox tells us the following:
       ! Particle k is in cell ii, and its the 'nc(ii)'th particle in cell ii
    end do
!    !$omp end parallel do

    !! loop over all cells
    ist(1)=1_ikind
    do ii=1, nct
       ist(ii+1)=ist(ii)+nc(ii)   ! ist(ii) is the starting index for each cell
       nc(ii)=0_ikind             ! erase the nc in each cell here 
    end do

    !! ip labels particle k(=ibox) based on its cell position 
    !! and the number of particles in that cell. The result is greater 
    !! efficencies when constructing mirror parts. 
    do k=nini,nend
       if(inbin(k))cycle  ! don't do this for particles in the bin
       ip(ist(cellpart(k))+nc(cellpart(k)))=k   
        nc(cellpart(k))=nc(cellpart(k))+1
    end do 

  end subroutine divide

!---------------------------------------------------------------------------
  subroutine  neighboring_list
    use common_parameter
    use common_2d
    use sphtools

    implicit none

    integer j1
    integer lx,ly,lx2,ly2

    ij_count=0
    ij_link=0 !neighbour list arrays

    do ly=1,ncy !Loop over all cells in y direction
       do lx=1,ncx !Loop over all cells in x direction
          j1 = lx + (ly-1)*ncx !Cell index   

          if(nc(j1).gt.0) then

! if the cell is not empty, then loop over it and over neighboring cells

             !      Cells in the same XY sheet

             lx2 = lx+1
             if(lx2.le.ncx)then  
                call neighbor_list_ij(j1,j1+1)         !East cell
             end if
             ly2 = ly+1
             if(ly2.le.ncy)then
                call neighbor_list_ij(j1,j1+ncx)       !North cell
             end if
             lx2 = lx+1
             ly2 = ly+1
             if(lx2.le.ncx .and. ly2.le.ncy)then
                call neighbor_list_ij(j1,j1+ncx+1)     !NE cell
             end if
             lx2 = lx-1
             ly2 = ly+1
             if(lx2.ge.1 .and. ly2.le.ncy)then
                call neighbor_list_ij(j1,j1+ncx-1)     !NW cell
             end if
          endif
       enddo
    enddo

    ! Current cell
    do j1=1,nsheet_m                  
       call neighbor_list_self(j1)
    enddo

  end subroutine neighboring_list

  !------------------------------------------------

  subroutine neighbor_list_ij(j1,j2)

    use common_parameter
    use common_2d
    use sphtools

    implicit none

    integer j1,j2,is1,ie1,is2,ie2,i,ii,j,jj,k
    double precision temp,rr2tmp


       if(nc(j2)==0 .or. nc(j1)==0)return        !If cell is empty, return
       is1 = ist(j1)
       ie1 = ist(j1+1)-1
       do ii=is1,ie1
          i = ip(ii)
          is2 = ist(j2)
          ie2 = ist(j2+1)-1
          do jj=is2,ie2
             j = ip(jj)
             rr2tmp=dot_product(r(i,:)-r(j,:),r(i,:)-r(j,:))
             if(rr2tmp.gt.0.d0 .and. rr2tmp .le.sup_size_2)then  !large kernel h 	!sup_size_2 = 9 h^2
                ! If particles within interaction range make the linking list
                if(i.le.npfb )then
                   ij_count(i)=ij_count(i)+1;k=ij_count(i)
                   ij_link(i,k)=j  !the k'th interaction of particle i is with j
                end if
                if(j.le.npfb) then
                   ij_count(j)=ij_count(j)+1;k=ij_count(j)
                   ij_link(j,k)=i
                endif
             endif
          enddo
       enddo
    return
  end subroutine neighbor_list_ij

  !------------------------------------------------

  subroutine neighbor_list_self(j1)
    use common_parameter
    use common_2d
    use sphtools

    implicit none

    integer j1,is1,ie1,is1n,ie1n,i,ii,j,jj,k

    real(rkind) :: rr2tmp

       if(nc(j1)==0)return  !!if the cell is empty, move on.
       is1 = ist(j1)
       ie1 = ist(j1+1)-1-1 !ist(j1+1)-1
       do ii=is1,ie1   ! loop through the particles in the cell for i
          i = ip(ii)
          is1n = ii+1!ist(j1)
          ie1n = ist(j1+1)-1
          do jj=is1n,ie1n ! loop through the particles in the cell for j
             j = ip(jj)
             rr2tmp=dot_product(r(i,:)-r(j,:),r(i,:)-r(j,:))
             if(rr2tmp .gt. 0.d0 .and. rr2tmp .le. sup_size_2)then  !large kernel h
                if(i.le.npfb ) then
                   ij_count(i)=ij_count(i)+1;k=ij_count(i)
                   ij_link(i,k)=j
                endif
                if(j.le.npfb ) then
                   ij_count(j)=ij_count(j)+1;k=ij_count(j)
                   ij_link(j,k)=i
                endif
             endif
          enddo
       enddo
    return
  end subroutine neighbor_list_self

! -------------------------------------------------
  subroutine  neighboring_list_parallel
    use common_parameter
    use common_2d
    use sphtools

    implicit none

    integer(ikind) :: i,ic,icx,icy
    integer(ikind) :: jc

    ij_count=0
    ij_link=0 !neighbour list arrays
    
    
    !! Loop over all particles (in parallel)
    !$omp parallel do private(ic,icx,icy,jc)
    do i=1,npfb 
       if(.not.inbin(i))then
       !! Which cell are we in?
       ic = cellpart(i)   
       icx = mod(ic,ncx)      !! Position in row of cells
       icy = (ic/ncx) +1   !! Which row...
          
       !! This cell
       jc = ic
       call neighbour_particle_cell(i,r(i,:),jc)       
       
       if(icx+1.le.ncx)then
          !! East cell
          jc = ic+1
          call neighbour_particle_cell(i,r(i,:),jc)       

          if(icy+1.le.ncy)then
             !! NE cell
             jc = ic + ncx + 1    
             call neighbour_particle_cell(i,r(i,:),jc)       
          end if
          
          if(icy-1.ge.1)then
             !! SE cell
             jc = ic - ncx + 1
             call neighbour_particle_cell(i,r(i,:),jc)                    
          end if
       end if
       
       if(icx-1.ge.1)then
          !! West cell
          jc = ic - 1
          call neighbour_particle_cell(i,r(i,:),jc)       

          if(icy+1.le.ncy)then
             !! NW cell
             jc = ic + ncx - 1
             call neighbour_particle_cell(i,r(i,:),jc)                              
          end if
          
          if(icy-1.ge.1)then
             !! SW cell
             jc = ic - ncx - 1
             call neighbour_particle_cell(i,r(i,:),jc)                              
          end if         
       end if
       
       if(icy+1.le.ncy) then
          !! North cell
          jc = ic + ncx
          call neighbour_particle_cell(i,r(i,:),jc)                              
       end if
       
       if(icy-1.ge.1)then
          !! South cell
          jc = ic - ncx
          call neighbour_particle_cell(i,r(i,:),jc)                              
       end if
       end if
    end do
    !$omp end parallel do
!    stop
  end subroutine neighboring_list_parallel

  subroutine neighbour_particle_cell(ii,ri,jc)
    integer(ikind),intent(in) :: ii,jc !! jc is cell
    real(rkind),dimension(:),intent(in) :: ri
    real(rkind) :: rr2tmp,qq,rad
    integer(ikind) :: j,jj,is,ie,kk
    real(rkind),dimension(dims) :: rij,gradw
    
    if(nc(jc).ne.0)then !! if the cell isn't empty
       is =ist(jc);ie=ist(jc+1)-1  !! Start and end indices of particles in cell jc
       do jj=is,ie       !! Loop over all particles in cell jc
          j=ip(jj)       !! j is regular index, jj is cell-ordered index

          rij(:) = ri(:)-r(j,:)
          rr2tmp=dot_product(rij,rij)  !! Distance squared
          if(rr2tmp .gt. 0.d0 .and. rr2tmp .le. sup_size_2)then
!          if(rr2tmp .le. sup_size_2)then              
             ij_count(ii)=ij_count(ii)+1;              !! Increment count
             kk=ij_count(ii)
             ij_link(ii,kk)=j                !! add to list  
          endif          
       end do
    end if
  
    return
  end subroutine neighbour_particle_cell


end module neighbour_finding

!--------------------------------------

