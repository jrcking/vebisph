
! --------------------------------------------------------------------

subroutine output
  use common_parameter
  use common_2d
  implicit none

  integer idp,nf,i,n_out
  character(70) fnm1

  ! File naming stuff
  nf = 22
  idp = n_dump

  ! STEP 1: Output particle data
!  if(mod(itime,10).eq.0)then   ! If it is time to output particle data
  if(time.gt.next_dump.or.output_everytime) then

     next_dump = next_dump + dt_out
     n_dump = n_dump + 1 

     ! Create the name of the output file
     if( idp .lt. 10 ) then 
        write(fnm1,'(A21,I1)') './data_directory/PART', idp
     else if( idp .lt. 100 ) then 
        write(fnm1,'(A21,I2)') './data_directory/PART', idp
     else if( idp .lt. 1000 ) then
        write(fnm1,'(A21,I3)') './data_directory/PART', idp
     else
        write(fnm1,'(A21,I4)') './data_directory/PART', idp
     end if

     ! Open the file
     open(nf+1,file=fnm1,status='unknown')

     ! Will we output mirror particles?
     if(output_mirrors) then
        n_out=npfb + nmirror
     else
        n_out=npfb 
     endif
     ! Output the data
     do i=1,n_out
        if(inbin(i)) cycle  ! don't waste space outputting parts in bin
#if const_mod==1        
        write(nf+1,*) r(i,1),r(i,2),u(i,1),u(i,2),theta(i),P(i),conc(i),a_out(i),n_surf(i),vort(i)
#else
        write(nf+1,*) r(i,1),r(i,2),u(i,1),u(i,2),theta(i),P(i),conc(i),a_out(i),n_surf(i),vort(i) &
                      ,tau_p(i,1,1),tau_p(i,1,2),tau_p(i,2,2)        
#endif        
     end do
     close(nf+1)
  
     !! write and push some time data
     write(94,*) time,n_out,itime,dt
     flush(94)

     ! STEP 2: Dump output for restarts
     open(100,file='TEMP_REC.dat')
     write(100,*)itime,time,n_dump
     write(100,*)np-n_inbin,npfb-n_inbin,nbS,np_inflow
     write(100,*)xb_min,xb_max,yb_min,yb_max, dx, av_conc


     do i=1,npfb
        if(inbin(i)) cycle
        write(100,*) r(i,1),r(i,2),u(i,1),u(i,2),theta(i),p(i),dv,n_surf(i)
#if const_mod!=1
           write(100,*) tau_p(i,1,1),tau_p(i,1,2),tau_p(i,2,1),tau_p(i,2,2)
#endif           
     enddo
     !! Indicate the inflow particles
     if(np_inflow.ne.0)then
        do i=1,np_inflow
           write(100,*) p_inflow(i)
        end do
     end if

     close(100)
  end if

  ! STEP 3: Output time data every step
  write(92,*) time,dt,npfb-n_inbin,maxval(abs(u(1:npfb,1))),maxval(abs(u(1:npfb,2))),maxval(abs(P(1:npfb)))
  write(93,*) LS_iters,LS_residual


  return
end subroutine output
