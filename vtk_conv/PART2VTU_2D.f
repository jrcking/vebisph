c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr. Alejandro Crespo, Dr. Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
c
c    This file is part of SPHYSICS.
c
c    SPHYSICS is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 3 of the License, or
c    (at your option) any later version.
c
c    SPHYSICS is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.


c       initial data for SPH models                                                                     72

      program PART2VTU_2D

      parameter(np_max = 270000,i_PART_counter_max=20000)
      parameter(num_phase_max=9)
      character chartemp*40, name_orig*40
      character name_vtu*40, name_vtu2*12, name_vtu3*9
      character np_string3*3, np_string4*4, np_string5*5
      character np_string6*6, np_string7*7, np_string8*8
      character frame_string1*1, frame_string2*2, frame_string3*3
      character frame_string4*4, frame_string5*5, frame_string6*6
      character supp*4,supp3*3,supp2*2,supp1*1, zero_string
      character string1*100,string2*100,string3*100,string4*100
      character chartemp2*100
      CHARACTER(LEN=10) :: FMT,FMT1
      CHARACTER(LEN=1)  :: TAB,DQ
      
      real xp(np_max),zp(np_max),up(np_max),wp(np_max)
      real p(np_max),rho(np_max),alpha(np_max),vort(np_max)
      real theta(np_max)
      real time(i_PART_counter_max), DT(i_PART_counter_max)
      integer nsurf(np_max)
      integer np_all(i_PART_counter_max), IT(i_PART_counter_max)
      real  DT1(i_PART_counter_max),DT2(i_PART_counter_max)  
      
      real np_phase_start(num_phase_max),np_phase(num_phase_max)
      real rho0_phase(num_phase_max),P0_phase(num_phase_max) 
      real gamma_phase(num_phase_max),viscos_phase(num_phase_max)
      real ST_coeff(num_phase_max),backgroundPressure
          
      TAB=CHAR(9)     
      FMT="(A)"
      FMT1="(2A)"
      DQ=CHAR(34)
      
       print*
       write(*,*) '<SPHYSICS>  Copyright (C) <2007>'
       write(*,*) 
     & '<Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, '
       write(*,*) 
     & 'Dr Benedict Rogers, Alejandro Crespo, '
       write(*,*) 
     & 'Muthukumar Narayanaswamy, Dr Shan Zou, and Dr Andrea Panizzo >'
       write(*,*) 'This program comes with ABSOLUTELY NO WARRANTY;    '
       write(*,*) 'This is free software, and you are welcome to      '
       write(*,*) 'redistribute it under conditions stated in         '
       write(*,*) 'the GPL License;                                   '

      print*
      print*,' ---     PART2VTU_2D.F    ---'
      print*,' ---   GENERATING GEOMETRY   ---'
      print*,' ---   Distributed under     ---'
      print*,' ---    the GPL License      ---'
      print*

      
c     %LOAD AND READ TIME,DT FROM FILE DT. THE FIRST HEADERLINE IN 
        open(unit=70,file='../data_directory/TIME_OUT',status='old')
        i_loop_finish = 0
        i_PART_counter = 0
        do while(i_loop_finish.eq.0)
          i_PART_counter = i_PART_counter + 1
          if(i_PART_counter.gt.i_PART_counter_max)then
            print*,'Number of entries in file DT exceeds max value'
            print*,'i_PART_counter.gt.i_PART_counter_max'
            print*,'Adjust i_PART_counter_max, i_PART_counter_max = ',
     &              i_PART_counter_max
            stop
          endif
          read(70,*,END = 76)time(i_PART_counter),
     &np_all(i_PART_counter),IT(i_PART_counter),DT(i_PART_counter)
           
          !Determine whether to exit loop
          if(i_loop_finish.eq.0)then
            i_loop_finish = i_loop_finish - 1
          endif
76          i_loop_finish = i_loop_finish + 1
          !print*
        enddo
        N_start = 1
        Nframes = i_PART_counter-2  !Why -2?
      write(6,*) "There are ",Nframes,"frames."

      !! JRCK addition...
      write(6,*) "Enter starting frame"
      read(*,*) N_start
    
      ngrab = N_start-2
      do iframe=N_start,Nframes
              
c       % READ IN THE PART FILE FOR EACH FRAME
        ngrab=ngrab+1
        if(ngrab.lt.10) then
           write(supp1,'(i0)') ngrab
           name_orig='../data_directory/PART'//supp1
        end if
        if(ngrab.ge.10.and.ngrab.lt.100) then
           write(supp2,'(i2)') ngrab
           name_orig='../data_directory/PART'//supp2
        end if
        if(ngrab.ge.100.and.ngrab.lt.1000) then
           write(supp3,'(i3)') ngrab
           name_orig='../data_directory/PART'//supp3
        end if
        if(ngrab.ge.1000.and.ngrab.lt.10000) then
           write(supp,'(i4)') ngrab
           name_orig='../data_directory/PART'//supp
        end if
        write(supp,'(i4.4)') ngrab

        name_vtu ='../paraview_files/PART'//supp//'.vtu'
        print*, 'iframe, name_orig, name_vtu ',
     &              iframe, ' ',name_orig, name_vtu 
        
        open(23,file=name_orig,status='old')
        open(24,file=name_vtu,status='unknown')
c       % READ POSITION, VELOCITY, DENSITY, PRESSURE, MASS AND VORTICITY DATA FOR ALL PARTICLES                      
         npp=0 ! keeps track of number of particles
         np = np_all(iframe+1)
         do i=1,np
             read(23,*,end=300) xp(i),zp(i),up(i),wp(i),theta(i),p(i),rho(i)
     &,alpha(i),nsurf(i),vort(i)
            npp=npp+1
         enddo
300    np=npp

        close (23)
        print*,'np ',np
                                                                
              
       
201     format(a40)
202     format(a100)
203     format(a25,i7,a17,i7,a2)
211     format(a21)
c     % OUTPUT TO FILE IN VTU FORMAT 
        if(np.lt.1000)then       
          write(np_string3,'(i3.3)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string3//DQ//' Numb
     &erOfCells='//DQ//np_string3//DQ//'>'
        elseif(np.lt.10000)then       
          write(np_string4,'(i4.4)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string4//DQ//' Numb
     &erOfCells='//DQ//np_string4//DQ//'>'
        elseif(np.lt.100000)then       
          write(np_string5,'(i5.5)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string5//DQ//' Numb
     &erOfCells='//DQ//np_string5//DQ//'>'
        elseif(np.lt.1000000)then       
          write(np_string6,'(i6.6)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string6//DQ//' Numb
     &erOfCells='//DQ//np_string6//DQ//'>'
        elseif(np.lt.10000000)then       
          write(np_string7,'(i7.7)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string7//DQ//' Numb
     &erOfCells='//DQ//np_string7//DQ//'>'
        elseif(np.lt.100000000)then       
          write(np_string8,'(i8.8)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string8//DQ//' Numb
     &erOfCells='//DQ//np_string8//DQ//'>'
        else
          print*,'Too many particles for np_string'
          stop  
        endif
        !print*,'np_string, np ',np_string, np 
        string1 = '<?xml version='//DQ//'1.0'//DQ//'?>'
        string2 = '<VTKFile type= '//DQ//'UnstructuredGrid'//DQ//'  vers
     &ion= '//DQ//'0.1'//DQ//'  byte_order= '//DQ//'BigEndian'//DQ//'>'
        string3 = ' <UnstructuredGrid>'
        write(24,211)string1
        write(24,202)string2
        write(24,202)string3
        write(24,202)string4
              
c       % WRITE IN PRESSURE DATA
        string1 = '   <PointData Scalars='//DQ//'Pressure'//DQ//' Vector
     &s='//DQ//'Velocity'//DQ//'>'
        string2 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Pressures'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        write(24,202)string2
        do ii=1,np
          write(24,*)p(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202)string3

c       % WRITE Temperature DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Theta'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        do ii=1,np
          write(24,*)theta(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3

c       % WRITE density DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Particle Density'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        do ii=1,np
          write(24,*)rho(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3

c       % WRITE alpha DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Alpha'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        do ii=1,np
          write(24,*)alpha(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3

c       % WRITE vort DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'w'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        do ii=1,np
          write(24,*)vort(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3

c       % WRITE U-VELOCITY DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'u'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        do ii=1,np
          write(24,*)up(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3

c       % WRITE W_VELOCITY DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'v'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        do ii=1,np
          write(24,*)wp(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3

c       % WRITE surface DATA
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//D
     &Q//'surface'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202)string1
        do ii=1,np
          write(24,*)nsurf(ii)
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3
              
c       % WRITE VELOCITY DATA
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Velocity'//DQ//' NumberOfComponents='//DQ//'3'//DQ//' format='
     &//DQ//'ascii'//DQ//'>'
        write(24,202) string1
        do ii=1,np
          write(24,*)up(ii),wp(ii),0.0
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3
        string4 = '   </PointData>'
        write(24,202) string4
              
c       % WRITE PARTICLE POSITION DATA
        string2 = '   <Points>'
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' NumberOfCo
     &omponents='//DQ//'3'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202) string2
        write(24,202) string1
        do ii=1,np
          write(24,*)xp(ii),zp(ii),0.0
        enddo
        string3 = '    </DataArray>'
        string2 = '   </Points>'
        write(24,202) string3
        write(24,202) string2
             
c       % WRITE CELL DATA. CELL IS OF TYPE VERTEX.        
        string2 = '   <Cells>'
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ/
     &/'connectivity'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202) string2
        write(24,202) string1
        do ii=1,np
          write(24,*)ii-1
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3
        
        
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ/
     &/'offsets'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202) string1
        do ii=1,np
          write(24,*)ii
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3
        
        
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ/
     &/'types'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(24,202) string1
        do ii=1,np
          write(24,*)1
        enddo
        string3 = '    </DataArray>'
        write(24,202) string3
        
        
        string1 = '   </Cells>' 
        string2 = '  </Piece>'
        string3 = ' </UnstructuredGrid>'
        string4 = '</VTKFile>'
        write(24,202) string1
        write(24,202) string2
        write(24,202) string3
        write(24,202) string4
      enddo
      close(24);
      stop
      end
