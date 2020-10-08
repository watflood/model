	subroutine precip_adjust()

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen 
        
!    This file is part of WATFLOOD (R)      
        
!    WATFLOOD(R) is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    WATFLOOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.

!    You should have received a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
     
!     this s/r created the precip adjustment factors
!     These are used in process_rain.f

!     rev. 9.5.69  Oct.  10/09  - NK: added xcount & ycount to error & paf files
!     rev. 9.8.30  Oct.  16/12  - NK: remove p(i,j)=0.0 from precip_adjust

      use area_watflood
	implicit none

	logical   ::    exists,newpafflg
	integer   ::    i,j,ios,block,ig,jg,iAllocate
      real*4, dimension(:,:),   allocatable :: grid_error
      real*4    ::    conv
      character*1 ::  answer
      
      allocate(grid_error(ycount,xcount),precadj(ycount,xcount),
     *      stat=iAllocate)
      if(iopt.eq.2)print*,'allocation done in precip_adjust '
      if(iopt.eq.2)print*
 
      nblock=0
      do i=ycount,1,-1
        do j=1,xcount
          precadj(i,j)=1.0
          grid_error=0.0
!     rev. 9.8.30  Oct.  16/12  - NK: remove p(i,j)=0.0 from precip_adjust
c          p(i,j)=0.0
        end do
      end do

      precflg=.false.
      INQUIRE(FILE='paf.r2s',EXIST=exists)
!     LOOK FOR A PAF.TXT FILE
!     IF IT IS NOT THERE, LOOK FOR AN ERROR.TXT FILE
      IF(exists)THEN
        precflg=.true.
        newpafflg=.true.
!       unit=99
        fln(99)='paf.r2s'
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call read_r2s_beta(99,fln(99))
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  write(51,*)'Precip adjustment factors  PAF'
        do i=1,ycount,1
	    do j=1,xcount
            precadj(i,j)=amin1(2.0,inarray(i,j))
            precadj(i,j)=amax1(0.5,inarray(i,j))
	    end do
          write(51,201)(precadj(i,j),j=1,xcount)
	  end do
        print*,'read & closed the paf.r2s file'
        precflg=.true.

      else   ! paf.txt does not exist
!       look for an error.txt file instead
!       IF THERE IS NO PAF.r2s FILE, LOOK FOR A ERROR.txt FILE
!       AND IF IT'S NOT THERE EITHER, PASS OVER THE ADJUSTMENT SECTION
        INQUIRE(FILE='error.txt',EXIST=exists)
        IF(exists)THEN
          open(99,file='error.txt',form='formatted',iostat=ios)
	    if(ios.ne.0)then
	      print*,'Problems opening error.txt file'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
 	      print*
	      stop 'Program aborted in precip_adjust @ 481'
          else
              print*,'++++++++++++++++++++++++++++'
              print*,'+     Opened error.txt     +'
              print*,'++++++++++++++++++++++++++++'
          endif
          if(iopt.eq.2)print*,'@ 2 in precip_adjust @ 481'
          precflg=.true.
          newpafflg=.true.
          write(51,202)precflg
	    write(51,*)'Grid error'
          read(99,203,iostat=ios)nblock
	    if(ios.ne.0)then
	      print*,'Probelms reading the first line in error.txt'
	      print*
	      stop 'Program aborted in precip_adjust @ 112'
	    endif
	    write(51,*)'nblock= ',nblock
	    
          do while(nblock.gt.0)
            block=float(nblock)
            write(51,204)nblock
            print*,'reading error.txt file'
 !     rev. 9.5.69  Oct.  10/09  - NK: added xcount & ycount to error & paf files
!           Note:  the precip data may be in mc2 format - 99 fields
            read(99,*)   ! xcount - not needed
	      read(99,*)   ! ycount - not needed
            do i=ycount,1,-1
!             this is not a GK format so file is north up            
c              read(99,*,iostat=ios)(p(i,j),j=1,xcount)
              read(99,*,iostat=ios)(grid_error(i,j),j=1,xcount)
              if(ios.ne.0)then
	          print*,'Problems reading data line ',i,' in error.txt'
	          print*
	          stop 'Program aborted in precip_adjust @ 112'
	        endif
            end do
            write(*,*)'Reading block',block
           do i=1,ycount
              write(51,2000)(grid_error(i,j),j=1,xcount)
              do j=1,xcount
c               grid_error(i,j)=grid_error(i,j)+p(i,j)/float(1+nblock**2)    !1.5
                precadj(i,j)=precadj(i,j)*block/
     *                   (block+grid_error(i,j)/100.0)
              end do
            end do
            write(51,*)'PAF - precipitation adjustment factor'
            do i= ycount,1,-1
                write(51,2001)(precadj(i,j),j=1,xcount)
            end do          
            read(99,203,iostat=ios)nblock      
            if(ios.ne.0)then
!             IF WE HAVE REACHED THE END OF THE FILE, OK 
!             THERE SHOULD BE A BLOCK INDEX NUMBER IF THERE 
!             ANOTHER DATA BLOCK
	        print*,'Problems reading data line ',i,' in error.txt'
	        print*
	        stop 'Program aborted in precip_adjust @ 112'
	      endif
          end do
          
!         SET LIMITS ON THE PRECIPITATION ADJUSTMENT FACTORS:
          do i=1,ycount
            do j=1,xcount
              precadj(i,j)=amin1(2.0,precadj(i,j))
              precadj(i,j)=amax1(0.5,precadj(i,j))
            end do
          end do
          if(iopt.eq.2)print*, '2 in rain'
          close(99)

!         Write a new PAF.r2s file
          author='CHARM***.exe'     
          name='PrecipAdjustmentFactor'
          coordsys_temp=coordsys1
!         GreenKenue uses LatLong - code below uses LATLONG
          if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
          if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
          zone_temp=zone1
          datum_temp=datum1
          xorigin_temp=xorigin
          yorigin_temp=yorigin
          xcount_temp=xcount
          ycount_temp=ycount
          xdelta_temp=xdelta
          ydelta_temp=ydelta
          attribute_name='PAF                                     '
          attribute_units='1                                      ' 
          attribute_type='None                                    '  
          source_file_name='error.txt'
          unit_conversion=1.0     
          fln(99)='newpaf.r2s'  
          do i=1,ycount
            do j=1,xcount
              outarray(i,j)=precadj(i,j)
            end do
          end do
          call write_r2s(99,fln(99))
        endif
      endif
      
d      do i=ycount,1,-1
d        write(*,201)(precadj(i,j),j=1,xcount)
d      end do
      
      return

!     FORMATS
 2000 format(<xcount>f5.0)
  201 format(<xcount>f5.2)
 2001 format(<xcount>f5.2)        !for mc2 files for map
  202 format(' Precflg= ',l1,'  Precipitation adjustment factors:')
  203 format(i5)
  204 format(' nblock=',i5)
  205 format(3i5,f10.2)
  206 format('    l    i    j       PAF  ycount,xcount:',2i5)   
  207 format(' Please note: limits of 0.5 to 2.0 set on PAF`s')
  208 format(10x,l1)
  217 format(6x,2i5,5x,a5,a30,i5)
  218 format(6x,2i5,5x,a5,a30)
  220 format(10x,i10)

	end subroutine precip_adjust
