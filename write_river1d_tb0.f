      SUBROUTINE write_flow1d_tb0(un,fn,nfg,ng,no_signf)

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
     
!     rev. 9.1.64  Jul. 11/08  - NK: Coded up for tb0

!***********************************************************************
! - THIS SUBROUTINE OUTPUTS and r2c file.

! - List of arguments:

!   I - un      int        unit number
!   I - fn      int        file number
!   I - itogo   int        no. of hours until next rainfall
!   R - unit_conversion    REAL*4     conversion factor (area2)
!   I - FLN     CHAR*12    file names
!   I - DIR     CHAR*12    directory name
!   I - ocflg   int        open_close flag 1 open file -1 close file
!   I - frmflg  int        frame flag  1 new frame -1 end frame  
!                          0 each call = 1 frame with frame marks
!***********************************************************************

      use area_watflood
	implicit none

      INTEGER       :: un,fn,ii,no_signf,ng,nfg,i,j,l,ios
	character(20) :: junk
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday

!     FIRST TIME THROUGH THIS SUBROUTINE ONLY
!     OPEN OUTPUT FILE AND WRITE HEADER

      if(no_signf.lt.0)then

!     write the header

!     FILE NAMES AND UNIT NUMBERS DIFFER BY 30
!      write(*,1400)fn,fln(fn)
! 1400 format(' opening fln(',i3,'):',a30,'---')
!      write(*,*)

      open(unit=un,file=filename(fn),status='unknown',iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file',filename(fn)(1:40)
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        print*,'or target directory does not exist'
        stop 'Program aborted in write_flow_1d.f @ 48'
      endif
      PRINT*
	print*,'Opened unit=',un,' filename=',fln(fn)(1:40)
      write(un,3005)'########################################'
      write(un,3005)':FileType tb0  ASCII  EnSim 1.0         '
	write(un,3005)'# DataType               Time Series    '
      write(un,3005)'#                                       '
      write(un,3005)':Application             WATFLOOD       '
	write(un,3005)':Version                 2.1.23         '
	write(un,3020)':WrittenBy          ',author
      call date_and_time(cday,time)
	write(un,3010)':CreationDate       ',
     *     cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
3010  format(a20,a4,'-',a2,'-',a2,2x,a2,':',a2)
      write(un,3005)'#                                       '
	write(un,3005)'#---------------------------------------'
      write(un,3005)'#                                       '
      write(un,3020)'#SourceFile         ',source_file_name
	write(un,3020)':Name               ',name
      write(un,3005)'#                                       '
!     rev. 9.9.62  Mar.  21/15  - NK: Change zone from character to integer
c	write(un,3004)':Projection         ',coordsys_temp
c	if(coordsys_temp.eq.'LATLONG   ')then
c	  write(un,3004)':Ellipsoid          ',datum_temp
c	endif
c	if(coordsys_temp.eq.'UTM       ')then
c	  write(un,3004)':Ellipsoid          ',datum_temp
c        write(un,3001)':Zone               ',zone_temp
c	endif
      write(un,3004)':StartTime          ',startdate,starttime
	write(un,3005)':DeltaT  24:00:00.000   //1 day          '
c      write(un,3005)'#                                       '
c      write(un,3003)':UnitConverson      ',unit_conversion
      write(un,3005)'#                                       '
	write(un,3005)':ColumnMetaData                         '
c      write(un,3021)'   :ColumnUnits     ',column_units
	write(un,3021)'   :ColumnType      ',(column_type(l),l=1,ng)
	write(un,3021)'   :ColumnName      ',(resname(l),l=1,ng)
c	if(coordsys_temp.eq.'LATLONG   ')then
c	  write(un,3022)'   :ColumnLocationX ',(xres(l),l=1,ng)
c	  write(un,3022)'   :ColumnLocationY ',(yres(l),l=1,ng)
c	else
c	  write(un,3023)'   :ColumnLocationX ',(xres(l),l=1,ng)
c	  write(un,3023)'   :ColumnLocationY ',(yres(l),l=1,ng)
c	endif
      write(un,3005)':EndColumnMetaData                      '
      write(un,3005)'#                                       '
     	write(un,3005)':endHeader                              '

      return

      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fix fix

!       this is just to read in the new format
!       still needs to be programmed for multiple classes

!     the header is written only the first time through
!     this part is called with each class

!     :Frame and :EndFrame lines are written only for time series




!       NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE 
!       The r2c grids are written upside down: 
!                 south on top, north on bottom    !!!!!!!!!!
        if(no_signf.eq.0)then
          do i=1,ycount_temp
            write(un,1300)(outarray(i,j),j=1,xcount_temp)
1300        format(9999(1x,f5.0))
          end do          
        elseif(no_signf.eq.1)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1301)(outarray(i,j),j=1,xcount_temp)
1301        format(9999(1x,f5.1))
          end do          
        elseif(no_signf.eq.2)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1302)(outarray(i,j),j=1,xcount_temp)
1302        format(9999(1x,f5.2))
          end do          
        elseif(no_signf.eq.3)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1303)(outarray(i,j),j=1,xcount_temp)
1303        format(9999(1x,f6.3))
          end do          
        elseif(no_signf.eq.4)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1304)(outarray(i,j),j=1,xcount_temp)
1304        format(9999(1x,f7.4))
          end do          
        elseif(no_signf.eq.5)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1305)(outarray(i,j),j=1,xcount_temp)
1305        format(9999(1x,f8.5))
          end do          
        elseif(no_signf.eq.6)then   ! swe, precip and temp
          do i=1,ycount_temp
            write(un,1306)(outarray(i,j),j=1,xcount_temp)
1306        format(9999(1x,f9.6))
          end do          
        elseif(no_signf.eq.8)then   ! flow
          do i=1,ycount_temp
            write(un,1308)(outarray(i,j),j=1,xcount_temp)
1308        format(9999(1x,e10.3))
          end do          
        elseif(no_signf.eq.9)then   ! drain, dsnow
          do i=1,ycount_temp
            write(un,1309)(outarray(i,j),j=1,xcount_temp)
1309        format(9999(1x,f6.2))
          end do          
        else                        ! init soil moisture
          do i=1,ycount_temp
            write(un,1307)(outarray(i,j),j=1,xcount_temp)
1307        format(9999(1x,e12.6))
          end do          
        endif


        close(unit=un,status='keep')
        write(51,*)'Closed unit ',un,' Filename=  ',fln(fn)
        write(*,*)'Closed unit ',un,' Filename=  ',fln(fn)

      return

! FORMATS
 3000 format(a10,i5)
 3001 format(a20,i16)
 3002 format(2a20)
 3003 format(a20,f16.7)
 3004 format(a20,a10,2x,2a10)
 3005 format(a40)
 3006 format(a3,a10)
 3007 format(a14,i5,a6,i5)
 3012 format(a9)
 3020 format(a20,a40)
 3021 format(a20,999(1x,a11))
 3022 format(a20,999f10.0)
 3023 format(a20,999f10.4)
     
       END SUBROUTINE write_flow1d_tb0









