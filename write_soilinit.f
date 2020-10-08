      SUBROUTINE write_soilinit()

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
     
!     rev. 9.3.05  Nov.  13/06  - NK: adder write_flowinit.for to flowinit.for

!     "author" must be specifief before calling this s/r

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      logical    :: exists
c	integer(4)    :: i,j,n,ii,ios,mm
	integer(4)    :: i,j,n,ii,ios,k
      CHARACTER(30) :: sys_command
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
!     write the init file for watroute in r2c format

!     rev. 10.2.50 Mar.  22/19  - NK: New resume files written to \resume\*.*
	fln(99)='resume\soil_init.r2c'
      inquire(FILE=fln(99),EXIST=exists)
      open(99,file=fln(99),status='unknown',iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file  soil_init.r2c'
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        stop 'Program aborted in write_soilinit.f @ 33'
      endif


!     writing a new initial conditions r2c file for watroute      

!     write the init file for runoff in r2c format
!     write the init file for runoff in r2c format
!     write the init file for runoff in r2c format
!     write the init file for runoff in r2c format
!     write the init file for runoff in r2c format

!     Jan 22/07  nk

      write(99,3005)'########################################'
      write(99,3005)':FileType r2c  ASCII  EnSim 1.0         '
      write(99,3005)'#                                       '
	write(99,3005)'# DataType               2D Rect Cell   '
      write(99,3005)'#                                       '
      write(99,3005)':Application             EnSimHydrologic'
	write(99,3005)':Version                 2.1.23         '
	write(99,3005)':WrittenBy               soilinit       '
      call date_and_time(cday,time)
	write(99,3010)':CreationDate       ',
     *       cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
3010  format(a20,a4,'/',a2,'/',a2,2x,a2,':',a2)
      write(99,3005)'#                                       '
        if(tbcflg.eq.'y')then
	    write(99,3011)'#start y/m/d/h:',
     *                   year_start,mo_start,day_start,hour_start
!     REV. 10.1.26 Mar.  23/16  - NK: Fixed comment for spinup period 
	    write(99,3011)'#end   y/m/d/h:',
     *                   year_now,month_now,day_now,hour_now
 3011 format(a15,4i10)
        endif
	write(99,3005)'#---------------------------------------'
      write(99,3020)':SourceFileName     ',fln(6)
      write(99,3005)'#                                       '
	write(99,3004)':Projection         ',coordsys1
!     rev. 9.9.62  Mar.  21/15  - NK: Change zone from character to integer
	if(coordsys1.eq.'UTM       ')then
          write(99,3001)':Zone               ',zone1
	    write(99,3004)':Ellipsoid          ',datum1
	endif
	if(coordsys1.eq.'LATLONG   ')then
	    write(99,3004)':Ellipsoid          ',datum1
      endif
      write(99,3005)'#                                       '
      write(99,3003)':xOrigin            ',xorigin
      write(99,3003)':yOrigin            ',yorigin
      write(99,3005)'#                                       '
      n=0
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'v         ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'d1        ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'d1fs      ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'uzs       ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'uzsfs     ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'snowc     ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'sca       ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'wcl       ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'def       ' 
      end do
	do ii=1,classcount
	  n=n+1
        write(99,3008)':AttributeName ',n,ii,'api       ' 
      end do
	  n=n+1
        write(99,3009)':AttributeName ',n,'tto       ' 
	  n=n+1
        write(99,3009)':AttributeName ',n,'ttomin    ' 
	  n=n+1
        write(99,3009)':AttributeName ',n,'ttomax    ' 

      write(99,3005)'#                                       '
      write(99,3001)':xCount             ',xcount
      write(99,3001)':yCount             ',ycount
      write(99,3003)':xDelta             ',xdelta
      write(99,3003)':yDelta             ',ydelta
      write(99,3005)'#                                       '
      write(99,3005)':EndHeader                              '

 3000 format(a10,i5)
 3001 format(a20,i16)
 3002 format(2a20)
 3003 format(a20,f16.7)
 3004 format(a20,a10,2x,a10)
 3005 format(a40)
 3006 format(a3,a10)
 3007 format(a14,i5,a6,i5)
 3008 format(a15,2i5,a10)
 3009 format(a15,i5,5x,a10)
 3012 format(a9)
 3020 format(a20,a40)
 	if(iopt.eq.2)print*, 'in flowinit at 1302'

!     initialize p() to male sure there is no junk in the unused grids
      do i=1,ycount
	  do j=1,xcount
	    p(i,j)=0.0
	  end do
	end do
!     p(1,1) is used to show the top left corner of each frame in the file. 
!     If there is a value, it's over written.
      k=0
!     Initial interception mm storage
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=v(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial surface storage mm bare
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=d1(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial surface storage mm snow
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=d1fs(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial upper zone storage mm bare
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=uzs(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial upper zone storage mm snow
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=uzsfs(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial snow water equivalent in mm swe 
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=snowc(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial snow covered area fraction sca
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=sca(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial free water content in snow mm
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=wcl(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial heat deficit mm
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=def(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial api mm
      do ii=1,classcount
        k=k+1
        p(1,1)=float(k)
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=api(n,ii)
        end do
        do i=1,ycount
          write(99,4009)(p(i,j),j=1,xcount)
        end do
	end do

!     Initial tto  
      k=k+1
      p(1,1)=float(k)
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=tto(n)
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

!     Initial ttomin  
      k=k+1
      p(1,1)=float(k)
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=ttomin(n)
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

!     Initial ttomax 
      k=k+1
      p(1,1)=float(k)
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        p(i,j)=ttomax(n)
      end do
	do i=1,ycount
        write(99,4009)(p(i,j),j=1,xcount)
      end do

      close(unit=99,status='keep')

 4001 format(999(' ',f5.0))
 4002 format(999(' ',f10.3))
 4003 format(999(' ',i5))
 4004 format(999(' ',f10.5))
 4005 format(999(' ',f10.7))
 4006 format(999(' ',f8.1))
 4007 format(999(' ',f8.0))
 4008 format(999(' ',f8.3))
 4009 format(999(' ',e10.3))
 4010 format(999(' ',f5.3))
 4011 format(999(' ',f5.1))

 4050 format(999(i4))
 4052 format(999(' ',i3))

      close(unit=99,status='keep')

      return

      end SUBROUTINE write_soilinit

