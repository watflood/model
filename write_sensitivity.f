      subroutine write_sens_result(rvrname,classname,delta_factor)

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
     
      use area_watflood
      USE EF_module
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      character(20) :: rvrname(12),classname(24)
      integer   :: ii,j,ios
      real*4    :: delta_factor


      open(unit=99,file='sensitivity.txt',status='unknown',
     *             iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file  sensitivity.txt'
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        stop 'Program aborted in write_sensitivity.f @ 23'
      endif

      write(99,*)'Routing parameters:'
      write(99,10000)'param ',(rivtype(ii),ii=1,nrvr)
      do j=1,12,2
        write(99,10000)rvrname(j)
        write(99,10001)int(delta_factor),(nrvr_array(ii,j),ii=1,nrvr)
        write(99,10002)int(delta_factor),(nrvr_array(ii,j+1),ii=1,nrvr)
      end do
      write(99,*)'Hydrological parameters:'
      write(99,10010)'param ',(nclass(ii),ii=1,classcount)
      do j=1,24,2
        if(j.ne.17)then
          write(99,*)classname(j)
          write(99,10003)int(delta_factor),
     *       (class_array(ii,j),ii=1,classcount)
          write(99,10004)int(delta_factor),
     *       (class_array(ii,j+1),ii=1,classcount)
        else
          write(99,*)classname(j)
          write(99,10005)(class_array(ii,j),ii=1,classcount)
          write(99,10006)(class_array(ii,j+1),ii=1,classcount)
        endif
      end do

      close(unit=99,status='keep')

	print*
	print*,'Please see `sensitivities.txt` in working directory' 
	print*,'for a summary of the sensitivities'
	print*

10000 format(a6,<nrvr>a10)
10001 format('-',i3,'%',<nrvr>f10.3)
10002 format('+',i3,'%',<nrvr>f10.3)
10010 format(a6,<classcount>a10)
10003 format('-',i3,'%',<classcount>f10.3)
10004 format('+',i3,'%',<classcount>f10.3)
10005 format('- 1dC',<classcount>f10.3)
10006 format('+ 1dC',<classcount>f10.3)


      return

      end subroutine write_sens_result