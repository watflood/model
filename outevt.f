      SUBROUTINE outevt(conv,scale,smc5,nhr,nhf)

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
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(12) :: date
      DIMENSION     :: smc5(16)
      integer  :: ios,nhr,nhf
      real*4   :: conv,scale,smc5

! WRITE A NEW EVENT.EVT FILE:
! SEE INPEVTA FOR VARIABLE DEFS

      open(unit=99,file=fln(99),status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,' problem opening ',fln(99),' in outevt.for @ line 15'
        print*,' ios=',ios,' possible problem: lile is locked'
        stop 'program aborted in outevt'
      endif

      write(99,5011,iostat=ios)date,snwflg

      if(ios.ne.0)then
        print*,' problem writing to',fln(99)
        print*,' ios=',ios
!        pause ' continue at your own risk!'
      endif

      write(99,5012)mo1
      write(99,5013)conv,scale
      write(99,5014)smc5  
      write(99,5015,iostat=ios)nhr
      print*,' ios=',ios
      write(99,5016)nhf
      write(99,6001)fln(1)
      write(99,6002)fln(2)
      write(99,6003)fln(3)
      write(99,6004)fln(4)
      write(99,6005)fln(5)
      write(99,6006)fln(6)
      write(99,6007)fln(7)
      write(99,6008)fln(8)
      write(99,6009)fln(9)
      write(99,6010)fln(10)
      write(99,6011)fln(11)
      write(99,6012)fln(12)
      write(99,6013)fln(13)
      write(99,6014)fln(14)
      write(99,6015)fln(15)
      write(99,6016)fln(16)
      write(99,6017)fln(17)
      write(99,6018)fln(18)
      write(99,6019)fln(19)

      close(unit=99,status='keep')

! THE DEFAULT IS TO RUN THE TEMPERATURE FILES:

! FORMATS

 5011 format('date & snow flag:         ',a14,a1,'nynnnnnnn')
 5012 format('month:                    ',i5)
 5013 format('rain conv. factor @ scale ',2f5.2)
 5014 format('ini moist each subbasin   ',5f5.2)
 5015 format('no hours of rain data     ',i5)
 5016 format('no hours of flow data     ',i5)
 6001 format('basin file name           ',a30,' unit= 31')
 6002 format('parameter file name       ',a30,'       32')
 6003 format('rain gage locn file name  ',a30,'       33')
 6004 format('stream gage locn file nam ',a30,'       34')
 6005 format('rain gage data file name  ',a30,'       35')
 6006 format('streamflow data file      ',a30,'       36')
 6007 format('reservoir release file    ',a30,'       37')
 6008 format('snow data file            ',a30,'       38')
 6009 format('radar met file            ',a30,'       39')
 6010 format('simple input met file     ',a30,'       40')
 6011 format('aes hourly radar rainfall ',a30,'       41')
 6012 format('clutter file              ',a30,'       42')
 6013 format('snow cover depln. curve   ',a30,'       43')
 6014 format('point temperatures        ',a30,'       44')
 6015 format('gridded temperatures      ',a30,'       45')
 6016 format('gridded min. temp.        ',a30,'       46')
 6017 format('gridded max. temp.        ',a30,'       47')
 6018 format('daily snow data file      ',a30,'       48')
 6019 format('spare                     ',a30,'       49')

      RETURN

      END SUBROUTINE outevt
