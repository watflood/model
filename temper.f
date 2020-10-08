      SUBROUTINE temper(n,ii,time,ttime)

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
     
!***********************************************************************

! THIS SUBROUTINE CALCULATES THE TEMPERATURE FOR THE CURRENT TIME STEP
!   -use general cos function, with no offset to represent
!    temperature variation between tmax and tmin
!   -tmax at hour 12
!   -tmin at hour 0

!   calling program = pmelt

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90

!***********************************************************************

      use area_watflood
	implicit none

      real*4      :: time,ttime,difft,sumt
      integer*4   :: hour,itime,n,ii

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

!     CALCULATE HOUR OF DAY      
      itime=int(time)
      hour=mod(itime,24)

!     1ST 12 HOURS OF DAY 
      if((hour.ge.0).and.(hour.le.12))then
         if(tmin1(n).eq.tmax(n))then
            ttime=tmin1(n)
         else   
            difft=tmin1(n)-tmax(n)
            sumt=tmin1(n)+tmax(n)
            ttime=-abs(difft)/2*cos(2*3.14159*hour/24)+sumt/2
         endif  
      else
!        2ND 12 HOURS OF DAY 
         if(tmin2(n).eq.tmax(n))then
            ttime=tmin2(n)
         else  
            difft=tmin2(n)-tmax(n)
            sumt=tmin2(n)+tmax(n)
            ttime=-abs(difft)/2*cos(2*3.14159*hour/24)+sumt/2
         endif  
      endif
      
!     DEBUGGING      
      if((n.eq.s(ipr,jpr)).and.(ii.eq.classcount).and.(numa.eq.0))then
         write(57,5000) difft,sumt
         write(57,5001) time,hour,n,tmin1(n),tmin2(n),tmax(n),ttime 
      endif


! FORMATS:

 5000 format(' diff,sum',2f6.1)	
 5001 format(' time,hour,n,tmin1,tmin2,tmax,ttime',
     *            2f6.1,i3,4f6.1)

      RETURN

      END SUBROUTINE temper

