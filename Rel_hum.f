      !-----------------------------------------------------------------
       SUBROUTINE REL_HUM(x1,rh,temp,dewpoint)
 
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
     
      !-----------------------------------------------------------------
      ! This subprogram is designed to accept temperature and dewpoint
      ! data and convert them to relative humidity data, using Newton's
      ! iterative Method.
      ! 
      ! rh :relative humidity
      ! new_rh_max :maximum relative humidity based on dewpoint
      !  temperature
	! dewpoint :dewpoint temperature
	! temp :temperature
	! x1 :USED IN OLD PROGRAM; TO BE REMOVED; PASSED BY REFERENCE AND 
      ! MUST ALSO BE REMOVED ONCE IT IS DETERMINED WHERE IT IS PASSED.
      !-----------------------------------------------------------------

      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

!     see Ken's thesis for e equation 4-12

      REAL ::rh,dewpoint,temp,new_rh_max,x1

      ! If the dewpoint is less than 0, then the maximum relative 
      ! humidity has to be adjusted to compensate for lower maximum  
      ! relative humidities at lower temperatures.  Otherwise, set to 
      ! 100%

      IF (dewpoint.LT.0)THEN
         new_rh_max=0.5*dewpoint+100
      ELSE
         new_rh_max=100
      END IF
       
!	rh=((exp(50*(1.31*dewpoint-36.2)/(4.21*dewpoint+1000)))/(/(6.11*
!     &exp((17.3*temp)/(temp+237.3))


	rh=((exp(50*(1.31*dewpoint-36.2)/(4.21*dewpoint+1000)))/(6.11*
     &       exp((17.3*temp)/(temp+237.3))))

      ! If the calculated rh is greater than the maximum rh at that 
      ! temperature, then it is reset to the maximum rh

      IF(rh.GT.new_rh_max)THEN
         rh=new_rh_max
      ! If it is less than 0%, then there is an error and 999.9999 is 
      ! outputted
      ELSE IF(rh.LT.0)then
         rh=999.9999
      END IF

      RETURN
      END SUBROUTINE REL_HUM
      !-----------------------------------------------------------------
