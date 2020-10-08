      SUBROUTINE process_hum(time,n)

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
! PROCESS_HUM written DEC/076 by Trish Stadnyk
!	- converts specific humidity from REMOiso output into relative
!       humidity required for isoWATFLOOD 
!	- This subroutine processes the specific humidity read by read_hum()
!***********************************************************************

      use area_watflood
	use areacg
	implicit none

      INTEGER :: n
      REAL*4  :: time,ea,es,pa,qa

!     change spech(n) into rel(h) based on
!     http://fish.cims.nyu.edu/project_aomip/forcing_data/atmosphere/humidity.html

      qa=spech(n)   ! specific humidty in g/kg
	pa=101325     ! atmospheric air pressure (Pa)

!     vapour pressure in air:
      ea=pa/(0.622/qa+0.378)
     
!     es - saturation vapour pressure (as in over oceans or water bodies)
      es=10**((0.7859+0.03477*tempv(n))/(1.0+0.00412*tempv(n))+2)

!     relative humidty (in percent) = ratio of vapour pressure to sat'n vap. pressure:
      relh(n)=ea/es


      write(103,1003)time,n,qa,pa,ea,es,relh(n)
 1003 FORMAT(f8.2,',',i8,10(f14.4))

	END SUBROUTINE process_hum