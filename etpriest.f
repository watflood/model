      SUBROUTINE etpriest(mon)

      
!***********************************************************************
!    Copyright (C) 1996 by Todd Neff 
        
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
! Program by: Todd A.M. Neff,  August 1996

! THIS SUBROUTINE CALCULATES PRIESTLY-TAYLOR PET FOR EACH ELEMENT &
! CALCULATES EVAPOTRANSPIRATION FROM PRIESTLY-TAYLOR EQUATION  
! ASCE MANUAL NO. 70 - REV. APRIL 1982 - PAGE 145

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90

!  alamb      - lambda (eqn 7.1) - latent heat of vapourization
!  delt       - delta (slope of saturation vap press curve - eqn 7.12)
!  gam        - psychrometric constant (eqn 7.15)
!  petn       - Priestley-Taylor PET 
!  press(mon) - mean monthly pressure (kPa)
!  radv       - vectored mean hourly net radiation (W/sq.m.)

! NOTE: petn is adjusted based on the time increment
!  i.e. for 1 hour it is multiplied by 3.6 since there are 3600 seconds
!       per hour and is divided by 1000 to conserve units

!***********************************************************************

      USE area_watflood
      implicit none

      real*4   :: gam,delt,ddg,petn
      integer*4 :: n,ii,mon

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

! adjust alpha throughout the year, just for mackrpn
! Jan, Feb, March =0.0, April=0.6, May=1.2,
! June,July,August=1.35, Sept=1.2, Oct=0.6, Nov,Dec=0.0

      alpha=0.942857+0.00285714*real(mon)**4-0.08*real(mon)**3+
     +0.682857*real(mon)**2-1.72*real(mon)
      alpha=max(0.0,alpha)
      alpha=min(1.35,alpha)
      gam=.001013/0.622*pres(mon)/alamb

      do n=1,naa

         if(tempv(n).ge.-5)then
            delt=25083./(tempv(n)+237.3)**2*exp
     *          (17.3*tempv(n)/(tempv(n)+237.3))
            ddg=delt/(delt+gam)
            petn=alpha*ddg*radv(n)/alamb/den*3.6
         else
            petn=0
         endif

c      if(n.eq.1152) print *, ddg,den

         if(petn.lt.0) petn=0

         do ii=1,classcount
            pet(n,ii)=petn
         end do

      end do

      RETURN

      END SUBROUTINE etpriest
