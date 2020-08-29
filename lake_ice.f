      SUBROUTINE lake_ice

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
     
!     REV. 10.1.16 Jan.  11/16  - NK: Added subroutine ice_factor.f

      use area_watflood
      implicit none

      integer     ::  n,jz
      logical     :: firstpass

      data firstpass/.true./


!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   
!     the equations were taken from the Lake Winnipeg model
      if(firstpass.eq.'y')then
!         initialize the  ice factor
        do n=1,naa
          dd_ice(n)=1000
          dd_thaw(n)=0.0
        end do
!     REV. 10.1.19 Jan.  15/16  - NK: Fixed initialization of ice_factr - moved from lake_ice > runof6
c!       default ice_fctr        
c        do n=1,naa
c          ice_fctr(n)=1.0
c        end do
        open(unit=549,file='ice_factor.txt',status='unknown')
        firstpass=.false.
      endif
      
c      if(iceflg.eq.'y'.and.mod(jz,24).eq.0)then
!       needed only once a day      
        if(jul_day_now.gt.275.or.jul_day_now.lt.91)then     
!         fall freeze up & winter           
          do n=1,naa
c            if(IBN(n).eq.5)then
              if(dd_ice(n).gt.0.0)then
!               dd_ice (freezing degree days) must be +ve          
c               lif=3.0981177*dd_ice(n)**-0.3324723
c                ice_fctr(n)=3.7222496*dd_ice(n)**-0.2747330
                ice_fctr(n)=3.7222496*dd_ice(n)**-0.4
c                ice_fctr(n)=amax1(1.0,ice_fctr(n))
!               in the winter is should not be less than 0.5
                ice_fctr(n)=amax1(0.25,ice_fctr(n))
                ice_fctr(n)=amin1(1.0,ice_fctr(n))
              endif
c              if(ice_fctr(n).lt.0.95)ice_fctr(n)=0.5
c            endif
          end do
        else  !spring thaw & summer
          do n=1,naa
c            if(IBN(n).eq.5)then
              if(dd_thaw(n).gt.0.0)then
c                ice_fctr(n)=1.0/(3.7222496*dd_thaw(n)**-0.2747330)
                ice_fctr(n)=1.0/(3.7222496*dd_thaw(n)**-0.2747330)
                ice_fctr(n)=1.0/(3.7222496*dd_thaw(n)**-.3)
c                ice_fctr(n)=amax1(1.0,ice_fctr(n))
!               in the summer it should not be more than 1.0
                ice_fctr(n)=amax1(0.25,ice_fctr(n))
                ice_fctr(n)=amin1(1.0,ice_fctr(n))
c              else
c                ice_fctr(n)=1.0
              endif
c            endif
          end do  ! NOTE: ice_fctr written in rte.txt fln(55)
        endif
c      endif  ! iceflg.eq.y
      
      end subroutine lake_ice
