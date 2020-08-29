      SUBROUTINE ice_factor

!***********************************************************************
!    Copyright (C) 2016 by Nicholas Kouwen  
        
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

      integer     :: i,j,l,n,jz
      logical     :: firstpass

      data firstpass/.true./

!     ice_factor is called if icerivflg and/or icelakeflg = y

!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   
!     the equations were taken from the Lake Winnipeg model
      if(firstpass.eq.'y')then
!     rev. 10.3.03 Jan.  21/20  = NK added dd_ice & dd_thaw to the resume.txt file      
          if(resumflg.ne.'y')then
!             initialize the  ice factor
              do n=1,naa
                  dd_ice(n)=1000
                  dd_thaw(n)=0.0
              end do
          endif
!     REV. 10.1.19 Jan.  15/16  - NK: Fixed initialization of ice_factr - moved from lake_ice > runof6
c!       default ice_fctr        
c        do n=1,naa
c          ice_fctr(n)=1.0
c        end do
        open(unit=549,file='ice_factor.txt',status='unknown')
        firstpass=.false.
      endif
      
!       Calculate the ice factors for eacg grid
!       needed only once a day      
        if(jul_day_now.gt.275.or.jul_day_now.lt.91)then     
!         fall freeze up & winter           
          do n=1,naa
              if(dd_ice(n).gt.0.0)then
!               dd_ice (freezing degree days) must be +ve   
!               equation for Waterhen:       
!               lif=3.0981177*dd_ice(n)**-0.3324723
                ice_fctr(n)=3.7222496*dd_ice(n)**-0.4
!               in the winter is should not be less than 0.25
                ice_fctr(n)=amax1(0.50,ice_fctr(n))
!               and should never be larger than 1.0
                ice_fctr(n)=amin1(1.0,ice_fctr(n))
              endif
          end do
        else  !spring thaw & summer
          do n=1,naa
              if(dd_thaw(n).gt.0.0)then
!               taken from Lake Winnipeg              
c               ice_fctr(n)=1.0/(3.7222496*dd_thaw(n)**-0.2747330)
                ice_fctr(n)=1.0/(3.7222496*dd_thaw(n)**-.3)
!               in the winter is should not be less than 0.25
                ice_fctr(n)=amax1(0.50,ice_fctr(n))
!               in the summer it should not be more than 1.0
                ice_fctr(n)=amin1(1.0,ice_fctr(n))
c              else
c                ice_fctr(n)=1.0
              endif
c            endif
          end do  ! NOTE: ice_fctr written in rte.txt fln(55)
       endif
      
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
      if(.not.icefactorfile)then
        if(icelakeflg.eq.'y')then
!         ice factors for lakes if the resrl\ice_factor.tb0 
!         file does NOT exists 
          do l=1,noresv
            i=ires(l)
            j=jres(l)
            n=s(i,j)
!           lake ice factor should probably not be less than 0.50
            lake_ice_factor(l,month_now)=amax1(0.90,ice_fctr(n))
          end do
c          print*,'lake ice factors are computed'
        endif
      endif

      
      end subroutine ice_factor
