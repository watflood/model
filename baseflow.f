      subroutine baseflow(n,dlz,sdlz,tdum)

!     s/r created May 5/03 NK
      
!***********************************************************************
!    Copyright (C) 2003 by Nicholas Kouwen  
        
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
      
!     rev. 9.4.08  May.  29/07  - NK: changed baseflow argument list

!     called in RUNOF6 and ROUTE

      USE area_watflood
	implicit none

      integer    :: n
      real*4     :: dlz,sdlz,tdum

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!       GROUNDWATER OUTFLOW:
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

!       DLZ = LOWER ZONE OUTFLOW IN MM

        dlz=flz2(n)*lzs(n)**pwr(n)
        if(dlz.gt.lzs(n))then
          dlz=lzs(n)
          lzs(n)=0.0
        else
          lzs(n)=lzs(n)-dlz
        endif

c	if(iopt.eq.2)print*,' checkpoint 5a in baseflow'

!       section to create runoff file for watroute with no 
!       groundwater component  nk  06/07/00
	  if(flz(1).ge.0.1E-12)then
          qlz(n)=dlz*tdum*frac(n)
        else
	    qlz(n)=0.0
	  endif

c	if(iopt.eq.2)print*,' checkpoint 5b in baseflow'

	  leakage=leakage+qlz(n)
!        qdrng=qdrng+rechrg(n)*tdum*frac(n)

        if(n.eq.nnprint)then
!       CALC LZ OUTFLOW
          sdlz=sdlz+dlz
        endif

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      return

      end subroutine baseflow