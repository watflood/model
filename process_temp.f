      SUBROUTINE process_temp(jz)

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
     
C***********************************************************************
C  PROCESS_TEMP - written Jun/06 by Dave Watson
C	- Derived from process routines written by Nick Kouwen inside rdtemp 
C	- This subroutine processes the temperatures loaded by read_temp_ef()
C***********************************************************************

!     rev. 9.5.63  Sep.  04/09  - NK: moved lapse rate from melt.f to process_temp.f
!     rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m

      USE area_watflood
	implicit none
	save

	integer i,j,n,i2,j2,m
	integer jz
	integer ntocorrect


      logical      :: offset,not_offset,remap,firstpass
      data firstpass/.true./

	ntocorrect = 0 !default value
      
!     remap the precip to a finer grid
      if(firstpass)then
        if(iyoffset.eq.0.and.jxoffset.eq.0)then
          not_offset=.true.
        else
          offset=.true.
        endif
        if(xcount2.ne.xcount.or.ycount2.ne.ycount)then
          remap=.true.
          m=int((xdelta2+0.000001)/xdelta)
          n=int((ydelta2+0.000001)/ydelta)
          print*,'temp grid:   xratio=',m,' yratio=',n
        endif
      endif

d     if(jz.le.1)then
d       if(iopt.ge.2)then
d	    print*,'iyoffset=   ',iyoffset
d	    print*,'jxoffset=   ',jxoffset
d         print*,'not_offset= ',not_offset
d         print*,'offset=     ',offset
d       endif
d     endif

      if(not_offset.and.remap)then
        do i=ycount,1,-1
          do j=xcount,1,-1
            ttemp(i,j)=ttemp((i+1)/n,(j+1)/m)
          end do
        end do
      endif

      if(offset)then
!       REMAP INTO SOUTH-WEST CORNER
        do i=1,ycount
           do j=1,xcount
             ttemp(i,j)=ttemp(i+iyoffset,j+jxoffset)
           end do
        end do
      endif
!     rev. 10.2.48 Feb.  25/19  - NK: Moved temp correction to read_temp from process_temp
!      moved to read_temp because it should be done on the read values
!      not over & over between reads
c!     RAISE OR LOWER THE TEMPERATURE FILED:
c      if(scaletem.gt.0.00001.or.scaletem.lt.-0.00001)then
c         do i=1,ycount
c            do j=1,xcount
c               ttemp(i,j)=ttemp(i,j)+scaletem
c            end do
c         end do
c      endif
      
!     rev. 9.01    Aug.   1/00  - added look up for minimum temperature
!                                 and function to calculate RH
!     reset the minimum temperature each midnight:
!     the offset for finding the min is a8 hours before midnight
!     a8 is entered in the parameter file.
!
!     only the centre grid is used ??????? fix fix ????
      if(mod(jz+int(a8),24).eq.0)then
        if(ttemp(yyy(naa/2),xxx(naa/2)).gt.-90.0)then
!         reset only if there is temperature data          
          do n=1,na
            tempvmin(n)=99.9
          end do
        endif
      endif

!     PUT INTO VECTOR FORMAT
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         if(ttemp(i,j).gt.-90.0)then
!          if missing data, leave temperature the same 
!          as last recorded value
!          store seed value for rh from previous calculation
!          initial value set in sub
!           x1=rh(n)
           tempv(n)=ttemp(i,j)
!          find minimum temperature for each day for rel humidity comp.
           tempvmin(n)=min(tempv(n),tempvmin(n))
         endif

!       rev. 9.5.63  Sep.  04/09  - NK: moved lapse rate from melt.f to process_temp.f
!       rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m

!        this is now done in tmp.exe
c        if(tlapse.ne.0.0)then
c!         ELEVATION IS IN HUNDREDS OF METERS
c!         LAPSE RATE IS IN DEGREE C / 1 m - USUALLY 0.5 DEGREE C
c          tempv(n)=tempv(n)-(elev(n)-elvref)*tlapse
c        end if
      end do        

!     rev. 9.7.11  Nov.  21/10  - NK: added monthly_climate_deltas.txt file
	if(climate_delta_flg)then
        do n=1,naa
          tempv(n)=tempv(n)+monthly_temperature_delta(month_now)
	  end do
	endif

      firstpass=.false.

      RETURN 

	END SUBROUTINE process_temp