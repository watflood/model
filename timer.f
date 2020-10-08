      SUBROUTINE timer(iz,jz,mz,clock,time,t,thr,dtmin,
     *                dtmax,div,m,ju,a66)

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
     
!     rev. 9.2.43  Jun.  21/06  - NK: fixed spikes in route
!     rev. 9.4.05  May.  04/07  - NK: revised timer for julian day calc.
!     rev. 9.5.21  Mar.  06/08  - NK: fixed dtmin for first time step each event
!     rev. 9.5.70  Oct.  11/09  - NK: fixed timer for r2c frames (use year_now)
!     rev. 10.1.66 Jan.  28/17  - NK: Fixed leap year in timer 
!                                 Also rewrote second portion of this s/r
!     rev. 10.1.88 May   23/17  - NK: Fixed Julian_day problems for iso R/W

!***********************************************************************

!  DEBUG INFORMATION IS WRITTEN TO UNIT 55 - RTE.LST

!     REV. 9.00    Mar.  2000 - TS: CONVERTED TO FORTRAN 90 

!  t			- time increment in seconds
!  thr		- time increment in hours
!  time		- time from the beginning of the event in hours
!  totaltime	- time from beginning in hours
!  mz			- the integer value of the hour
!				smallest time interval allowed will be a66 seconds (param file)
!  dtmin		- found in 'route' and is equal to the smallest travel
!				through a square found in the computation during the time inc
!				being considred.

!***********************************************************************

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

        INTEGER        :: mohours_ly(24),ju_mon_ly(24)
        INTEGER        :: mohours(24),ju_mon(24),
     *                    m,iz,ik,jz,mz,ju,i
        real(4)        :: time,t1
        real*4         :: a66,t,div,thr,clock,dtmax,dtmin,tmv(3)
        logical        :: firstpass,next_leapyear

      DATA mohours/744,672,744,720,744,720,744,744,720,744,720,744,
     *             744,672,744,720,744,720,744,744,720,744,720,744/
      DATA ju_mon/1, 32, 60, 91,121,152,182,213,244,274,305,335,
     *            1, 32, 60, 91,121,152,182,213,244,274,305,335/

!     rev. 10.1.66 Jan.  28/17  - NK: Fixed leap year in timer 
      DATA mohours_ly/744,696,744,720,744,720,744,744,720,744,720,744,
     *             744,696,744,720,744,720,744,744,720,744,720,744/
      DATA ju_mon_ly/1, 32, 61, 92,122,153,183,214,245,275,306,336,
     *            1, 32, 61, 92,122,153,183,214,245,275,306,336/
c     *          366,397,425,456,486,517,547,578,609,639,670,700/
      data firstpass/.true./
      

!     DURING PRECIP, EACH EVEN HOUR IS INCREMENTED 
!     AFTER PRECIP, EACH KT HOUR WILL BE INCORPORATED
!     M IS EQUAL TO 1 IN THE 'RESET' SUBROUTINE AT THE STRAT    
!          OF EACH NEW HYDROGRAPH
!     THE TIME INCREMENT IS ALLOWED TO DECREASE RSMCDLY BUT TO INCRE.
!          Y SLOWLY - STAETEMNT 2 - TO PREVENT HYRDRAULIC INSTABILITY

!     a temporary fix for a66 getting clobbered
!      a66=900.0

      if(iopt.gt.1)then
         write(55,6006)
         write(55,6005)
     *         iz,jz,mz,itogo,t1,dtmin,dtmax,clock,time,t,thr
      endif

c      if(m.le.1)then
      if(firstpass)then
!        THIS SECTION IS ONLY USED AT THE START OF EACH EVENT

!     rev. 9.2.43  Jun.  21/06  - NK: fixed spikes in route
         t=3600.0
c	   time=1.0  ! done in sub now
         
!         time reset in sub =1.00 for new event
!         This is the main thing - to reset to 1.0 so the first time frame
!         is read in any r2c input files          
          jz=int(time)   
          mz=jz+1

          itogo=1
          div=t/2.0
          thr=t/3600.
          iz=-1
          iz=0
          clock=0.0
          dtmax=3600.0
          do ik=1,3
              tmv(ik)=a66
          end do
c!         time reset in sub =1.00 for new event
c!         This is the main thing - to reset to 1.0 so the first time frame
c!         is read in any r2c input files          
c          jz=int(time)   
c          mz=jz+1
         
!     rev. 10.1.93 Aug   17/17  - NK: allow year1 etc. to be passed for each event
!         no change here !!!!!!!
          if(firstpass)then
              year_now=year1
              month_now=mo1
              if(mod(year_now,4).eq.0)then
                  leapyear=.true.
              else
                  leapyear=.false.
              endif
              if(leapyear)then
                  ju=ju_mon_ly(mo1)+day1-1    ! changed May4/07  nk
              else
                  ju=ju_mon(mo1)+day1-1    ! changed May4/07  nk
              endif
              ju=max(ju,1)
              jul_day_now=ju           ! added Apr. 06/15 NK

!     rev. 10.1.66 Jan.  28/17  - NK: Fixed leap year in timer 
!     rev. 10.1.67 Feb.  18/17  - NK: Ignore start year in subsequent event files 
              if(leapyear)then
                  day_now=ju-ju_mon_ly(mo1)+1
              else
                  day_now=ju-ju_mon(mo1)+1
              endif
              hour_now=mod(jz,24)
              hour_now_fews=hour_now
              if(hour_now.eq.24)then
                  hour_now_fews=0
              else
                  hour_now_fews=hour_now
              endif

              if(iopt.ge.3)then
                  write(55,6003)
                  write(55,6005)
     *            iz,jz,mz,itogo,t1,dtmin,dtmax,clock,time,t,thr
              endif
          endif
         
c         firstpass=.false.

c         RETURN      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----------------

      else
         
         iz=jz
!        WHEN IZ.NE.JZ A NEW RAIN WILL BE READ IN
!        (this is probably no longer needed because this s/r is alled only hourly)         
         jz=int(time)  ! time is incremented in sub
         mz=jz+1

!        MIN TRAVEL TIME IS DETERMINED IN ROUTE - SHOULD NOT BE .LT. a66
!        SO MAX TIME STEP IN SECONDS:
         dtmax=amax1(dtmin,a66)
!        ROUND THIS OFF SO THAT IT CAN BE DIVIDED INTO 3600:       
c         dtmax=float(int(dtmax)/int(a66))*a66
         dtmax=float(int(dtmax/a66))*a66

!        CURRENT TIME IS TIME & JZ IN HRS OR T IN SECONDS
!        CALCULATE THE TIME UNTIL THE NEXT RAINFALL
!        THE EXACT TIME UNTIL THE NEXT RAIN IS: 

         t1=float(jz+itogo-1)-time
         t1=max(1.0,t1)
!        THE NEXT TIME STEP IS MAX OF DTMAX AND TIME TO NEXT RAINFALL:
!        FOR SNOW MELT, THE MAX TIME STEP IS 1 HOUR
         t=min(t1*3600.0,dtmax,kt*3600.0)
         div=t/2.0
         thr=t/3600.
c         iz=jz
c!        WHEN IZ.NE.JZ A NEW RAIN WILL BE READ IN
c!        (this is probably no longer needed because this s/r is alled only hourly)         
c         jz=int(time)  ! time is incremented in sub
c         mz=jz+1

!     rev. 9.4.05  May.  04/07  - NK: revised timer for julian day calc.
!     rev. 10.1.74 Apr.  01/17  - NK: Changed timer to fix 1 day-off problem 
c         if(mod(jz-1,24).eq.0)then
        if(mod(jz,24).eq.0)then
	     ju=ju+1
c          write(778,7777)' ',ju   !,jz,hour_now,day_now,month_now,year_now,leapyear

           !          reset the julian day on Jan 1 each year
!          First see if next year is a leap year

!         see if we're in a leapyear:
           if(mod(year_now,4).eq.0)then
              leapyear=.true.
           else
              leapyear=.false.
           endif           
           
           if(leapyear)then
	        if(ju.gt.366)then
                ju=1
                year_now=year_now+1
              endif
	     else
	        if(ju.gt.365)then
                ju=1
                year_now=year_now+1
               endif
           endif
c           ju=ju+1
         endif

         jul_day_now=ju

!     rev. 9.5.70  Oct.  11/09  - NK: fixed timer for r2c frames (use year_now)
!        FIND OUT WHAT MONTH WE ARE IN:
         do i=12,1,-1
           mo=i  
           if(leapyear)then
             if(ju.ge.ju_mon_ly(mo))then
               month_now=mo
               go to 1000
             endif
           else
             if(ju.ge.ju_mon(mo))then
               month_now=mo
               go to 1000
             endif
           endif
         end do
!        Do not change this go to - 
!        need to jump out of the loop before it's done
 1000    if(mo.gt.12)month_now=month_now-12

         If(leapyear)then
           day_now=ju-ju_mon_ly(month_now)+1
         else
           day_now=ju-ju_mon(month_now)+1
         endif

        hour_now=mod(jz,24)
        if(hour_now.eq.24)then
            hour_now_fews=0
        else
            hour_now_fews=hour_now
        endif
        
!       take this out and it won't read next data
        if(hour_now.eq.0)then
          hour_now=24
        endif
        
c        if(mod(hour_now,1).eq.0)then
c          write(778,7777)
c     *     'b',ju,jz,hour_now,day_now,month_now,year_now,leapyear
c7777      format(a1,6I5,5x,l1)
c        endif
      
      endif

      if(iopt.ge.3)then
!         if(mz-mz/12*12.eq.0)write(55,6003)
         if(jz-mz/12*12.eq.0)write(55,6003)
         write(55,6005)
     *         iz,jz,mz,itogo,t1,dtmin,dtmax,clock,time,t,thr
      endif

!     rev. 10.1.64 Jan.  26/17  - NK: Added XML output file 
c      if(mod(hour_now,24).eq.0)then
      if(hour_now.eq.1)then
        if(month_now.le.9.and.day_now.le.9)then
          write(yyyymmdd12(ju),6011)year_now,month_now,day_now
6011    format(i4,'-0',i1,'-0',i1)      
        elseif(month_now.gt.9.and.day_now.le.9)then
          write(yyyymmdd12(ju),6012)year_now,month_now,day_now
6012    format(i4,'-',i2,'-0',i1)      
        elseif(month_now.le.9.and.day_now.gt.9)then
          write(yyyymmdd12(ju),6013)year_now,month_now,day_now
6013    format(i4,'-0',i1,'-',i2)      
        elseif(month_now.gt.9.and.day_now.gt.9)then
          write(yyyymmdd12(ju),6014)year_now,month_now,day_now
6014      format(i4,'-',i2,'-',i2)      
        endif
c        print*,yyyymmdd12(ju)(1:12),' in timer'
      endif
! FORMATS
      
!     rev. 10.1.90 Jul.  27/17  - NK: Added date_now for i/o files
      if(month_now.le.9.and.day_now.le.9)then
          write(date_now,20001)year_now,month_now,day_now
20001     format(i4,'-0',i1,'-0',i1)
      elseif(month_now.le.9.and.day_now.gt.9)then
          write(date_now,20002)year_now,month_now,day_now
20002     format(i4,'-0',i1,'-',i2)
      elseif(month_now.gt.9.and.day_now.le.9)then
          write(date_now,20003)year_now,month_now,day_now
20003     format(i4,'-',i2,'-0',i1)
      elseif(month_now.gt.9.and.day_now.gt.9)then
          write(date_now,20004)year_now,month_now,day_now
20004     format(i4,'-',i2,'-',i2)
      endif

      if(firstpass)then
          startdateXML=date_now
          starttimeXML='00:00:00'
      endif

 6003 format('   iz   jz   mz  itogo   t1    dtmin   dtmax   clock
     *    time    t       thr')
 6005 format(4i5,f8.2,2f8.0,2f8.2,f8.0,f8.2)
 6006 format(' into timer')

      firstpass=.false.
      
      RETURN
      END SUBROUTINE timer

