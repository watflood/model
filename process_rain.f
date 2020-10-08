      SUBROUTINE process_rain(conv,scale)

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
C  PROCESS_RAIN - written Jun/06 by Dave Watson
C     - Derived from process routines written by Nick Kouwen inside rdrain 
C     - This subroutine processes the rain loaded by read_rain_ef()
C***********************************************************************

!     rev. 9.3.06  Dec.  17/06  - NK: added precip adjustment for bias
!     rev. 9.3.07  Dec.  29/06  - NK: added sum_precip for whole domain
!     rev. 9.5.17  Feb.  28/08  - NK: moved scale snow from sub to process rain
!     rev. 9.5.30  May.  26/08  - NK: conv back in read_rain & process_rain arg. list
!     rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg

      use area_watflood
      implicit none
      save

      integer i,j,n,ios,i2,j2,m
      real*4   :: conv,scale,conv1
      logical  :: exists,offset,not_offset,remap,firstpass
c      real*4, dimension(:,:),   allocatable ::  array
      
      data firstpass/.true./
      
!     remap the precip to a finer grid
      if(.not.netCDFflg.and.firstpass)then
        if(iyoffset.eq.0.and.jxoffset.eq.0)then
          not_offset=.true.
        else
          offset=.true.
        endif
        if(xcount2.ne.xcount.or.ycount2.ne.ycount)then
          remap=.true.
          m=int((xdelta2+0.000001)/xdelta)
          n=int((ydelta2+0.000001)/ydelta)
          print*,'precip grid: xratio=',m,' yratio=',n
        endif
c		allocate(array(ycount,xcount),stat=iAllocate)
      endif

      if(not_offset.and.remap)then
        do i=ycount,1,-1
          do j=xcount,1,-1
!           avoid div by 0 for i=1 & j=1          
            p(i,j)=p((i+n-1)/n,(j+m-1)/m)
          end do
        end do
      endif
            
!     rev. 9.3.06  Dec.  17/06  - NK: added precip adjustment for bias
!     adjust bias found in radar precip using long-term accumulation
!     this takes precedence over any values from the paf files read by
!     precip_adjust.f - see below at *****
      if(firstpass)then
        inquire(FILE='under_by.csv',EXIST=exists)
        if(exists)then
          print*,'file under_by.csv  found'
          open(unit=99,file='under_by.csv',status='old',iostat=ios)
          if(ios.ne.0)then
            print*,'Problem opening file  under_by.csv'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*
            stop 'Program aborted in process_rain @ 32'
          endif
          precflg=.true.
          do i=ycount,1,-1
            read(99,*)(precadj(i,j),j=1,xcount)
            write(*,*)i,(precadj(i,j),j=1,xcount)
          end do
          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            if(precadj(i,j).ne.0.0)then
              precadj(i,j)=1.0/precadj(i,j)
            else
              precadj(i,j)=1.0
            endif
          end do
        else
        
!         This needs to be done only at the start of a run. 
!         Process_rain is done every time step.
          if(pafflg.eq.'y')then
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call precip_adjust()    ! *****
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif
        endif
      endif

!       REV. 8.98   July  15/99 -  MET GRID SHIFTING FOR WEATHER MODELS  
!       REMAP INTO SOUTH-WEST CORNER:
!     rev. 9.3.07  Dec.  29/06  - NK: added sum_precip for whole domain

      if(offset)then
        do i=1,ycount
          do j=1,xcount
            p(i,j)=p(i+iyoffset,j+jxoffset)
          end do
        end do
      endif

C NICK from DAVE - do we need the following? or would precip already be adjusted
!       THIS SECTION MULTIPLIES THE PRECIP BY THE PRECIP ADJUSTMENT
!       FACTOR
!       FRANK S JAN 98 (IF PREADJ=0.0 THEN THERE IS NO CHANGE)
      if(precflg)then
            do n=1,naa
                  i=yyy(n)
                  j=xxx(n)
                  p(i,j)=p(i,j)*precadj(i,j)
            end do
            if(iopt.eq.2)print*,'grid adjusted'
      endif

C NICK from DAVE-I would like to get away from this guage/radar differentiation
C we can set the entire MET file to gauge or radar in the header..but should we assume the data is 
C already converted

      if(conv.ne.1.0.or.scale.ne.1.0)then

            if(source.ne.'Gauge')then
!           RADAR IS SCALED IF THE OPTION IS USED BUT ONLY THE
!           RADAR DATA                
                  conv1=conv*scale
            else
!           WHEN MISSING RADAR IS REPLACED BY GAUGE DATA, 
!           IT IS NOT SCALED
                  conv1=conv
            endif
            do n=1,naa
              i=yyy(n)
              j=xxx(n)
              p(i,j)=p(i,j)*conv1
            end do
      endif         
      
!     rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg
!     DN to NICK -- this is the only change required here!
!     (all the meat is in the sub)
!     rev. 2014-08 - DN: hook for precip corrections
c      if(fcstflg.eq.'y')then
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call forecast_correct_precip()    ! correction happens in the sub
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       2014-08 -- DN
!       modifies p(i,j) in forecast mode (if applicable) based on adjustment factors
!       specified in forecast parameter file
        !print *, 'in process_rain@151, checking forecast mode flag'
        !print *, 'in process_rain@152, checking forecast mode flag'
        if(fcst_mode)then
          !print *, 'in process_rain@153, correcting precip!'
          !print *,'fcst_hr0,hour_now,hour1,hour_start:'
          !print *, fcst_hr0, hour_now, hour1, hour_start
          !print *, fcst_hr0, hour_now+24*jul_day_now, fcst_snow_adj
          if(fcst_hr0.gt.(hour_now+24*jul_day_now))then
          !         DN 2015-02-07 -- these are temporary debug statements
            !print*,'correcting precip at hour:',hour_now+24*jul_day_now
            !print*,hour_now+24*jul_day_now,'year:',year_now
            !print *, 'tempv: ',tempv(1),'rainsnowtemp: ',rainsnowtemp
            !print *, 'adj: ',fcst_snow_adj,fcst_rain_adj
            !print *, 'p_ij: ',p(1,1),'adjusted: ',p(1,1)*fcst_snow_adj
            !print *, p(1,1)*fcst_rain_adj
            do n=1,naa
!             which variable to get <current_hour> from reliably?
              i=yyy(n)
              j=xxx(n)
            ! use global rainsnowtemp from evt file to separate snow from rain
            ! often set to 0
              if(tempv(n).lt.rainsnowtemp)then    
                p(i,j)=p(i,j)*fcst_snow_adj
                !print *, 'adjusting snow! at ',hour_now+24*jul_day_now
              else
                p(i,j)=p(i,j)*fcst_rain_adj
                !print *, 'adjusting rain! at ',hour_now+24*jul_day_now
              end if
            end do
          end if
        end if
c      endif   ! fcstflg
!     end rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg

!     rev. 9.2.36  Mar.  30/06  - NK: Scaleallsnow changed to precip snow
!     rev. 9.5.17  Feb.  28/08  - NK: moved scale snow from sub to process rain
!     this little ditty was just a temp fix to scale precip when it is snow.
      if(scalesnw.gt.0.and.snwflg.eq.'y')then
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          if(tempv(n).lt.0.0)then    ! could pick a larger temp.
            p(i,j)=p(i,j)*scalesnw
          endif
        end do
      endif
      
!     rev. 10.4.22 Apr.  22/20  = NK fixed event precip scale factor  "readscale"
c      if(.not.scale.eq.1.0.and.scale.gt.0.01)then
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
            p(i,j)=p(i,j)*scale
        end do
c      endif

!     rev. 9.7.11  Nov.  21/10  - NK: added monthly_climate_deltas.txt file
	if(climate_delta_flg)then
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
            p(i,j)=p(i,j)*monthly_precipitation_delta(month_now)
	  end do
c	  print*,'month=',month_now,' delta_T=',
c               	  monthly_temperature_delta(month_now)
	endif

c        write(955,*)'in process rain time=',totaltime
c        do i=ycount,1,-1
c          write(955,99955)(p(i,j),j=1,xcount)
c99955     format(999f5.1)          
c        end do
c        write(955,*)totaltime

99999 firstpass=.false.

      return

      END SUBROUTINE process_rain