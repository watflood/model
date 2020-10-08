      subroutine write_wfo(t,jz,iz)

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen and Dave Watson (NRC)
        
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
     
!     rev. 9.1.44  Jun.  11/03  - Added Cumulative precip to the wfo file
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1

      use area_watflood
      use areacg
      
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

        REAL(4)		::     wfo_spec_version_number,tot1,t
        INTEGER       :: nratio,nj,jan,iallcnt1,i,j,ios,n,ii,jj,inta,
     *                  iallocate,l,npick,jz,iz
      INTEGER :: wfo_write_timestamp
      INTEGER :: wfo_write_attribute_data

	DATA iallcnt1/0/


!     NOTE: FOR MONTHLY RUNS, DIMENSIONS CAN BE CHANGED FROM 
!           3,8784,500  TO  12,366,3000

!>>>>>>>>>>>>>  AB: STUFF FOR ENSIM
      INTEGER(4) :: wfo_yy,wfo_mm,wfo_dd,wfo_hh,wfo_mi,wfo_ss,
     *                   wfo_ms
      INTEGER(4) :: wfo_seq


!       WRITE STUFF FOR ENSIM:
!       WRITE STUFF FOR ENSIM:
!       WRITE STUFF FOR ENSIM:
!       WRITE STUFF FOR ENSIM:
!       WRITE STUFF FOR ENSIM:
!       but not if we're optimizing 

!     rev. 9.1.51  Jan.  28/04  - NK: added iz.ne.jz conditional to ENSIM output  
!     Before this revision stuff was written to the WFO file at intervals < 1 hr
!     i.e. at the routing time steps and not hourly or greater.
!     Now stuff can not be reported at smaller than 1 hr dt's
!     rev. 9.2.20  Oct.  28/05  - NK: WFO_SPEC - reporting start & finish times 


!     Accumulate the precip during the current reporting time step
!           (not the modelling time step)
      if(wfo_pick(2).eq.1)then
        do n=1,naa
          ii=yyy(n)
          jj=xxx(n)
          wfo_sum_p(n)=wfo_sum_p(n)+p(ii,jj)
        end do
      endif

!     rev. 9.1.50  Jan.  14/04  - NK: version number added to the wfo_spec.txt file
!     Total precip accumulation since the run began
!     This is modified to be able to accumulate for a period
!     shorter than the model run so the total precip for one event
      if(wfo_pick(3).eq.1)then
!     can be shown in ENSIM
        do n=1,naa
          ii=yyy(n)
          jj=xxx(n)
          wfo_cum_p(n)=wfo_cum_p(n)+p(ii,jj)
        end do
      endif
      
!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
      if(wfo_pick(8).eq.1)then
!     can be shown in ENSIM
        do n=1,naa
          wfo_sum_qsyn(n)=wfo_sum_qsyn(n)+qo2(n)
        end do
      endif
!     write wfo file only at chosen interval

      if(mod(jz,min(ireport,mhtot)).eq.0.and.iz.ne.jz)then
!       need the min(ireport,mhtot) for the case where the user
!       specifies ireport to be greated than the event length
!       which happens often because months are different lengths

        wfo_seq=wfo_seq+1
!       SET TIME YR,MONTH,DAY,HR,MIN,SEC,MILLISEC FOR THIS STEP
c        wfo_yy=0
c        wfo_mm=0
c        wfo_dd=0
c        wfo_hh=int(totaltime)
c!       wfo_hh=int(time)
c        inta=wfo_hh
c        wfo_mi=int((totaltime-float(inta))*60.0)
c        wfo_ss=0
c        wfo_ms=0
        
!     REV. 10.1.33 Jun   20/16  - NK: Change the time stamp in the watflood.wfo file
!     year_now,month_now,day_now,hour_now
      wfo_yy=year_now
      wfo_mm=month_now 
      wfo_dd=day_now   
      wfo_hh=hour_now
      
!     rev. 10.1.68 Feb.  18/17  - NK: Made midnight 00 instead of 24 
      if(wfo_hh.eq.24)wfo_hh=0
      wfo_mi=0
      wfo_ss=0
      wfo_ms=0
      
c      WRITE(123,*)wfo_yy,wfo_mm,wfo_dd,wfo_hh       ! <<<<<<<<<<<<<<<<<<<<<<<<<<< take out
      
!       WRITE THE TIMESTAMP
        if(wfo_write_timestamp(wfo_seq,wfo_seq,wfo_yy,wfo_mm,
     *    wfo_dd,wfo_hh,wfo_mi,wfo_ss,wfo_ms).ne.1) then
          write (*,'(A)') ' '
          write (*,'(A)') '  *** FATAL ERROR ***     '
          write (*,'(A)') ' Unable to Write Timestamp'
          STOP 'Program terminated in sub @780'
        end if


!               ATTRIBUTE 1 - TEMP
!               I'M DOING THIS B/C TEMP IS NOT DYNAMICALLY ALLOCATED SO
!               IT IS TOO LARGE AND SO WE NEED TO USE OUTWFO
                if(wfo_pick(1).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=ttemp(ii,jj)
                  end do
                  if(wfo_write_attribute_data
     *                    (xcount,ycount).ne.1)then
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!               ATTRIBUTE 2 - PRECIP
!               I'M DOING THIS B/C P IS NOT DYNAMICALLY ALLOCATED SO IT
!               IS TOO LARGE AND SO WE NEED TO USE OUTWFO
                if(wfo_pick(2).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
!                   summed for reporting time step only
                    outwfo(jj,ii)=wfo_sum_p(n)
!                   reset the accumulator to 0.0 each time is is reported
                    wfo_sum_p(n)=0.0
                  end do
                  if(wfo_write_attribute_data
     *                  (xcount,ycount).ne.1)then
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!     rev. 9.1.44  Jun.  11/03  - Added Cumulative precip to the wfo file
!               ATTRIBUTE 3 - Cum. PRECIP
!               I'M DOING THIS B/C P IS NOT DYNAMICALLY ALLOCATED SO IT
!               IS TOO LARGE AND SO WE NEED TO USE OUTWFO
                if(wfo_pick(3).eq.1)then
                  do n=1,naa
                   ii=yyy(n)
                   jj=xxx(n)
!                  summed for reporting time step only
                   outwfo(jj,ii)=wfo_cum_p(n)
                  end do
                  if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif
                
!               ATTRIBUTE 4 - LOWER ZONE STORAGE
                if(wfo_pick(4).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=lzs(n)
                  end do
                  if(wfo_write_attribute_data
     *             (xcount,ycount).ne.1)then 
!                    WRITE ERROR
                     GO TO 98010
                  end if
                endif
c!               ATTRIBUTE 5 - groundwater outflow
c!               trcflg must be 'y' for this
c!               checked in rd_wfo_spec.for
c                if(wfo_pick(5).eq.1)then
c!     rev. 9.9.25  Sep.  02/14  - NK: Finally fixed the error when nbasin=0
c                  if(trcflg.eq.'y')then
c                    do n=1,naa
c                      ii=yyy(n)
c                      jj=xxx(n)
c                      l=nbasin(ii,jj)
c                      if(l.ne.0)then
c                        outwfo(jj,ii)=isoout2gw(n,l)/t
c                      else
c!                       this can happen if a watershed does not have a gauge.                      
c                        outwfo(jj,ii)=-1.0
c                      endif
cc                      outwfo(jj,ii)=isoout2gw(n,l)/t
c                    end do
cc                  else
cc                    do n=1,naa
cc                      ii=yyy(n)
cc                      jj=xxx(n)
cc                      l=nbasin(ii,jj)
cc                      if(l.ne.0)outwfo(jj,ii)=-999.0
cc                      outwfo(jj,ii)=isoout2gw(n,l)/t
cc                    end do
c                    
cc                  endif    !trcflg
c                 if(wfo_write_attribute_data
c     *             (xcount,ycount).ne.1)then 
c!                    WRITE ERROR
c                     GO TO 98010
c                  end if
c                endif

!               ATTRIBUTE 5 - groundwater outflow
                if(wfo_pick(5).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=qlz(n)
                  end do
                  if(wfo_write_attribute_data
     *                      (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif
                
!               ATTRIBUTE 6 - GRID RUNOFF
                if(wfo_pick(6).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=qr(n)
                  end do
                  if(wfo_write_attribute_data
     *                      (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!               REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
!               ATTRIBUTE 7 - OBSERVED GRID OUTFLOW
                if(wfo_pick(7).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=wfo_qhyd(n,jz)
                  end do
                  if(wfo_write_attribute_data
     *              (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!               ATTRIBUTE 8 - COMPUTED GRID OUTFLOW
                if(wfo_pick(8).eq.1)then
!                 REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
!                 calculate the mean flow for this time step
                  do n=1,naa
                    wfo_qsyn(n)=wfo_sum_qsyn(n)/float(ireport)
!                   reset to zero for next time step                    
                    wfo_sum_qsyn(n)=0.0
                  end do

                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=wfo_qsyn(n)
                  end do
                  if(wfo_write_attribute_data
     *              (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!               ATTRIBUTE 9 - WEIGHTED SWE
                if(wfo_pick(9).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=totsnw(n)
                  end do
                  if(wfo_write_attribute_data
     *             (xcount,ycount).ne.1)then 
!                    WRITE ERROR
                     GO TO 98010
                  end if
                endif

!     rev. 9.1.35  Dec. 26/02  - Added wetland & channel heights to the wfo file
!               ATTRIBUTE 10 - wetland internal water depth 
                if(wfo_pick(10).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=hwet2(n)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!               ATTRIBUTE 11 - channel depth
                if(wfo_pick(11).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=hcha2(n)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!               rev. 9.1.21  Jun.  28/02  - Added wetland storage & outflow to the wfo file
!               ATTRIBUTE 12 - wetland storage
                if(wfo_pick(12).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
c                    outwfo(jj,ii)=wstore2(n)
                    outwfo(jj,ii)=bankfull(ii,jj)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                    end if
                endif

!               ATTRIBUTE 13 - wetland OUTFLOW
                if(wfo_pick(13).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=qowet2(n)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!               ATTRIBUTE 14 - weighted evapotranspiration
                if(wfo_pick(14).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=eloss(n)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
!               ATTRIBUTE 15 - bankfull
                if(wfo_pick(15).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
c                    outwfo(jj,ii)=bankfull(ii,jj)
                    outwfo(jj,ii)=totuzs(n)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
                endif

     
                    
                    
c            write(56,*)'time= ',totaltime
c            do i=ycount,1,-1
c              write(56,56000)(bankfull(i,j),j=1,xcount)
c56000         format(<xcount>f8.0)              
c            end do
     
                    
                    
                    

!     rev. 9.9.70  Jun.  12/15  - NK: Add del_rain, and dSTRconc2  to the wfo file
                if(frcflg.eq.'y')then
!               ATTRIBUTE 12 - dlt_rain
c                if(wfo_pick(12).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=dlt_rain(n)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
c                endif

!               ATTRIBUTE 13 - dSTRconc
c                if(wfo_pick(13).eq.1)then
                  do n=1,naa
                    ii=yyy(n)
                    jj=xxx(n)
                    outwfo(jj,ii)=dSTRconc2(n)
                  end do
                  if(wfo_write_attribute_data
     *            (xcount,ycount).ne.1)then 
!                   WRITE ERROR
                    GO TO 98010
                  end if
c                endif
                endif



                npick=15      ! change this if items are added above
!               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~` 

!               ATTRIBUTE 13 - DEPRESSION STORAGE
                do i=1,classcount
                  if(wfo_pick(npick+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=d1(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR!
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 14 - DEPRESSION STORAGE (SNOW)
                do i=1,classcount
                  if(wfo_pick(npick+(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=d1fs(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 15 - SNOW WATER EQUIVALENT
                do i=1,classcount
                  if(wfo_pick(npick+2*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=snowc(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 16 - SNOW COVERED AREA
                do i=1,classcount
                  if(wfo_pick(npick+3*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=sca(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 17 - UPPER ZONE STORAGE
                do i=1,classcount
                  if(wfo_pick(npick+4*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
!     rev. 10.4.21 Apr.  21/20  = NK Add UZS deficit to wfo file = UZS(class=classcount)
                      if(i.ne.classcount)then
                        outwfo(jj,ii)=uzs(n,i)
                      else
                        outwfo(jj,ii)=totuzs(n)
                      endif    
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 18 - UPPER ZONE STORAGE (SNOW)
                do i=1,classcount
                  if(wfo_pick(npick+5*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=uzsfs(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               rev  9.1.29  Oct.  24/02  - Added q1, qint & drng to wfo file
!               ATTRIBUTE 19 - surface flow q1
                do i=1,classcount
                  if(wfo_pick(npick+6*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=q1(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 20 - surface flow q1fs (SNOW)
                do i=1,classcount
                  if(wfo_pick(npick+7*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=q1fs(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 21 - interflow qint 
                do i=1,classcount
                  if(wfo_pick(npick+8*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=qint(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 22 - interflow qintfs (SNOW)
                do i=1,classcount
                  if(wfo_pick(npick+9*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=qintfs(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 23 - recharge drng
                do i=1,classcount
                  if(wfo_pick(npick+10*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=drng(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!               ATTRIBUTE 24 - recharge drngfs (SNOW)
                do i=1,classcount
                  if(wfo_pick(npick+11*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=drngfs(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do

!     rev. 9.1.81  Apr.  04/05  - NK: added sublimation,et and etfs to wfo file
!               ATTRIBUTE 25 - evapotranspiration
!               refer to section below "reset et sums"
                do i=1,classcount
                  if(wfo_pick(npick+12*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=sum_pet(n,i)
!     *		                   /float(min(ireport,mhtot))
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do
!               ATTRIBUTE 26 - evapotranspiration (SNOW)
!               refer to section below "reset et sums"
                do i=1,classcount
                  if(wfo_pick(npick+13*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=sum_et(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do
!               ATTRIBUTE 27 - sublimation (SNOW)
!               refer to section below "reset et sums"
                do i=1,classcount
                  if(wfo_pick(npick+14*(classcount)+i).eq.1)then
                    do n=1,naa
                      ii=yyy(n)
                      jj=xxx(n)
                      outwfo(jj,ii)=sum_sublim(n,i)
                    end do
                    if(wfo_write_attribute_data
     *                (xcount,ycount).ne.1)then 
!                     WRITE ERROR
                      GO TO 98010
                    end if
                  endif
                end do
      endif

c!     "reset et sums"
c      if(mod(jz,min(ireport,mhtot)).eq.0)then
c        do ii=1,classcount
c          do n=1,naa  
c            sum_pet(n,ii)=0.0
c            sum_et(n,ii)=0.0   taken out Mar. 07/08 -nk- we want sum!!!
c            sum_sublim(n,ii)=0.0   taken out Jan. 17/11 -nk- we want sum!!!
c          end do
c        end do
c      endif


      return


!         fix fix  cound be a function
!     this should be fixed the proper way in wfo_write ...nk
!     THIS IS THE WFO ATTRIBUTE WRITING ERROR
98010 write (*,'(A)') ' '
      write (*,'(A)') '  *** FATAL ERROR ***     '
      write (*,'(A)') ' Unable to Write Attribute Data'
      STOP ' Program aborted in sub @1185'



      end subroutine write_wfo
