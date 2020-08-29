      SUBROUTINE melt(time,jan,ju)

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen and John Donald
        
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
     
!***************************************************************************
! PROGRAM BY: John DOnald - SEPT. 1991
!
! THIS SUBROUTINE ESTIMATES SNBOW ACCUMULATION & CALCULATES THE SNOWMELT
!
!
!     REV. 7.84 - Dec.  16/96 -  changed pmelt so that snowmelt only
!                                occurs on snow covered area
!     REV. 7.9  - Dec   18/96 -  took out call to read temperatures
!     REV. 8.22 - Mar.  15/96 -  glacier MF 2X when new snow=gone
!     REV. 8.24 - Apr.  07/96 -  added glacier melt multiplier gladjust
!     REV. 8.94a- Feb.  02/99 -  reset heat deficit to 0.0 on Sept.01
!     rev. 9.06    Feb.  15/01  - fixed deficit calc in melt
!     rev. 9.08 - Mar.  26/01 -  checked limits on heat def.
!     rev. 9.1.04  Oct.   4/01  - added A7 for weighting old/new sca in melt
!     rev. 9.1.19  Jun.  22/02  - Added A9 as the max heat deficit/swe ratio
!     rev. 9.5.63  Sep.  04/09  - NK: moved lapse rate from melt.f to process_temp.f
!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to par file
!     rev. 9.9.13  Apr.  04/16  - NK: Fix water balance
!     rev. 9.9.78  Sep.  16/15  - NK: Fixed wcl in melt.f
!
!	Latest edit Oct. 5/0
!
!     unit 31 = fln(11) - snowfall file
!     unit 32 = fln(12) - temperature file
!     unit 33 = fln(13) - wind file
!
!***************************************************************************

      use area_watflood
	implicit none

      CHARACTER(1) :: new,firstpass
      INTEGER      :: jan,ju,iallcnt2,iallocate,n,ii,j
	REAL*4       :: time,ttime,taold,temmmp                

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

!     ADDED TO MODULE MELTAONLY.FOR
!      REAL         :: excess1,snowf,snowt,raint    ! NK - added SEPT 20/2000

      DATA new/'f'/firstpass/'y'/
	DATA iallcnt2/0/

d      if(iopt.eq.2)print*,'in melt @ 0'

c      if(jan.eq.1.and.nnn.le.0)then
      if(firstpass.eq.'y')then
	  firstpass='n'

	  do n=1,naa     !  added March 15/07 nk
	    do ii=1,classcount
		    sum_sublim(n,ii)=0.0
		    sublim(n,ii)=-1.0    ! added May 12/13 nk
            fexcess(n,ii)=-1.0
	      ati(n,ii)=0.0
	      water(n,ii)=0.0
	    end do
	  end do

	  do ii=1,classcount
!         initialize qh   -  not used??  nk 29/12/04
          qh(ii)=-999.   ! fix - not used anywhere ?
	    qn(ii)=0.0     ! fix - not used anywhere ?
	    qe(ii)=0.0     ! fix - not used anywhere ?
	    qp(ii)=0.0     ! fix - not used anywhere ?
          refrz(ii)=0.0     ! fix - not used anywhere ?
	    qsnow(ii)=0.0     ! fix - not used anywhere ?
          qrain(ii)=0.0     ! fix - not used anywhere ?
!     rev. 10.1.99 Oct   08/17  - NK: Added error check for # sdc classes in Melt.f
          if(nsdc(ii).gt.2)then
              print*,'Error - nsdc in the par file is > 2'
              print*,'for class # ii'
              print*,'Check other classes'
              stop 'program aborted in melt @ 78'
          endif
        end do

      endif

d      if(iopt.eq.2)print*,'in melt @ 1'

!     NOTE:  SNOWC ADJUSTED FOR CHANGE IN SCA IN RUNOF5 <<<

!     FOLLOWING ADDED FOR HAMLIN CHANGES JAN 97

!     READ SNOW COVER AND MET DATA
!     RDSNOW WILL READ THE INITIAL SNOW COVER FOR EACH LAND CLASS
!     FOR JD'S OLD FILES, IT WILL ALSO READ THE MIN AND MAX TEMPS       
!     FOR THE NEW SYSTEM, THE MIN AND MAX TEMP ARE NOT USED AND HOURLY
!     TEMPS ARE EXPECTED             
!     HOWEVER, IF LONGER DELTA T'S ARE USED, A STEP FUNCTION IS ASSUMED

!     if this call is deleted, be sure to set new='t'  !!!!!
c      if(new.eq.'f')then
c        call rdsnow(time,jan,new)
c      endif
      new='t'

d      if(iopt.eq.2)print*,'in melt @ 2'

!     * * * * * * * * * * * * * *
!     SNOWMELT CALCULATION
!     * * * * * * * * * * * * * *
!     FOR EACH GRID IN THE BASIN, CALCULATE THE SNOWMELT:
!     TO PREVENT HEAT DEFICIT CARRY OVER FROM ONE YEAR TO THE NEXT:
!     FEB. 02.99   NK

	if(ju.eq.244)then
        do n=nastart,naend
          do ii=1,classcount
            def(n,ii)=0.0
          end do
        end do
      endif

!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!     rev. 9.8.28  Oct.  12/12  - NK: fixed heat deficit reset for resume
      if(dds_flag.eq.-1)then
 	  if(id.eq.1.and.time.le.1.001)then
          do n=nastart,naend
            do ii=1,classcount
              def(n,ii)=0.0
            end do
          end do
        endif
      endif

d      if(iopt.eq.2)print*,'in melt @ 3'

      do n=nastart,naend
        do ii=1,classcount
          if(aclass(n,ii).gt.0.0)then
!           WHAT'S THE TEMPERATURE?
      
!           RDTEMP READS THE TEMPERATURE ARRAY IN DELTA t INCREMENTS
!           THE OUTPUT VARIABLE IS TEMPV(n) - THE TEMPERATURE VECTOR

cc            if(new.eq.'t')then
c!             TEMPERATURE ARRAY FOR THIS DELTA T WAS READ IN MELT
c!             BUT HAS TO BE DEFINED FOR THIS GRID:
              ta(n,ii)=tempv(n)
cc            else
cc              call temper(n,ii,time,ttime)
cc              ta(n,ii)=ttime
cc            endif
c      
c            if(new.eq.'f')then
c!             IS IT SNOWING? CONVERT TO HOURLY RATE.
c              snowf(n,ii)=dsnow(n)/float(idt31)
c              if(snowf(n,ii).lt.0.0) snowf(n,ii)=0.0
c            else
c!             WE ARE USING MET FILES FOR TOTAL PRECIP
              snowf(n,ii)=0.0                       ! get rid of this variable
c            endif

!  --------------------------------------------------------------
!  WHAT TO DO? --> USE 0 DEG. C AS SNOW/RAIN BOURNDARY
!              --> IDENTIFY SNOW AND RAIN BY TEMPERATURE REGARDLESS
!                  OF WHETHER IT'S MEASURED AS RAIN OR SNOW.
!
!                  NB: SNOW AND RAIN CAN OCCUR ON THE SAME DAY
!                       (SNOW IS DAILY, RAIN IS HOURLY) 
!  -------------------------------------------------------------
!     TEMPERATURE GREATER THAN TBASE DEG C --> ALL PRECIP IS RAIN

!     REV. 7.74   May.  23/96 -  INCLUDE LAPSE RATE & ELVREF AS PART OF
!                                  .TMP FILE
            
!            if(n.eq.s(ipr,jpr).and.ii.eq.2)print*,fmadjust

!           ADJUST THE MELT FACTOR - INCREASE WITH DEGREE DAYS
            if(fmadjust.gt.0.0)then
              fmadj(n)=fmadjust*(tto(n)-ttomin(n))/100.0
              fmadj(n)=amin1(fmahigh,fmadj(n))
              fmadj(n)=amax1(fmalow,fmadj(n))
            else
              fmadj(n)=1.0
            endif

!	REV. 8.21 - Mar.  15/96 -    RAIN/SNOW CHOICE TIED TO BASE TEMP  

!		THIS FOLLOWING LINE WAS CHANGED SO THAT TA NOW ONLY HAS TO BE
!		GREATER THAN THE CURRENT BASE TEMPERATURE - FS: MARCH/97 

            if(tempv(n).gt.base(ii))then
c            if(tempv(n).gt.2)then
              deld(n,ii)=0.0
!             TOTAL PRECIP = SNOWF+R(N,II)
              raint(n,ii)=snowf(n,ii)+r(n,ii)
              snowt(n,ii)=0
              qtot(ii)=(fm(ii)*fmadj(n))*(tempv(n)-base(ii))
            else
!             TEMPERATURE LESS THAN OR EQUAL TO 0 DEG.C --> ALL PRECIP IS SNOW
!             CAN MELT IF MELT FACTOR FOR CLASS II IS LESS THAN O DEG.C
              raint(n,ii)=0 
              snowt(n,ii)=snowf(n,ii)+r(n,ii)

!             check if snow cover heat def. satisfied - if yes, ati=0
              if(def(n,ii).le.0.0) ati(n,ii)=0.0

!             CHANGE IN NEG. HEAT CONT. DUE TO ENERGY EXCHANGE
              deld(n,ii)=fmn(ii)*(ati(n,ii)-tempv(n))

!             CHANGE IN NEG. HEAT CONT. DUE TO NEW SNOW AT TA (TA -VE HERE)
              deld(n,ii)=deld(n,ii)-(snowf(n,ii)*tempv(n)/160)

!             UPDATE ATI        
              ati(n,ii)=ati(n,ii)+tipm(ii)*(tempv(n)-ati(n,ii))

!             LIMIT ATI TO -52.8 C (FROM NWSRFS)
              ati(n,ii)=amax1(-52.8,ati(n,ii))
           
!             AVAILABLE HEAT IN TERMS OF mm OF MELT:
              qtot(ii)=(fm(ii)*fmadj(n))*(tempv(n)-base(ii))
            endif

!            if(n.eq.naa/2.and.ii.eq.1)print,tto(n),ttomin(n),fmadj(n),fm(ii)
 
!           NB: MY SDC'S GO THRU 0,0 BUT NWSRFS DON'T - WHY??  

!           INITIALIZE SCA
            oldsca(n,ii)=sca(n,ii)
!            sca(n,ii)=0.0

!           CALC SNOW DEPTH & SCA
            if(snowt(n,ii).le.0.0)then
!             there is no new snow and the sca is adjusted
              dsno(n,ii)=snowc(n,ii)/rho(ii)
              
!             INTERPOLATE SCA FROM SDC
              if(dsno(n,ii).gt.sdcd(nsdc(ii),ii))then
!             See if depth of snow > 100 cover depth:
                sca(n,ii)=1.0
              elseif(dsno(n,ii).le.0.0)then
                sca(n,ii)=0.0
              else
                do j=1,nsdc(ii)-1
!                THESE CONDITIONALS OVERLAP  
                 if((dsno(n,ii).ge.sdcd(j,ii)).and.
     *              (dsno(n,ii).le.sdcd(j+1,ii)))then
	             top(n,ii)=dsno(n,ii)-sdcd(j,ii)
	             bot(n,ii)=sdcd(j+1,ii)-sdcd(j,ii)
	             sca(n,ii)=
     *                sdcsca(j,ii)+top(n,ii)/bot(n,ii)*
     *                (sdcsca(j+1,ii)-sdcsca(j,ii))
                   sca(n,ii)=amax1(0.0,sca(n,ii))
                   sca(n,ii)=amin1(1.0,sca(n,ii))
	           endif
                end do
              endif     
            else                                  ! SNOWT > 0.0
!             FOR FRESH SNOW, SNOWT > 0.0 AND SCA IS SET TO 1.0
!             AND THE SDC CURVES ARE IGNORED
!             THIS GETS RID OF INSTABILITY
              sca(n,ii)=1.0
            endif


!     rev. 0.1.04  Oct.   4/01  - added A7 for weighting old/new sca in melt
!     NOTE:
!     this results in slightly different evaporation from soil and
!     is why the new spl9 was different from the old spl8 (0.50)
!     so now this is a new parameter

!           force gradual changes in sca
            sca(n,ii)=amax1(a7*oldsca(n,ii),sca(n,ii))

!           ADDED BY FS - DEC/96 OTHERWISE SCA MIGHT BE 10E-45
            if(sca(n,ii).lt.0.000001)sca(n,ii)=0.00

!           ADJUST THE SNOW COVER  
C            if(sca(n,ii).ne.oldsca(n,ii))then
            if(sca(n,ii).ne.oldsca(n,ii).and.oldsca(n,ii).gt.0.0)then
              temmmp=wcl(n,ii)
              if(sca(n,ii).ge.0.000001)then
                snowc(n,ii)=snowc(n,ii)*oldsca(n,ii)/sca(n,ii)
                wcl(n,ii)=wcl(n,ii)*oldsca(n,ii)/sca(n,ii)
                water(n,ii)=water(n,ii)*oldsca(n,ii)/sca(n,ii)
!                jf=1
!                if(ii.eq.4)
!     *          if(writeflg(64))write(64,6201)temmmp,wcl(n,ii),oldsca(n,ii),sca(n,ii),jf
              else
	          wcl(n,ii)=0.0
!	          jf=2
!                if(ii.eq.4)
!     *          if(writeflg(64))write(64,6201)temmmp,wcl(n,ii),oldsca(n,ii),sca(n,ii),jf
              endif
            endif

!           CHECK IF SCA HAS CHANGED IN THIS GRID SINCE LAST TIME STEP
            if(sca(n,ii).lt.oldsca(n,ii))then
!             SCA IS GOING DOWN AND SHOULD BE LESS THAN 1
!             REDUCE D1 AND UZS TO REFLECT GREATER BARE GROUND AREAS
!             BASED ON CONSERVATION OF WATER 
!             ADJUST D1 AND UZS TO COMPENSATE FOR CHANGE IN SNOW COVERED GRND
              d1(n,ii)=d1(n,ii)+
     *                (oldsca(n,ii)-sca(n,ii))*(d1fs(n,ii)-d1(n,ii))/
     *                                          (1.0-sca(n,ii))
              uzs(n,ii)=uzs(n,ii)+
     *                (oldsca(n,ii)-sca(n,ii))*(uzsfs(n,ii)-uzs(n,ii))/
     *                                          (1.0-sca(n,ii))
!     fix fix  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
c              snowc(n,ii)=snowc(n,ii)*sca(n,ii)/oldsca(n,ii)
c              wcl(n,ii)=wcl(n,ii)*oldsca(n,ii)/sca(n,ii)

            elseif(sca(n,ii).gt.oldsca(n,ii))then
!             SCA IS INCREASING (AND SHOULD BE LARGER THAN 0)
              d1fs(n,ii)=d1fs(n,ii)+
     *                 (sca(n,ii)-oldsca(n,ii))*(d1(n,ii)-d1fs(n,ii))/
     *                                           sca(n,ii)
              uzsfs(n,ii)=amax1(uzsfs(n,ii),0.000001) !fixed ts 02/12/03

              uzsfs(n,ii)=uzsfs(n,ii)+
     *                 (sca(n,ii)-oldsca(n,ii))*(uzs(n,ii)-uzsfs(n,ii))/
     *                                           sca(n,ii)

!     fix fix  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
c              snowc(n,ii)=snowc(n,ii)*sca(n,ii)/oldsca(n,ii)
c              wcl(n,ii)=wcl(n,ii)*oldsca(n,ii)/sca(n,ii)
            endif

!     This section works fine when the snow is continuous but totally screws up when
!     you read in new snow course data on the fly.
!     It was fixed with the  .and. .....  added to the conditional

c            if(sca(n,ii).le.0.00001)then
            if(sca(n,ii).eq.0.0.and.oldsca(n,ii).gt.0.0)then
!             THERE IS NO SNOW BUT WHEN IT RETURNS, WATER WILL COME FROM
!             RECEEDING PACK
              d1fs(n,ii)=0.0 

c              uzsfs(n,ii)=0.02
              uzsfs(n,ii)=0.00

	        sca(n,ii)=0.0   ! new 22/04/02
!     fix fix  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            endif

c            if(sca(n,ii).ge.0.99999)then
            if(sca(n,ii).eq.1.0.and.oldsca(n,ii).lt.1.0)then
!             THERE IS NO BARE GROUND, BUT WHEN BARE GROUND OPENS UP, WATER
!             WILL COME FROM PACK 
              d1(n,ii)=0.0
!              uzs(n,ii)=0.0001
              sca(n,ii)=1.0    ! new 22/04/02
!     fix fix  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
            endif

!           ADD NEW SNOW TO COVER    

!           WHEN SCA = 0.0, ie. FIRST SNOW ON BARE GROUND, LET IT COVER THE
!           WHOLE GROUND. IN THE NEXT TIME IT WILL BE DISTRIBUTED

            snowc(n,ii)=snowc(n,ii)+snowt(n,ii)

!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
!           sublim_factor is always set to max out at 0.5
!           so sublim < snowc
            if(sublimflg.eq.1)then
	        sublim(n,ii)=sublim_factor(ii)*snowt(n,ii)
	      else
!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to par file
!             sublimate only when there is snow that is not melting 
	        if(snowc(n,ii).ge.sublim_rate(ii).and.def(n,ii).gt.0.0)then
	          sublim(n,ii)=sublim_rate(ii)
	        else
	          sublim(n,ii)=0.0
	        endif
	      endif

            snowc(n,ii)=snowc(n,ii)-sublim(n,ii)


!           eloss is used in the watbal only

            if(ii.eq.classcount.or.nclass(ii).eq.'water     ')then
!             do nothing - water class is not included in water balance
	      else
              eloss(n)=eloss(n)+sublim(n,ii)*aclass(n,ii)*sca(n,ii)
	      endif

            sum_sublim(n,ii)=sum_sublim(n,ii)+sublim(n,ii)
     *                         *aclass(n,ii)*sca(n,ii) !added Jul. 24/12 nk

!           UPDATE SNOW COVER HEAT DEFECIT
            def(n,ii)=def(n,ii)+deld(n,ii)
            if(def(n,ii).lt.0.0)def(n,ii)=0.0

!           LIMITS ON HEAT DEFICIT OF SNOW COVER
!           NWSRFS--> MAX. DEFICIT = .33*SWE
!           OTHERWISE THERE MIGHT BE TOO MUCH HEAT DEFICIT FOR ALL THE 
!           SNOW TO MELT - FEB 2/99 NK  

!     rev. 9.1.19  Jun.  22/02  - Added A9 as the max heat deficit/swe ratio
            if(def(n,ii).gt.a9*snowc(n,ii))then
              def(n,ii)=a9*snowc(n,ii)
            endif

!           REV. 7.77   Jul.  02/96 -     FIXED SNOW REDISTRIBUTION

            if(snocap(ii).gt.0.0)then
!             REDISTRIBUTION IS TURNED OFF IF SNOCAP = -VE
!             CHECK CAPACITY - REDIST TO RECEIVING CALSS IDUMP(ii)
!             FIRST CHECK CAPACITY:
              if(snowc(n,ii).gt.snocap(ii).and.idump(ii).ne.ii)then
!              CALCULATE EXCESS SNOW: 
               extra(ii)=snowc(n,ii)-snocap(ii)

!              REV 8.99H  MAY 9/00  NK
!              ADJUST FOR GRU SIZE:
               if(aclass(n,idump(ii)).gt.0.0)then
                 extra(ii)=extra(ii)*aclass(n,ii)/aclass(n,idump(ii))
               endif
!              CHECK IF RECEIVING CLASS HAS ROOM:

               if(snocap(idump(ii)).ge.snowc(n,idump(ii))+extra(ii))then
                 snowc(n,ii)=snocap(ii)
!                REDISTRIBUTE CLASS=IDUMP(ii)
                 snowc(n,idump(ii))=snowc(n,idump(ii))+extra(ii)
                 extra(ii)=0.0
               else
!                SNOW STAY'S PUT AND EXCESS IS SET BACK TO 0
                 extra(ii)=0.0
               endif
              endif
            endif

! NOTE: NO HEAT REDISTRIBUTION. OK???

!           ADJUST QTOT, RPOBG, ROSN FOR SCA
!           R(n,ii) CHANGED TO RAINT
!           RAINT INCLUDES SNOW & RAIN PRECIP
            robg(ii)=raint(n,ii)
            rosn(ii)=raint(n,ii)

!           CALCULATE QNET
!           qtot is heat input -ve or +ve
            qnet(ii)=qtot(ii)-def(n,ii)
            if(qnet(ii).lt.0)then
!             temp below base temp & heat is lost -> increase def
              def(n,ii)=def(n,ii)-qtot(ii)    ! ???????????????????????
!     rev. 9.08 - Mar.  26/01 -  checked limits on heat def.
              if(def(n,ii).lt.0.0)def(n,ii)=0.0
              qnet(ii)=0.0
            else
              def(n,ii)=0.0
            endif

!           CALCULATE SNOWMELT
            if(qnet(ii).gt.snowc(n,ii))then
!             MELT ALL SNOW     
              smelt(ii)=snowc(n,ii)
            else       
!             CAN'T MELT ALL SNOW
              if(qnet(ii).gt.0.0)then
!               RIPE PACK, NO HEAT DEFECIT 
                smelt(ii)=qnet(ii)
              else  
!               NO SNOW MELT       
                smelt(ii)=0.0
              endif              
            endif

!           NO MATTER, WE HAVE TO MAKE SURE THE WATER BALANCE IS OK
!           EVEN IF THE HEAT CALCS ARE WRONG!

!           BEGINNING OF HAMLIN CHANGES JAN/97
!           HEAT AND WATER BALANCE OF SNOW COVER      
!           --> SUBTRACT ANY EXCESS HEAT DEFICIT NOT YET SATISIFED
!           - COULD HAPPEN IF LIGHT RAIN ON VERY COOL DAY
!           def(n,ii) -  the heat deficit in mm melt
!           wcl(n,ii) -  free water in the pack I think???

!           water     - all available liquid water
            water(n,ii) = smelt(ii)+rosn(ii)+wcl(n,ii)-def(n,ii)

!           wlmax     -  the max amount the pack can hold
            wlmax(n,ii) = snowc(n,ii)*whcl(ii)

            snowc(n,ii)=snowc(n,ii)-smelt(ii)

            if(def(n,ii).gt.water(n,ii)) then
!             DEFICIT IS GREATER THAN THE LIQUID WATER
!             SNOWPACK IS NOT RIPE ===> SATISFY THE HEAT DEFICIT FIRST
!             FREEZE ALL THE LIQUID WATER AND ADD IT TO THE SNOW PACK
!             THERE IS NO FREE WATER AND NO EXCESS
!             THIS MEANS THAT THE MELT WAS UNDONE
!     rev. 9.9.78  Sep.  16/15  - NK: Fixed wcl in melt.f
c              wcl(n,ii)=0.0
              excess(ii)=0.0
!             THIS IS CHANGED. MELT WATER IS NOT ADDED BACK TO THE PACK
!             BUT THE -VE HEAT HAS BEEN TAKEN OUT OF DEF SO IT IS PUT BACK
              snowc(n,ii)=snowc(n,ii)+rosn(ii)+wcl(n,ii)             !nk
              wcl(n,ii)=0.0 !TH: wcl not being added to pack, moved wcl=0 to after snowc recalc        
!              def(n,ii)=def(n,ii)+rosn(ii)+wcl(n,ii)                 !nk
!     rev. 9.06    Feb.  15/01  - fixed deficit calc in melt
              def(n,ii)=def(n,ii)-rosn(ii)-wcl(n,ii)                 !nk

!             added line below March 26/01  nk
!     rev. 9.08 - Mar.  26/01 -  checked limits on heat def.
              if(def(n,ii).lt.0.0)def(n,ii)=0.0

            else
!             DEFICIT IS LESS THAN LIQUID WATER
!             SNOWPACK IS RIPE AND DEFICIT IS WIPED OUT
!             WHAT EVER DEFICIT THERE WAS, FROZE SOME WATER AND IT IS
!             ADDED TO THE PACK
!     rev. 9.9.78  Sep.  16/15  - NK: Fixed wcl in melt.f
c              water(n,ii)=water(n,ii)-def(n,ii)  !TH: ran into v. small numerical difference in isosnow, def removed twice.
!             REMOVE DEFICIT BY REFREEZING SOME LIQUID WATER & ADD TO PACK
              snowc(n,ii)=snowc(n,ii)+def(n,ii)
              def(n,ii)=0.0
!             CHECK IF WATER HOLDING CAPACITY HAS BEEN EXCEEDED (ie.RUNOFF) 
              if(water(n,ii).gt.wlmax(n,ii)) then
!               NOT ALL THE LIQUID WATER CAN LEAVE THE PACK 
!               WATER AVAILABLE FOR RUNOFF
                excess(ii)=water(n,ii)-wlmax(n,ii)
                wcl(n,ii)=wlmax(n,ii)
              else
!               WATER HOLDING CAPACITY HAS BEEN EXCEEDED (ie.NO RUNOFF)  
                wcl(n,ii)=water(n,ii)
                excess(ii)=0.0
              endif
            endif
            snowc(n,ii)=amax1(0.0,snowc(n,ii))
            wlmax(n,ii)=snowc(n,ii)*whcl(ii)
!           CHECK IF SNOWPACK IS RIPE THEN ATI=0.0 (ANDERSON'73 3-6(a))
            if(def(n,ii).le.0.0)then
              ati(n,ii)=0.0           
            endif
       
!           ++++++++++++++++++++++++++++++++++++++  FS BEGIN
!           THIS IS WHERE JD PUT TOGETHER SNOWPACK RUNOFF(EXCESS)
!           AND RAIN ON BARE GROUND (ROBG), SO FOR MY RUNOF5.FOR TO
!           WORK THERE IS A NEW VARIABLE CALLED FEXCESS(ii) WHICH 
!           HOLDS THE SNOWPACK RUNOFF.  THIS IS HOW I WOULD DO IT:
            fexcess(n,ii)=excess(ii)

!           REV. 8.22 - Mar.  15/96 - GLACIER MF2X WHEN NEW SNOW=GONE
!           REV. 8.24 - Apr.  07/96 - ADDED GLACIER MELT MULT USE

!           IF WE ARE LOOKING AT THE GLACIER CLASS 
!           THEN IF THERE IS NO SNOW ON TOP OF THE GLACIER 
!           (ie. snowc.lt.1.0) THEN FORCE THERE TO BE RUNOFF 2 TIMES > 
!           THE POTENTIAL MELT (NICK CAME UP WITH THE 2 TIES FIGURE)
!           FS: FEB/97
            if(sca(n,ii).le.1.0.and.
     *         nclass(ii)(1:7).eq.'glacier'
     *                 .and.qtot(ii).gt.0.0)then
              glmelt(n)=qtot(ii)*gladjust
            endif

!           ++++++++++++++++++++++++++++++++++++++  FS END
!           * * * * *  END OF HAMLIN CHANGES JAN 97  * * * * *   

          endif                         ! if(aclass(n,ii).gt.0.0)
        end do                          ! do ii=1, classcount

!       SUM TOTAL SNOWMELT AND RAINFALL CONTRIBUTION TO RUNOFF
!       FROM EACH GRU (BEFORE STREAMFLOW ROUTING)

!       -----------------  OUTPUT RESULTS ---------------------

!      if(iopt.eq.2)print*,'in melt @ 4'
   
!       FRANK CHANGES OPTION TO GO TO SNOUT1

c        if(maxn.eq.2001.and.time.eq.7000)then
c          call snout1(time,n)
c        endif

!       if(iopt.eq.2)print*,'in melt @ 4a'

!     rev. 10.1.57 Dec.  06/16  - NK: Added snwNN.txt files for iopt > 0
        if(iopt.ge.2)then
          if(numa.eq.0)then
            if(n.eq.nnprint)then
!             WRITE HEADER
              ii=iiprint           ! ii must be .le. classcount!!!!!
              ii=11           ! ii must be .le. classcount!!!!!
              if(jan.le.1)then
                do ii=1,classcount
                  write(800+ii,6000)ii
                  write(800+ii,6003)n
                  write(800+ii,6001) 
                end do  
              endif
!             WRITE TO SNW.csv FOR LAND CLASS II=IIOUT       
              do ii=1,classcount
                write(800+ii,6002) 
     *	        time,tempv(n),raint(n,ii),snowt(n,ii),
     *          sca(n,ii),qn(ii),qe(ii),qh(ii),
     *          qp(ii),qtot(ii),def(n,ii),qnet(ii),robg(ii),
     *          rosn(ii),smelt(ii),excess(ii),refrz(ii),wcl(n,ii),
     *          snowc(n,ii),qsnow(ii),qrain(ii),r(n,ii),water(n,ii),
     *          glmelt(n),sublim(n,ii),snowf(n,ii),fexcess(n,ii)
              end do
            endif
          endif
!         if(iopt.eq.2)print*,'in melt @ 5'
        endif
      end do                             ! do n=1, naa


c	if(totaltime.gt.8760)then
c	  do n=1,na
c	    do i=1,ii
c  	      if(snowc(n,ii).eq.0.0)then
c              penalty(n)=1.0
c	      endif
c	    end do
c	  end do
c	endif


d      if(iopt.eq.2)print*,'in melt @ return'

! FORMATS

 6000   format('land class =',i5)    
 6001   format(
     *       '   time,    temp, rainfall, snowfall,    sca,     qn,  ',
     *       '     qe,      qh,       qp,     qtot,    def,   qnet, ',
     *       '   robg,    rosn,    smelt,   excess,  refrz,    wcl, ',
     *       ' snowc ,  qsnow ,   qrain ,    r(  ),  water,  glmelt,',
     *       '  sublm,   snowf,  fexcess')
 6002   format(27(f8.2,',')) 
 6003   format(' vector element=',i5) 
 5700 format(' ',16f5.1)
 6201 format(4e12.5,3x,i1)


      RETURN

      END SUBROUTINE melt

