      SUBROUTINE aet(time,t,mon,ju)

!***********************************************************************
!    Copyright (C) 1996 by Nicholas Kouwen and Todd Neff 
        
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
! PROGRAM BY: Todd A.M. Neff,  MARCH 1996
!
! Subprogram aet.for for the Watflood program runof6.for
! Calculates actual evapotranspiration (AET) indexed to
! potential evapotranspiration (PET)
! Index is a function of degree days (tto(n)) and the
!  value of soil moisture (uzsi - upper zone storage
!  indicator)
!
!
!
!     This s/r was modified by ts and then nk added his changes back in
!     April 07

!     REV. 8.89  - Nov.  30/98  -  simplified uzs parameters
!     REV. 9.00  - Mar.   2000  -  TS: CONVERTED TO FORTRAN 90
!     REV. 9.03    Nov.   2000  - TS: ADDED WATFLOOD WETLAND ROUTING
!     rev. 9.1.26  Sep.  11/02  - fixed wetland evaporation re: uzsi
!     rev. 9.4.03  Apr.  18/07  - NK: For water ev(n,ii)=pet(n,ii)*fpet(ii)
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
!     rev. 9.9.13  Apr.  04/16  - NK: Fix water balance
!     rev. 10.1.56 Nov.  30/16  - NK: Fixed evt in AET.f to account for sca
!
!     changes made to include c&g model stuff  nk  April. 26/07
!
!  fpet2 - reduction in PET for low soil temperatures
!  full  - saturation capacity = retn(ii)/fcap (mm)
!  pwp   - permanent wilting point = fraction of the field capacity
!  sat   - soil moisture at saturation 
!  intev - interception evaporation
!  ftall - reduction in PET for tall vegetation
!          fpetmo(month,ii)=ftall(ii)*h(month,ii)/hmax     set in rdpar.for
!  ev    - is the AET
!  evt   - total AET so far - used for plotting only
!  fcap  - field capacity (% / 100)
!  ffcap - permanent wilting point (% / 100)
!  spore - saturation capacity (% / 100)
!
!***********************************************************************

	USE areacg
      use area_watflood
      implicit none
      save
      integer  :: n,ii,nnaa,mon,i,j,ios,ju
      real*4   :: uzsi,t,time,eratio,chk
      character*1 :: firstpass

! start nk addition apr. 19/07
      data firstpass/'y'/
      if(firstpass.eq.'y')then
!       write the header for evap.txt
	  if(iopt99)write(74,7403)nnprint
	  if(iopt99)write(74,7400)
!       this used to be in soilinit but needed here when read_soilinit is used
!     rev. 9.8.34  Oct.  23/12  - NK: Added sums to the resume.txt file
        if(resumflg.eq.'n')then
          do n=1,naa
	      do ii=1,classcount
              evt(n,ii)=0.0
              sum_pet(n,ii)=0.0    ! added May 14/13  nk
              sum_et(n,ii)=0.0     ! added May 14/13  nk
            end do
	    end do
	  endif
        if(frcflg.eq.'y')then
          do n=1,naa
	      do ii=1,classcount
	        icgflg(n,ii)=0
            end do
	    end do
        endif
	endif
! end nk addition apr. 19/07

!d      if(iopt.eq.2)print*,' checkpoint 10 in aet'

! ARRAY ALLOCATIONS FIRST TIME THROUGH:
      !       TS: RESET EACH TIME WE RETURN TO A GRID (for craig_gordon.for)
      if(frcflg.eq.'y'.and.isocg.eq.0)then
        allocate(dele(na,classcount),isoEconc(na,classcount),dela(na),&
       estar(na),dell(na),alphastar(na),ekin(na,classcount),&
       delstar(na,classcount),estar2H(na),alphastar2H(na),dela2H(na),&
       dell2H(na),ekin2H(na,classcount),delstar2H(na,classcount),&
       dele2H(na,classcount),iso2HEconc(na,classcount),stat=iall)
	  if(iall.ne.0) STOP 'Error allocating C&G arrays'

!       INITIALIZATION:
        do n=1,na
	    do ii=1,classcount
            ekin(n,ii)=0.0
	    end do
	  end do
	  isocg=1    ! Craig & Gordon allocation flag
	endif

      do n=nastart,naend
        i=yyy(n)
        j=xxx(n)
        if(p(i,j).le.0.0.or.firstpass.eq.'y')then
          do ii=1,classcount
!* * * * *  IMPERVIOUS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! TS: ADDED IMPERVIOUS EVAP CALCULATION: MAR 27/06
            if(ii.eq.classcount)then
! ** USE DS(n,ii) INSTEAD OF D2 TO INDICATE SURFACE STORAGE
! ** ON IMPERVIOUS AREA
              if(d2(n).gt.0.0)then
!               H2O TO EVAPORATE ON SURFACE
                ev(n,ii)=pet(n,ii)
                chk=d2(n)-ev(n,ii)
                if(chk.ge.0.0)then
!                 THERE WAS ENOUGH H2O, SUBTRACT OFF SURFACE STORAGE
                  d2(n)=d2(n)-ev(n,ii)
                else
!                 NOT ENOUGH H2O, TAKE ALL SURFACE STORAGE
                  ev(n,ii)=d2(n)
                  d2(n)=0.0
                endif
	        else
!               NO SURFACE STORAGE
                ev(n,ii)=0.0
	        endif
	        sum_et(n,ii)=sum_et(n,ii)+ev(n,ii)
	        sum_pet(n,ii)=sum_pet(n,ii)+pet(n,ii)

!c            elseif(nclass(ii).eq.'water     ')then    !elseif(ii.eq.classcount.or.nclass(ii).eq.'water     ')then
            elseif(ii.eq.ii_water)then    
            
!  * * * * *  WATER  * * * * * * * * * *  WATER  * * * * * * * * * *  WATER  * * * * * 
!             OPEN WATER EVAP IS EQUAL TO POTENTIAL RATE

!              DO NOTHING - STUFF MOVED TO LAKE_EVAP
!              DO NOTHING - STUFF MOVED TO LAKE_EVAP
!              DO NOTHING - STUFF MOVED TO LAKE_EVAP

!* * * * *  WETLAND  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! TS: MOVED FROM BELOW SO EV(N,WETLAND) IS CALCULATED AS SPECIAL CASE
! ALL OTHER CLASSES GET EV(N,II) CALC (27/03/06)
!            elseif(ii.eq.classcount-2.or.nclass(ii).eq.'wetland   '.and.
!c            elseif(ii.eq.classcount-2.and.
!c     *                 wetflg.eq.'y'.and.theta(n).gt.0.0)then
!     rev. 9.9.13  Apr.  04/16  - NK: Fix water balance
            elseif(ii.eq.classcount-2.and.wetland_flag(n))then
!                this assumes that the wetland class is before the water class
!                added theta to conditional Dec. 23/02 nk
!     rev. 10.1.48 Nov.  08/16  - NK: addet fpet(ii_water) to the wetland evaporation
! TH: trying to reduce wetland evap with water factor
                 ev(n,ii)=pet(n,ii)*fpet2(n)*fpetmo(mon,ii)*fpet(ii_water)

!     rev. 9.1.26  Sep.  11/02  - fixed wetland evaporation re: uzsi
!                DO THE WATER BALANCE FOR WETLANDS:

                 if(hwet2(n)*1000.0*theta(n).gt.ev(n,ii))then
!                   WSTORE2 IS REDUCED BY EVAPORATION - NO CONSTRAINT
!                   ftall has already been taken into account
	              qswevp(n)=ev(n,ii)*wetwid(n)*rl(n)/1000./t
!cc                    evt(n,ii)=evt(n,ii)+ev(n,ii)                  
!     rev. 10.1.56 Nov.  30/16  - NK: Fixed evt in AET.f to account for sca
                    evt(n,ii)=evt(n,ii)+ev(n,ii)*(1-sca(n,ii))

                 else
!                   half(?) AVAILABLE FREE WATER IS REMOVED
                    ev(n,ii)=hwet2(n)*1000.0*theta(n)/2.0
	              qswevp(n)=ev(n,ii)*wetwid(n)*rl(n)/1000./t
!cc                   evt(n,ii)=evt(n,ii)+ev(n,ii)
!     rev. 10.1.56 Nov.  30/16  - NK: Fixed evt in AET.f to account for sca
                    evt(n,ii)=evt(n,ii)+ev(n,ii)*(1-sca(n,ii))

                 endif

!* * * * *  ALL OTHER CLASSES  * * * * * * * * * * * * * * * * * * * * * * * * * *
! SO LONG AS THAT CLASS EXISTS IN THIS GRID	      
            else
              if(aclass(n,ii).gt.0.0)then
!               no UZS for wetlands 
                if(uzs(n,ii).ge.retn(ii))then
                   uzsi=1.
                elseif(uzs(n,ii).le.ffcap(ii))then
                   uzsi=0.
                else
                   uzsi=sqrt((uzs(n,ii)-ffcap(ii))/(retn(ii)-ffcap(ii)))
                endif

!               CALCULATE THE MAX EVAPORATION FOR THIS TIME STEP
!               BASED ON PET, UZS AND TEMPERATURE INFLUENCE
!!!!! SOME PEOPLE ARGUE THAT CANOPU EVAP IS IN ADDITION TO SOIL EVAP
!               CALC EVAP (NOT WATER, WETLAND OR IMPERVIOUS)
                ev(n,ii)=pet(n,ii)*fpet2(n)*uzsi*fpetmo(mon,ii)

                if(uzs(n,ii)-ffcap(ii).ge.ev(n,ii))then
!                  UZS IS REDUCED BY EVAPORATION - NO CONSTRAINT
                   uzs(n,ii)=uzs(n,ii)-ev(n,ii)
!cc                   evt(n,ii)=evt(n,ii)+ev(n,ii)
!     rev. 10.1.56 Nov.  30/16  - NK: Fixed evt in AET.f to account for sca
                    evt(n,ii)=evt(n,ii)+ev(n,ii)*(1-sca(n,ii))
                elseif(uzs(n,ii).le.ffcap(ii))then
!                  UZS LT WILTING POINT SO NO EVAPORATION
                   ev(n,ii)=0.0
                else
!                  ALL FREE WATER IS REMOVED
                   ev(n,ii)=uzs(n,ii)-ffcap(ii)
                   uzs(n,ii)=ffcap(ii)
!cc                   evt(n,ii)=evt(n,ii)+ev(n,ii)
!     rev. 10.1.56 Nov.  30/16  - NK: Fixed evt in AET.f to account for sca
                    evt(n,ii)=evt(n,ii)+ev(n,ii)*(1-sca(n,ii))
                endif
! 
              else      
!              ACLASS(n,ii) IS 0.0 (no area)
!              BY-PASS EVAP-SEPARATION CODE TO SAVE TIME
               if(frcflg.eq.'y')icgflg(n,ii)=0
               GOTO 100

	        endif   ! ACLASS.GT.0.0

            endif     ! CALCULATING EV(n,ii) BASED ON CLASS

!           for the glake model, if lake evaporation file is read, 
!           this sum is done in sub
!           added Mar. 12/08 -nk-
            if(.not.rd_evp_flg)then
	        sum_et(n,ii)=sum_et(n,ii)+ev(n,ii)
	        sum_pet(n,ii)=sum_pet(n,ii)+pet(n,ii)
	      endif

!     Frac model     Frac model     Frac model     Frac model     Frac model     Frac model

!     rev. 10.1.47 Nov.  08/16  - NK: Major changes in the ISO part of AET.f
            if(frcflg.eq.'y')then
! TH: modded this section for use in isoWATFLOOD
! *******************************************************************************
! PUT E-ET SEPARATION HERE.
! NEED TO IDENTIFY FRACTION THAT WILL BE APPLYING C&G TO, EVCG(n,ii)
! LEAVE AS A FUNCTION OF LANDCLASS STILL.

! NEED TO CONSIDER THAT LAND COVER TYPES CHANGE W/ H2O-SHED
! WATER       = 100% E
! IMPERMEABLE = 100% E
! BARREN      = 100% E
! INTERCEPT   = 100% E -- negligable cause storages so small
! OTHERS...   = X% E + Y% T

!              TS: RESET EACH TIME WE RETURN TO A GRID (for craig_gordon.for)
               ekin(n,ii)=0.0
	         devcg(n,ii)=0.0

               if(fpet2(n).lt.1.0)then
!                GROUND NOT WARM ENOUGH, VEG NOT DEVELOPED ENOUGH TO TRANSPIRE
!                THEREFORE ALL ET=E
	           icgflg(n,ii)=0  !0 this way it is always =1
                 evcg(n,ii)=ev(n,ii)

	         else   !  fpet2.ge.1.0
	           icgflg(n,ii)=1
!                THERE IS A MIX OF E AND T, BUT MOSTLY T
!                T IS ASSUMED 100% UNLESS THERE'S SOIL MOISTURE
!                HOW WE TREAT THE PARTITIONING DEPENDS ON LANDCOVER
                 if(ii.ge.classcount-1)then
!                WATER OR IMPERMEABLE EVAP SEPARATION
!                DOESN'T LOOP THROUGH FOR classcount, SO CALC WHEN ii=classcount
                   evcg(n,ii)=ev(n,ii)

	           !     WETLAND EVAP SEPARATION
!                WSAT=wetland degree of saturation (from 0 to 1)
                   elseif(ii.eq.classcount-2.and.wetland_flag(n))then  !  TS: fen, not bog
                      wsat(n)=hwet2(n)/(wcap(n)/wetarea(n)/theta(n))

                    if(wsat(n).ge.1.0)then
!                     WETLAND IS FLOODED, OPEN WATER: AET=PET so eratio=1
                      evcg(n,ii)=ev(n,ii)

	              else
!                     WETLAND IS NOT FLOODED: AET<PET
!                      acg(ii)=1.
!	                bcg(ii)=2.0
                      eratio=acg(ii)*wsat(n)**bcg(ii)
!                      evcg(n,ii)=intev(n,ii)+eratio*ev(n,ii)  ! TS: affect on MB?
                      evcg(n,ii)=eratio*ev(n,ii)

	              endif 
                   

	           else 
	             eratio=acg(ii)*uzsi**bcg(ii)
                   evcg(n,ii)=eratio*ev(n,ii)

                 endif   ! LANDCLASS TYPE
            
               endif     ! fpet2.ge.1.0

!              TS - called from ISOinter now - Apr 3/07
!              ISOTOPE FRACTIONATION:  
!               nnaa=s(ipr,jpr)
!               call craig_gordon(n,ii,time)

            endif        ! FRACTIONATION: frcflg.eq.'y'
! *******************************************************************************
100   CONTINUE            ! if(aclass.eq.0.0)

          end do   ! ii=1,classcount

! start nk addition apr. 19/07
! trish check:  is this needed? should there be more?
        else       ! precip this hour - no et
	    do ii=1,classcount
            ev(n,ii)=0.0
          end do
! end nk addition apr. 19/07

        endif      ! p(i,j).le.0.0
      end do       ! n=nastart,naend


!d 	if(iopt.eq.2)print*,' checkpoint 20 in aet'

!     rev. 9.9.13  Apr.  04/16  - NK: Fix water balance
      do n=nastart,naend
         do ii=1,classcount-2   
            if(wetland_flag(n).and.ii.ne.classcount-2&
              .or.wetland_flag(n)==.false.)&
              eloss(n)=eloss(n)+ev(n,ii)*aclass(n,ii)*(1.0-sca(n,ii))
         end do
!        water not added in here - done in lake_evap
         eloss(n)=eloss(n)+ev(n,classcount)*aclass(n,classcount)*(1.0-sca(n,classcount))
      end do

      if(frcflg.eq.'y')then
!        SKIP WETLANDS, BUT ADD OPEN WATER EVAPORATION:
        do n=nastart,naend
          do ii=1,classcount
	      if(wetland_flag(n).and.ii.ne.classcount-2&
              .or.wetland_flag(n)==.false.)& 
               devcg(n,ii)=ev(n,ii)*aclass(n,ii)*(1-sca(n,ii))   !TS was evcg(n,ii)
          end do
!            devcg(n,classcount)=evcg(n,classcount)*aclass(n,classcount)
!     *                               *(1-sca(n,classcount))
        end do
      endif 
!        FRACTIONATION WEIGHTED DEPTH OF EVAPORATION PER GRID
! trish check:  should this be inside ii loop?

!d	if(iopt.eq.2)print*,' checkpoint 30 in aet'

      if(iopt.ge.1)then
        if(aclass(nnprint,iiprint).lt.0.0)then
          nnaa=s(ipr,jpr)
!         write(74,7401)time,fpet2(nnaa),uzsi,tto(nnaa),intev(nnaa,1),
!     *         ev(nnaa,1),pet(nnaa,1),uzs(nnaa,1),p(yyy(nnaa),xxx(nnaa))
          if(iopt99)write(74,7401,iostat=ios)time,fpet2(nnprint),uzsi,&
              tto(nnprint),intev(nnprint,iiprint),ev(nnprint,iiprint),&
              pet(nnprint,iiprint),&
              uzs(nnprint,iiprint),p(yyy(nnprint),xxx(nnprint))
          if(ios.ne.0)then
            print*,'output conversion error in aet @ 394'
          endif
         
          if(frcflg.eq.'y')write(95,7402)time,(ev(nnprint,ii),(evcg(nnprint,ii)-intev(nnprint,ii)),ii=1,classcount)
         else
           if(firstpass.eq.'y')then   ! write this once only
             if(iopt99)write(74,*)'No area for land cover',iiprint
             if(iopt99)write(95,*)'No area for land cover',iiprint
           endif
         endif
      endif

      firstpass='n'

!d 	if(iopt.eq.2)print*,' checkpoint 40 in aet'

! FORMATS:

 7400 format('    time  fpet2(n)  uzsi   tto(n)  intev ev(n.ii) pet(n.ii) uzs(n.ii) p(n,ii)')

 7401 format(f8.0,6f8.2,f12.2,f8.2)
 7402 format(f8.0,2(<classcount>f8.2))
 7403 format('OutputForGrid# ',i8,' class # ',i8)
 7500 format(f8.0,i5,<classcount+2>(f12.4))


      RETURN

      END SUBROUTINE aet


