      SUBROUTINE tracerGW(iz,jz,time,t,jan,tdum,jjz,jjzold)

!***********************************************************************
!    Copyright (C) 2016 by Nicholas Kouwen, Tricia Stadnyk and Tegan Holmes  
        
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
     
!*****************************************************************************
! GROUNDWATER TRACER
! This subroutine is designed to route isotope tracers through the
! WATFLOOD model to track the amount of runoff from groundwater.  
! It will do so based on a mixed cell model (continuity).
! Adapted from LFL's sediment transport subroutine.
!
! Created:  May 2003    - NK
! Modified: Summer 2003 - TS
!
!     REV. 9.1.42  May   31/03 - GW Tracer module added - first try
!     REV. 9.1.43  June  18/03 - Expanded tracer module
!     REV. 9.1.44  June  19/03 - TS: Added sub-basin tracer =0
!     REV. 9.1.45  July  29/03 - TS: Added landcover tracer =1
!     REV. 9.1.46  Aug.  18/03 - TS: Added rain tracer (qstream) =3
!     REV. 9.1.47  Aug.  19/03 - TS: Added sub-basin/landcover tracer =2
!     REV. 9.1.48  Aug.  26/03 - TS: Added flow-type & snowmelt tracers =4/5
!     REV. 9.1.49  Nov.  23/03 - TS: Added wetlands to GW Tracer + Wetland Tracer
!     REV. 9.1.50  Jan.  27/04 - TS: Separated tracer codes
!     REV. 9.8.77  July  2013  - TH: Added lake routing, diversions & nudging
!     REV.         Feb.  4/16  - TH: Replaced convergence loop with direct solution
!
!*****************************************************************************

      use area_watflood
      implicit none	
	
!     FROM RUN-TO-RUN THOUGH, TRACER MASSES ARE CONSERVED.
!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      INTEGER :: n,m,lll,l,iz,jz,i,j,jan,inxt,jnxt,lnxt
	INTEGER :: rbin,resnum,jjz,jjzold
      REAL*4  :: t,isooutold,decay,time
     	REAL*4  :: tdum,roe

!     changed ro to roe to avoid conflict with wq  nk  May 22/07

      roe=1.0

!     AT VERY START OF WATFLOOD RUN, SET res(n) EQUAL TO l IF LAKE, 0 OTHERWISE
!     LOOP FROM FIRST GRID, n=1, TO LAST NON-OUTLET GRID, naa
      res_set: if(index==0)then
        do n=1,naa
          res(n)=0
        end do

        do n=1,naa
          lll=next(n)
          do l=1,noresv
            if(yyy(n)==ires(l).and.xxx(n)==jres(l)) res(n)=l
          end do
        end do
      end if res_set



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *	
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                                                               * *
!     TRACER 100 => GW TRACER (orig: nick kouwen)               * *
!                                                               * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     INITIALIZE AND RESET VALUES FOR EACH TIME STEP
      initialize: do n=1,naa
        i=yyy(n)          !integer holder for y-coord
        j=xxx(n)          !integer holder for x-coord
        l=nbasin(i,j)     !basin # of current grid
	  lll=next(n)       !n # of next grid in flow path
        resnum=ireach(n)  !lake # of current grid
	  rbin=ireach(lll)  !lake # of next grid in flow path

        if(rbin/=0) then
          isolakeGW(rbin)=0.0  !isolakeGW is the GW inflow to a lake
          isolakest(rbin)=0.0  !isolakest is the total open water 
!                                   evaporation from a lake
	  endif

!       SHIFT VALUES UP TIME STEP
        if(l/=0)then
          isoin1GW(n,l)=isoin2GW(n,l)
          isoin2GW(n,l)=0.0
          isoout1GW(n,l)=isoout2GW(n,l)
          isoout2GW(n,l)=0.0
          isostore1GW(n,l)=isostore2GW(n,l)
          isostore2GW(n,l)=0.0
	    if(wetflg=='y')then
            isoin1wet(n,l)=isoin2wet(n,l)
            isoin2wet(n,l)=0.0
            isoout1wet(n,l)=isoout2wet(n,l)
            isoout2wet(n,l)=0.0
            isowstore1(n,l)=isowstore2(n,l)
            isowstore2(n,l)=0.0
	    endif
	  endif
      end do initialize

!     IF THERE ARE DIVERSIONS, GIVE THE DIVERTED FLOW AN ASSUMED TRACER MASS      
      divert: do m=1,nodivert
        n=gridgive(m)  
!       gridgive holds the grid number the diversion flows into
        i=yyy(n)
        j=xxx(n)
        l=nbasin(i,j)
        isoin2GW(n,l)=isoin2GW(n,l)+0.75*qr(n)*t
!       qr(n) holds qdivert for this time step
      end do divert

!       ^^^^^^^^^^^  TRACER ROUTING  ^^^^^^^^^^^^^^^^^^^^^^^^^^
!       isoin    = isotope mass in kg
!       isoout   = isotope mass out in kg
!       isostore = isotope in storage in kg
!       isoconc  = isotope concentration in the unit in kg/m^3
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
      route: do n=1,naa
        lll=next(n)
        i=yyy(n)
        j=xxx(n)
        l=nbasin(i,j)
        resnum=ireach(n)
	  rbin=ireach(lll)
        inxt=yyy(lll)          !y-coord of next grid
        jnxt=xxx(lll)          !x-coord of next grid
	  lnxt=nbasin(inxt,jnxt) !basin # of next grid in flow path
	  
	  
!	  WHEN IN LAKE, ADD OPEN WATER EVAPORATION FOR GRID 
!       TO TOTAL EVAPORATION FOR THE LAKE 
        if(resnum>0)
     *    isolakest(resnum)=isolakest(resnum)+strloss(n)

!       RUN FOR LAND AND LAKE OUTLET ONLY
        if(ireach(n)<=0.0 .OR.
     *     resnum>0.and.res(n)>0) then

!       WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
!       IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!       IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
        in_basin: if(slope(n)>0.0.and.l/=0)then

          positive_storage: if(store2(n)>0.0)then
             
!       RESET HOLDERS FOR OLD ISOOUT             
             isooutold=1.0e+25
	       oldiso=isoout2GW(n,l)

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!            ROUTE TRACER THROUGH WETLAND FIRST (IF THERE ARE WETLANDS):
             wetlands: if(wetland_flag(n)==.true.)then

!              SET INITIAL CONCENTRATION AS 1 kg/m^3]
!              UNITS:  = [KG]
!              TRACE MASS OF GW ENTERING WETLAND FROM QLZ
               isoin2wet(n,l)=isoin2wet(n,l)+1.0*qlz(n)*t

               

               wet_pos_outflow: if(qowet2(n)>=0.0)then
!                SPECIFY INITIAL MASS FOR CHANNEL FROM WETLAND:
                 call tracewet(iz,jz,time,n,l,t)
                 isoin2GW(n,l)=isoin2GW(n,l)+isooutwet

               else   ! qowet2(n) is -ve
!                SPECIFY INITIAL MASS IN WETLAND (+MASS FROM CHANNEL)
!                NOTE: -VE CAUSE QOWET2 IS -VE BUT MASS IS ADDED!

                 isoin2wet(n,l)=isoin2wet(n,l)
     *                       -isoconcGW(n,l)*qowet2(n)*t
                 isoin2GW(n,l)=isoin2GW(n,l)
     *                        +isoconcGW(n,l)*qowet2(n)*t
                 call tracewet(iz,jz,time,n,l,t)
     
	         endif wet_pos_outflow

	       else  !if wetland_flag(n)=='n' ->no wetland routing
!              SET INITIAL CONCENTRATION AS 1 kg/m^3
!              UNITS:  = [KG]
!              ADD THE MASS OF GW ENTERING CHANNEL FROM QLZ
	         isoin2GW(n,l)=isoin2GW(n,l)+1.0*qlz(n)*t
	         
!	      AT LAKE OUTLET, MAKE INFLOW EQUAL THE ACCUMULATED INFLOWS TO WHOLE LAKE
	      if(resnum>0.and.res(n)>0) isoin2GW(n,l)=isolakeGW(resnum)
	        
	       endif wetlands

!              UPDATE THE CONCENTRATION IN THE CHANNEL
 
               if(resnum>0.and.res(n)>0)then !if at lake outlet grid
                 isostore2GW(n,l)=(isostore1GW(n,l)+(isoin1GW(n,l)+
     *              isoin2GW(n,l)-isoout1GW(n,l))/2.)
     *              /(1.0+(qo2(n)*coeff(n)+isolakest(resnum))*t
     *              /2.0/store2(n))
	           isostore2GW(n,l)=amax1(isostore2GW(n,l),0.00001)
	           isoconcGW(n,l)=isostore2GW(n,l)/store2(n)
                 isoout2GW(n,l)=isoconcGW(n,l)*qo2(n)*t*coeff(n)
     *                       +isoconcGW(n,l)*isolakest(resnum)*t
     	           isoout2GW(n,l)=amax1(isoout2GW(n,l),0.00001)
     	           
               else
                 isostore2GW(n,l)=(isostore1GW(n,l)+(isoin1GW(n,l)+
     *            isoin2GW(n,l)-isoout1GW(n,l))/2.)
     *              /(1.0+(qo2(n)*coeff(n)+strloss(n))*t
     *              /2.0/store2(n))
	           isostore2GW(n,l)=amax1(isostore2GW(n,l),0.00001)
	           isoconcGW(n,l)=isostore2GW(n,l)/store2(n)
                 isoout2GW(n,l)=isoconcGW(n,l)*qo2(n)*t*coeff(n)
     *                       +isoconcGW(n,l)*strloss(n)*t
     	           isoout2GW(n,l)=amax1(isoout2GW(n,l),0.00001)
     	         end if
             
          else   ! follow what's done in rerout
!           TS: added May 21/08 s.t tracer works when store2<0	
            if(isoin2GW(n,l)<0.0)then
              isoout2GW(n,l)=isoin2GW(n,l)/2.0
            else
              isoout2GW(n,l)=0.0
	      endif
            isostore2GW(n,l)=isostore1GW(n,l)+(isoin1GW(n,l)+
     *            isoin2GW(n,l)-isoout1GW(n,l)-isoout2GW(n,l))/2.
	      isostore2GW(n,l)=amax1(isostore2GW(n,l),0.00001)

          end if positive_storage
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   

!         UPDATE MASS BEING INPUT INTO NEXT GRID CELL IF THE NEXT GRID IS NOT A "DUMMY GRID"
!         SINCE ISOOUT INCLUDES EVAP LOSSES, REMOVE EVAP SO IT IS NOT PASSED TO NEXT GRID
          if(lnxt/=0) then !IF NO LAKES, CONTINUE TO PASS
            if(nopt(l)==2.and.n==iflowgrid(l))then
            !IF IT'S A NUDGED STATION, USE THE NUDGED FLOW
              isoin2GW(lll,lnxt)=isoin2GW(lll,lnxt)
     *            +qhyd(l,iz+1)*t*isoconcGW(n,l)
            else
              if(resnum>0.and.res(n)>0)then
               isoin2GW(lll,lnxt)=isoin2GW(lll,lnxt)
     *           +isoout2GW(n,l)-isoconcGW(n,l)*isolakest(resnum)*t
              else 
               isoin2GW(lll,lnxt)=isoin2GW(lll,lnxt)
     *           +isoout2GW(n,l)-isoconcGW(n,l)*strloss(n)*t
            end if
           endif
	    endif

!         ADD MASS OUTFLOW OF CURRENT GRID AS INFLOW TO LAKE (NEXT GRID) (cms)
!         BUT ONLY IF ON THE EDGE OF A LAKE (NOT IN A LAKE)
          if(rbin>0) then
           !IF IT'S A NUDGED STATION, USE THE NUDGED FLOW
           if(nopt(l)==2.and.n==iflowgrid(l))then
            isolakeGW(rbin)=isolakeGW(rbin)
     *            +qhyd(l,iz+1)*t*isoconcGW(n,l)
           else 
            isolakeGW(rbin)=isolakeGW(rbin)+isoout2GW(n,l)
           endif
	    endif
	    	    

!         CALCULATE MASSES FOR MASS BALANCE CHECK AT END OF RUN:
          massin(n)=(isoin1GW(n,l)+isoin2GW(n,l))/2.
          massout(n)=(isoout1GW(n,l)+isoout2GW(n,l))/2.
          masstore(n)=(isostore2GW(n,l)-isostore1GW(n,l))
          ISOdelta(n)=massin(n)-massout(n)

 !        CHECK OVERALL MASS BALANCE:
 !        MASS IN - MASS OUT = MASS STORED (IF EVERYTHING WORKS!) 
          MBcheck: if(iz/=jz.and.n==nnprint)then
            sqerr=(ISOdelta(n)-masstore(n))**2
            write(91,8000)time,massin(n),massout(n),ISOdelta(n),
     *                    masstore(n),sqerr

            wetMB: if(wetflg=='y')then  !MB check for wetlands 
              sqerr=(wISOdelta(n)-wmasstore(n))**2
              write(94,8000)time,wmassin(n),wmassout(n),wISOdelta(n),
     *                      wmasstore(n),sqerr
	      endif wetMB
	    endif MBcheck

         end if in_basin           ! SLOPE>=0 AND L.NE.0

	  endif
          
      end do route            ! GRID NO. LOOP
      
!     WRITE OUT DATA TO TRACER.CSV FILE, MAKE SURE WE ARE WITHIN
!     CALCULATION AREA OF BASIN (L NOT= 0)
!      if(jan.eq.3)then    ! Prints every 0.25 hours
      if(iz/=jz)then ! Prints every hour
        do l=1,no
!         STORES THE GRID NO. AS FXN OF THE GAUGE(BSN) NO.
          nn(l)=s(iy(l),jx(l))
	  end do

!       WRITE TO THE FILE ONLY WHEN YOU ARE AT THE FLAGGED OUTPUT
!       GRID=NNPRINT. 
!       UNITS: isoout2IBN(n,ll)=[kg]/[s]=[kg/s]=[m^3/s]
!       PRINTS OUT EVERY SUB-BASIN @ SUB-BASIN GAUGE LOCATION
        if(wetflg=='y') write(93,9101)time,
     *    (l,nn(l),qlz(nn(l)),qowet2(nn(l)),qswevp(nn(l)),
     *     isoout2wet(nn(l),l)/(roe*t),isoconcwet(nn(l),l),l=1,no)

!       added by NK Mar. 17/06
        if(mod(iz,kt)==0)then
          write(90,99887)
     *      totaltime,(qo2(nn(l))*isoconcGW(nn(l),l),l=1,no)
	  endif

!        i=yyy(nnprint)
!        j=xxx(nnprint)
!        l=nbasin(i,j)
!        if(wetflg.eq.'y'.and.dds_flag.ne.1) write(204,20400)time,
!     *    (isoin2GW(nnprint,l)-isoout2wet(nnprint,l))/(roe*t),
!     *    isoin2GW(nnprint,l)/qo2(nnprint)/(roe*t),
!     *    isoout2wet(nnprint,l)/(roe*t),isoconcwet(nnprint,l),
!     *    isowstore2(nnprint,l),
!     *    isoout2GW(nnprint,l)/(roe*t),isoconcGW(nnprint,l),
!     *    isostore2GW(nnprint,l)

      endif

!      filename(90)='..\simout\tracer.csv'
!      filename(91)='..\simout\tracerMB.csv'     ! added Oct.30/03 TS
!      filename(93)='..\simout\tracerWET.csv'    ! added Dec.01/03 TS
!      filename(94)='..\simout\tracerWETMB.csv'  ! added Dec.01/03 TS


! FORMATS:
 8000 FORMAT(f10.3,256(',',f15.5))
c 9000 FORMAT(f10.3,<no>(2(',',I15),4(',',f15.5)))
c 9100 FORMAT(f10.3,<no>(2(',',I15),5(',',f15.5)))
c 9101 FORMAT(f10.3,<no>(2(',',I15),5(',',f15.5)))
!               fixed for general format  nk  Jan. 17/01
 9000 FORMAT(f10.3,<no>(2(',',I15),4(',',g15.5)))
 9100 FORMAT(f10.3,<no>(2(',',I15),5(',',g15.5)))
 9101 FORMAT(f10.3,<no>(2(',',I15),5(',',g15.5)))
10000 FORMAT('   time   ,',
     *<no>('         basin#,          grid#,',
     *'        qlz    ,       qswevp ,        ISOoutQ,',
     *'        ISOconc,'))
11000 FORMAT('    time  ,        IN     ,       OUT     ,',
     *'        DELTA  ,       STORE   ,       SSE     ,')
20400 FORMAT(f10.3,999(',',f15.5))
99887 FORMAT(f10.3,<no>((',',f15.5)))
12345 FORMAT(f10.3,<no>(3(',',f15.5))) 



      RETURN
      END SUBROUTINE tracerGW      
