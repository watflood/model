      SUBROUTINE tracer0(iz,jz,time,t,jan,tdum)

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
! SUB-BASIN TRACER
! This subroutine is designed to route isotope tracers through the
! WATFLOOD model to separate out flows from individual sub-basins.  
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
!     REV. 9.1.50  Jan.  27/04 - TS: Separated tracer codes & fixed IBN tracking
!     REV. 9.1.51  Apr.  25/05 - TS: Added IF(lnxt.eq.ll) statement **not tested**
!
!*****************************************************************************

c      USE area1
c	USE area2
c      USE area3
c      USE area4
c      USE area5
c      USE area6
c	USE area10
c      USE areaet
cc	USE areamelt
c      USE areatrc
c	USE areawet

      use area_watflood
	implicit none

!*****NOTE: THIS PROGRAM ASSUMES (FOR NOW) ZERO BACKGROUND CONCENTRATIONS
!*****OF THE TRACERS.  FROM RUN-TO-RUN THOUGH, TRACER MASSES ARE CONSERVED.
!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      INTEGER :: n,lll,l,iz,jz,ii,jj,i,j,jan,ios,inxt,jnxt,lnxt,
     *           ll,llold
      REAL*4  :: wi,wold,t,woldn,woldp,isooutold,decay,time
     	REAL*4  :: tdum



      if(index.eq.0)then
        do n=1,na
          res(n)=0
        end do

        do n=1,na
          lll=next(n)
          do l=1,noresv
            if(yyy(n).eq.ires(l).and.xxx(n).eq.jres(l)) res(n)=l
          end do
        end do
      endif
            
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                                                               * *
!     TRACER 0 => SUB-BASIN TRACER                              * *
!                                                               * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     INITIALIZE AND RESET VALUES FOR EACH TIME STEP
      do n=1,na
!       INITIALIZE TRACERS:
        isoin1IBN(n)=isoin2IBN(n)
        isoin2IBN(n)=0.0
        isoout1IBN(n)=isoout2IBN(n)
        isoout2IBN(n)=0.0
        isostore1IBN(n)=isostore2IBN(n)
        isostore2IBN(n)=0.0
        if(wetland_flag(n))then
          isoin1wet(n,l)=isoin2wet(n,l)
          isoin2wet(n,l)=0.0
          isoout1wet(n,l)=isoout2wet(n,l)
          isoout2wet(n,l)=0.0
          isowstore1(n,l)=isowstore2(n,l)
          isowstore2(n,l)=0.0
	  endif
      end do

!     ^^^^^^^^^^^ CHANNEL ROUTE THE TRACER ^^^^^^^^^^^^^^^^^^
!     isoin    = isotope mass in kg
!     isoout   = isotope mass out in kg
!     isostore = isotope in storage in kg
!     isoconc  = isotope concentration in the unit in kg/m^3
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        do n=1,naa
          lll=next(n)
          inxt=yyy(lll)
          jnxt=xxx(lll)
	    lnxt=nbasin(inxt,jnxt)

          i=yyy(n)
          j=xxx(n)
	    ll=nbasin(i,j)

!         RUN TRACER CODE ONLY IF NOT IN A LAKE!!
          if(ireach(n).le.0.0 .OR.
     *       ireach(n).gt.0.and.res(n).gt.0) then


!           CHECK TO ENSURE WE'RE IN THE SAME BASIN AS PREVIOUS...
!           IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!           IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
            if(llold.eq.ll.and.ll.ne.0)then
 
!            WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
             if(slope(n).gt.0.0)then

               if(store2(n).ne.0.0)then
                  isooutold=1.0e+25

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                ROUTE TRACER THROUGH WETLAND FIRST (IF THERE ARE WETLANDS):
                 if(wetland_flag(n)==.true..and.nbasin(i,j).eq.l)then

!                  SET INITIAL CONCENTRATION AS 1 kg/m^3] UNITS:  = [KG]
!                  TRACE MASS OF TRACER ENTERING WETLAND, ADDING TRACER
!                  ONLY IF GRID=N IS IN THE CURRENT SUB-BASIN=L:
                   isoin2wet(n,ll)=isoin2wet(n,ll)+1.0*qr(n)*t   

                   call tracewet(iz,jz,time,n,ll,t)

                   if(qowet2(n).gt.0.0)then
!                    SPECIFY INITIAL MASS FOR CHANNEL FROM WETLAND:
                     isoin2IBN(n)=isoin2IBN(n)+isooutwet

                   else   ! qowet2(n) is -ve
!                   SPECIFY INITIAL MASS IN WETLAND (+MASS FROM CHANNEll)
!                   NOTE: -VE CAUSE QOWET2 IS -VE BUT MASS IS ADDED!
                    isoin2wet(n,ll)=isooutwet-isoconcIBN(n)*qowet2(n)*t
                    isoin2IBN(n)=isoin2IBN(n)+isoconcIBN(n)
     *                                           *qowet2(n)*t
	             endif

	           else
!                  EITHER NO WETLANDS, OR NOT IN RIGHT SUB-BASIN THEREFORE
!                  NO TRACER IS BEING ADDED ANYWAY...SKIP WETLAND ROUTING.
!                  TRACE MASS OF TRACER ENTERING CHANNEL FROM SUBBASIN=L
                   isoin2IBN(n)=isoin2IBN(n)+1.0*qr(n)*t
	           endif

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! WETLANDS AND QOWET2 IS -VE, SO REBALANCE THE WETLAND CONC'N WITH CHANNEL
!                INFLOW. THEN ROUTE WHATEVER'S LEFT IN THE CHANNEL.
                 if(wetland_flag(n)==.true..and.qowet2(n).lt.0
     *                               .and.ll.eq.llold)then ! QOWET2 IS -ve
!                  CHANNEL FEEDS THE WETLAND, THEREFORE MIX CHANNEL CONC'N
!                  WITH WETLAND AND UPDATE WETLAND CONC'N FOR NEXT TIME STEP.
!                  ROUTE TRACER FROM CHANNEL TO WETLAND:
                   isowetold=1.0e+25
!                  ROUTE THE TRACER THROUGH THE WETLAND STORAGE:
                   call tracewet(iz,jz,time,n,ll,t)
                 endif   ! QOWET IS -VE

! NO WETLANDS AND QOWET2 IS +VE (OR ALREADY REROUTED THRU WETLAND)
!                CHANNEL ROUTING OR WETLAND->CHANNEL ROUTING
!                ROUTE TRACER THROUGH CHANNEL - CONVERGENCE LOOP:

!                UPDATE THE TRACER MASS IN CHANNEL n
                 isostore2IBN(n)=isostore1IBN(n)+
     *                (isoin1IBN(n)+isoin2IBN(n)-
     *                 isoout1IBN(n)-isoout2IBN(n))/2.
	           isostore2IBN(n)=amax1(isostore2IBN(n),0.00001)

!                CONVERGENCE LOOP:
                 do ijk=1,50

!                  UPDATE THE CONCENTRATION
                   isoconcIBN(n)=isostore2IBN(n)/store2(n)

!                  COMPUTE THE TRACER MASS LEAVING THE GRID...
!                  EVAP LOSS TO PRESERVE FLOW BALANCES.
!                  THIS SHOULD BE BACK DATED ONE TIME STEP
                   isoout2IBN(n)=isoconcIBN(n)*qo2(n)*t
     *                            -isoconcIBN(n)*strloss(n)*t
	             isoout2IBN(n)=amax1(isoout2IBN(n),0.00001)

!                  UPDATE THE TRACER MASS IN CHANNEL n
                   isostore2IBN(n)=isostore1IBN(n)+
     *                  (isoin1IBN(n)+isoin2IBN(n)-
     *                   isoout1IBN(n)-isoout2IBN(n))/2.
	             isostore2IBN(n)=amax1(isostore2IBN(n),0.00001)

!                  CHECK FOR CONVERGENCE - WE NEED THIS TO ITERATE ON ISOOUT2()
                   if(abs(isoout2IBN(n)-isooutold).lt.0.01*isooutold) 
     *                GOTO 01
                   isooutold=isoout2IBN(n)
                 end do
   01            CONTINUE
          
               else
                 isostore2IBN(n)=0.0
                 isoout2IBN(n)=0.0
               end if                       ! if(store2(n).ne.0.0)
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  
!              ISOTOPES ARE CONSERVATIVE!
               decay = 0.0

!              UPDATE MASS BEING INPUT INTO NEXT GRID CELL IF THE NEXT 
!              GRID IS NOT A "DUMMY GRID" AND ITS IN THE SAME BASIN
               if(lnxt.ne.0.and.lnxt.eq.ll)then
                  isoin2IBN(lll)=(1-(decay/100))*
     *                           (isoin2IBN(lll)+isoout2IBN(n))
	         else
                  isoin2IBN(lll)=0.0
	         endif

!              ACCUMULATE MASSES FOR MASS BALANCE CHECK AT END OF RUN:
               massin(n)=(isoin1IBN(n)+isoin2IBN(n))/2.
               massout(n)=(isoout1IBN(n)+isoout2IBN(n))/2.
               masstore(n)=(isostore2IBN(n)-isostore1IBN(n))
               ISOdelta(n)=massin(n)-massout(n)
 

 !             CHECK OVERALL MASS BALANCE:
 !             MASS IN - MASS OUT = MASS STORED (IF EVERYTHING WORKS!) 
               if(iz.ne.jz.and.n.eq.nnprint)then
                 sqerr=(ISOdelta(n)-masstore(n))**2.
                 write(91,8000)time,massin(n),massout(n),ISOdelta(n),
     *                         masstore(n),sqerr
               endif

             endif      ! SLOPE>=0


            else
!             DON'T COUNT CONTRIBUTION FROM THIS GRID CAUSE IT DRAINS
!             TO A DIFFERENT BASIN (NOT CONTRIBUTING TO SUB-BASIN OUTLET)
	        CONTINUE
	      endif      ! BASIN CHECK


          else




!           IN A LAKE, DON'T RUN TRACER!
	      CONTINUE 





          endif

!         SAVE THE BASIN NUMBER FOR THE NEXT ITERATION
          llold=ll

!	if(n.eq.39) print*,n,ll,llold,lnxt,isoout2IBN(n)
!	if(n.eq.39) pause

	end do       ! GRID NO. LOOP
     


!     WRITE OUT DATA TO TRACER.CSV FILE
!      if(jan.eq.3)then    ! Prints every 0.25 hours
      if(iz.ne.jz)then     ! Prints every hour
!       WRITE TO THE FILE ONLY WHEN YOU ARE AT THE FLAGGED OUTPUT
!       GRID=NNPRINT. 
!       UNITS: isoout2IBN(n,ll)=[kg]/[s]=[kg/s]=[m^3/s]
!       PRINTS @ GRID NNPRINT FLOWS FOR BASIN L ONLY...
!!	  nnprint=39
        ii=yyy(nnprint)
        jj=xxx(nnprint)
	  ll=nbasin(ii,jj)
        write(90,9001)time,
     *    ll,nnprint,qo2(nnprint),qr(nnprint),
     *    isoout2IBN(nnprint)/t,isoconcIBN(nnprint)
      endif


! FORMATS:
 8000 FORMAT(f10.3,256(',',f15.5))
 9000 FORMAT(f10.3,<no>(2(',',I15),4(',',f15.5)))
 9001 FORMAT(f10.3,2(',',I15),4(',',f15.5))




      RETURN
      END SUBROUTINE tracer0