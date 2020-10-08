      SUBROUTINE tracer5(iz,jz,time,t,jan,tdum,jjz,jjzold)

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
! SNOW-MELT TRACER
! This subroutine is designed to route isotope tracers through the
! WATFLOOD model to track the amount of runoff from e4ach flow-type includ. melt H20.  
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
!
!*****************************************************************************


      use area_watflood
	implicit none


!*****NOTE: THIS PROGRAM ASSUMES (FOR NOW) ZERO BACKGROUND CONCENTRATIONS
!*****OF THE TRACERS.  FROM RUN-TO-RUN THOUGH, TRACER MASSES ARE CONSERVED.
!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      INTEGER :: n,lll,l,iz,jz,ii,jj,i,j,jan,ios,inxt,jnxt,lnxt
	INTEGER :: rbin,jjz,jjzold
      REAL*4  :: wi,wold,t,woldn,woldp,decay,time,
     *           GWchk,GWfschk,SWchk,SWfschk,IFchk,IFfschk,isooutold1,
     *           isooutold2,isooutold3,isooutold4,isooutold5,isooutold6,
     *           isooutold7,isooutold8
     	REAL*4  :: tdum,isoratio,roe
c     	REAL*4  :: tdum,isoratio,ro


      roe=1.0    ! not used in this s/r  nk

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
!     TRACER 5 => SNOWMELT TRACER (SW,IF,GW)                      * *
!                                                               * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *	
      if(itrace.eq.5)then

!       INITIALIZE AND RESET VALUES FOR EACH TIME STEP
        do n=1,na
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)
	    lll=next(n)
	    if(lll.ne.0) rbin=ireach(lll)

          if(rbin.ne.0 .and. iz.ne.jz) then
!           RESET TO '0' ONLY IF TIME HAS CHANGED
!           OTHERWISE, KEEP ACCUMULATING THE MASS!
!  	      print*,'initializing isolakeGW',jjz,jjzold,rbin
! 	      pause
            isolakeGW(rbin)=0.0
	    endif

	    isosumQ(n)=0.0
	    isosumQfs(n)=0.0

	    if(l.ne.0)then
            isoin1SW(n,l)=isoin2SW(n,l)
            isoin2SW(n,l)=0.0
            isoout1SW(n,l)=isoout2SW(n,l)
            isoout2SW(n,l)=0.0
            isostore1SW(n,l)=isostore2SW(n,l)
            isostore2SW(n,l)=0.0
            isoin1SWfs(n,l)=isoin2SWfs(n,l)
            isoin2SWfs(n,l)=0.0
            isoout1SWfs(n,l)=isoout2SWfs(n,l)
            isoout2SWfs(n,l)=0.0
            isostore1SWfs(n,l)=isostore2SWfs(n,l)
            isostore2SWfs(n,l)=0.0
            isoin1IF(n,l)=isoin2IF(n,l)
            isoin2IF(n,l)=0.0
            isoout1IF(n,l)=isoout2IF(n,l)
            isoout2IF(n,l)=0.0
            isostore1IF(n,l)=isostore2IF(n,l)
            isostore2IF(n,l)=0.0
            isoin1IFfs(n,l)=isoin2IFfs(n,l)
            isoin2IFfs(n,l)=0.0
            isoout1IFfs(n,l)=isoout2IFfs(n,l)
            isoout2IFfs(n,l)=0.0
            isostore1IFfs(n,l)=isostore2IFfs(n,l)
            isostore2IFfs(n,l)=0.0
            isoin1GW(n,l)=isoin2GW(n,l)
            isoin2GW(n,l)=0.0
            isoout1GW(n,l)=isoout2GW(n,l)
            isoout2GW(n,l)=0.0
            isostore1GW(n,l)=isostore2GW(n,l)
            isostore2GW(n,l)=0.0
            isoin1GWfs(n,l)=isoin2GWfs(n,l)
            isoin2GWfs(n,l)=0.0
            isoout1GWfs(n,l)=isoout2GWfs(n,l)
            isoout2GWfs(n,l)=0.0
            isostore1GWfs(n,l)=isostore2GWfs(n,l)
            isostore2GWfs(n,l)=0.0
	      isoin1LZS(n,l)=isoin2LZS(n,l)
            isoin2LZS(n,l)=0.0
	      isoout1LZS(n,l)=isoout2LZS(n,l)
	      isoout2LZS(n,l)=0.0
	      isoLZS1(n,l)=isoLZS2(n,l)
	      isoLZS2(n,l)=0.0
            isoin1LZSfs(n,l)=isoin2LZSfs(n,l)
            isoin2LZSfs(n,l)=0.0
            isoout1LZSfs(n,l)=isoout2LZSfs(n,l)
            isoout2LZSfs(n,l)=0.0
            isoLZS1fs(n,l)=isoLZS2fs(n,l)
	      isoLZS2fs(n,l)=0.0
	      if(wetland_flag(n))then
              isoin1wet(n,l)=isoin2wet(n,l)
              isoin2wet(n,l)=0.0
              isoout1wet(n,l)=isoout2wet(n,l)
              isoout2wet(n,l)=0.0
              isowstore1(n,l)=isowstore2(n,l)
              isowstore2(n,l)=0.0
              isoin1fswet(n,l)=isoin2fswet(n,l)
              isoin2fswet(n,l)=0.0
              isoout1fswet(n,l)=isoout2fswet(n,l)
              isoout2fswet(n,l)=0.0
              isowstore1fs(n,l)=isowstore2fs(n,l)
              isowstore2fs(n,l)=0.0
              isoin1SWwet(n,l)=isoin2SWwet(n,l)
              isoin2SWwet(n,l)=0.0
              isoout1SWwet(n,l)=isoout2SWwet(n,l)
              isoout2SWwet(n,l)=0.0
              isowstore1SW(n,l)=isowstore2SW(n,l)
              isowstore2SW(n,l)=0.0
              isoin1SWfswet(n,l)=isoin2SWfswet(n,l)
              isoin2SWfswet(n,l)=0.0
              isoout1SWfswet(n,l)=isoout2SWfswet(n,l)
              isoout2SWfswet(n,l)=0.0
              isowstore1SWfs(n,l)=isowstore2SWfs(n,l)
              isowstore2SWfs(n,l)=0.0
              isoin1IFwet(n,l)=isoin2IFwet(n,l)
              isoin2IFwet(n,l)=0.0
              isoout1IFwet(n,l)=isoout2IFwet(n,l)
              isoout2IFwet(n,l)=0.0
              isowstore1IF(n,l)=isowstore2IF(n,l)
              isowstore2IF(n,l)=0.0
              isoin1IFfswet(n,l)=isoin2IFfswet(n,l)
              isoin2IFfswet(n,l)=0.0
              isoout1IFfswet(n,l)=isoout2IFfswet(n,l)
              isoout2IFfswet(n,l)=0.0
              isowstore1IFfs(n,l)=isowstore2IFfs(n,l)
              isowstore2IFfs(n,l)=0.0
	      endif
          endif
	  end do

!       ^^^^^^^^^^^  TRACER ROUTING  ^^^^^^^^^^^^^^^^^^^^^^^^^^
!       isoin    = isotope mass in kg
!       isoout   = isotope mass out in kg
!       isostore = isotope in storage in kg
!       isoconc  = isotope concentration in the unit in kg/m^3
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
        do n=1,naa
          lll=next(n)
          inxt=yyy(lll)
          jnxt=xxx(lll)
	    lnxt=nbasin(inxt,jnxt)
	    rbin=ireach(lll)
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)

!         RUN TRACER CODE ONLY IF NOT IN A LAKE
          if(ireach(n).le.0.0 .OR.
     *       ireach(n).gt.0.and.res(n).gt.0) then

!       ^^^^^^^^^^^^^ GW ROUTE THE TRACER ^^^^^^^^^^^^^^^^^
!         TO GET MASS OF TRACER ENTERING LZS FROM DRNG
          if(lzs(n).ne.0.0.and.l.ne.0)then
            isooutold7=1.0e+25
            isooutold8=1.0e+25
		  oldiso7=isoout2LZS(n,l)
            oldiso8=isoout2LZSfs(n,l)

!           CALC MASS OF TRACER GOING INTO GW STORAGE, SET INITIAL
!           CONCENTRATION AS 1 kg/m^3 * m^3/s * s = KG:
            isoin2LZS(n,l)=isoin2LZS(n,l)+1.0*qdrng(n)*t
            isoin2LZSfs(n,l)=isoin2LZSfs(n,l)+1.0*qdrngfs(n)*t

!           UPDATE THE TRACER MASS IN GW STORAGE
            isoLZS2(n,l)=isoLZS1(n,l)+(isoin1LZS(n,l)+
     *           isoin2LZS(n,l)-isoout1LZS(n,l)-isoout2LZS(n,l))/2.
	      isoLZS2(n,l)=amax1(isoLZS2(n,l),0.00001)

            isoLZS2fs(n,l)=isoLZS1fs(n,l)+(isoin1LZSfs(n,l)+
     *       isoin2LZSfs(n,l)-isoout1LZSfs(n,l)-isoout2LZSfs(n,l))/2.
	      isoLZS2fs(n,l)=amax1(isoLZS2fs(n,l),0.00001)

!           CONVERGENCE LOOP:
            do ijk=1,50

!             UPDATE CONCENTRATION USING INSTANT MIXING!! 
!             CONC = KG/m^3
!             LZS CONV: mm / 1000m/mm * m^2 = m^3
              isoconcLZS(n,l)=isoLZS2(n,l)/(lzs(n)/1000.*frac(n)*al**2)
              isoconcLZSfs(n,l)=isoLZS2fs(n,l)
     *                                    /(lzs(n)/1000.*frac(n)*al**2)
!	        if(isoconcLZS(n,l).lt.0.0001) isoconcLZS(n,l)=0.0
!	        if(isoconcLZSfs(n,l).lt.0.0001) isoconcLZSfs(n,l)=0.0

!             COMPUTE THE TRACER MASS LEAVING THE LZS
!             THIS SHOULD BE BACK DATED ONE TIME STEP
!             KG/m^3 * m^3/s * s = KG
!             ASSUME: NO EVAP IN SOIL COLUMN!!
              isoout2LZS(n,l)=isoconcLZS(n,l)*qlz(n)*t
              isoout2LZSfs(n,l)=isoconcLZSfs(n,l)*qlz(n)*t

!             TS: ADDED RELAXATION FOR ITERATING TO CODE TO HELP CONVERGENCE (04/1/06)
              isowt=amax1(.5,real(ijk/51.))
              isoout2LZS(n,l)=(1.0-isowt)*isoout2LZS(n,l)+isowt*oldiso7
	        oldiso7=isoout2LZS(n,l)
              isoout2LZSfs(n,l)=(1.0-isowt)*isoout2LZSfs(n,l)
     *                          +isowt*oldiso8
	        oldiso8=isoout2LZSfs(n,l)

!             UPDATE THE TRACER MASS IN GW STORAGE
              isoLZS2(n,l)=isoLZS1(n,l)+(isoin1LZS(n,l)+
     *             isoin2LZS(n,l)-isoout1LZS(n,l)-isoout2LZS(n,l))/2.
	        isoLZS2(n,l)=amax1(isoLZS2(n,l),0.00001)

              isoLZS2fs(n,l)=isoLZS1fs(n,l)+(isoin1LZSfs(n,l)+
     *       isoin2LZSfs(n,l)-isoout1LZSfs(n,l)-isoout2LZSfs(n,l))/2.
	        isoLZS2fs(n,l)=amax1(isoLZS2fs(n,l),0.00001)

!             CHECK FOR CONVERGENCE - WE NEED THIS TO ITERATE ON ISOOUT2()
              if(abs(isoout2LZS(n,l)-isooutold7).lt.0.00001
     *          .and.abs(isoout2LZSfs(n,l)-isooutold8).lt.0.00001) 
     *          GOTO 110
              isooutold7=isoout2LZS(n,l)
              isooutold8=isoout2LZSfs(n,l)

            end do

  110       CONTINUE

          elseif(l.ne.0)then
           isoLZS2(n,l)=0.0
           isoout2LZS(n,l)=0.0
           isoLZS2fs(n,l)=0.0
           isoout2LZSfs(n,l)=0.0
          endif



!       ^^^^^^^^^^^^^ CHANNEL ROUTE THE TRACER ^^^^^^^^^^^^^^^^^
!         ROUTED QLZ CONTRIBUTION ADDED TO CHANNEL AND ROUTED AGAIN.
!         WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
!         IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!         IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
          if(slope(n).gt.0.0.and.l.ne.0)then
            if(store2(n).ne.0.0)then
               isooutold1=1.0e+25
               isooutold2=1.0e+25
               isooutold3=1.0e+25
               isooutold4=1.0e+25
               isooutold5=1.0e+25
               isooutold6=1.0e+25
		     oldiso1=isoout2GW(n,l)
               oldiso2=isoout2SW(n,l)
               oldiso3=isoout2IF(n,l)
		     oldiso4=isoout2GWfs(n,l)
               oldiso5=isoout2SWfs(n,l)
               oldiso6=isoout2IFfs(n,l)

!  * * * * * * * * * PUT WETLAND CODE HERE  * * * * * * * * * * * * * * * * 
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!              ROUTE TRACER THROUGH WETLAND FIRST (IF THERE ARE WETLANDS):
               if(wetland_flag(n)==.true.)then

!                SET INITIAL CONCENTRATION AS 1 kg/m^3]
!                UNITS:  = [KG]
!                TRACE MASS ENTERING WETLAND FROM FLOWPATH
                 isoin2wet(n,l)=isoin2wet(n,l)+isoout2LZS(n,l)
                 isoin2SWwet(n,l)=isoin2SWwet(n,l)+1.0*sumq1(n)*t
                 isoin2IFwet(n,l)=isoin2IFwet(n,l)+1.0*sumqint(n)*t
                 isoin2fswet(n,l)=isoin2fswet(n,l)+isoout2LZSfs(n,l)
                 isoin2SWfswet(n,l)=isoin2SWfswet(n,l)+1.0*sumq1fs(n)*t
                 isoin2IFfswet(n,l)=isoin2IFfswet(n,l)
     *                                              +1.0*sumqintfs(n)*t

!      if(n.eq.46)write(801,'(3f10.3)')
!     *   time,isoin2IFwet(n,l),isoin2IFfswet(n,l)
!      if(n.eq.67) pause 'routing Tracer 5 through wetland @1'

                 call tracewet(iz,jz,time,n,l,t)

                 if(qowet2(n).gt.0.0)then
!                  SPECIFY INITIAL MASS FOR CHANNEL FROM WETLAND:
                   isoin2GW(n,l)=isoin2GW(n,l)+isooutwet1
                   isoin2SW(n,l)=isoin2SW(n,l)+isooutwet2
                   isoin2IF(n,l)=isoin2IF(n,l)+isooutwet3
                   isoin2GWfs(n,l)=isoin2GWfs(n,l)+isooutwet4
                   isoin2SWfs(n,l)=isoin2SWfs(n,l)+isooutwet5
                   isoin2IFfs(n,l)=isoin2IFfs(n,l)+isooutwet6

!      if(n.eq.67) pause 'routing Tracer 5 thru wetland @2, qowet2 +ve'

	           else   ! qowet2(n) is -ve
!                  SPECIFY INITIAL MASS IN WETLAND (+MASS FROM CHANNEL)
!                  NOTE: -VE CAUSE QOWET2 IS -VE BUT MASS IS ADDED!
                   isoconcGW(n,l)=isostore2GW(n,l)/store2(n)
                   isoconcSW(n,l)=isostore2SW(n,l)/store2(n)
                   isoconcIF(n,l)=isostore2IF(n,l)/store2(n)
                   isoconcGWfs(n,l)=isostore2GWfs(n,l)/store2(n)
                   isoconcSWfs(n,l)=isostore2SWfs(n,l)/store2(n)
                   isoconcIFfs(n,l)=isostore2IFfs(n,l)/store2(n)

                   isoin2wet(n,l)=isooutwet1-isoconcGW(n,l)*qowet2(n)*t
                   isoin2SWwet(n,l)=isooutwet2
     *                             -isoconcSW(n,l)*qowet2(n)*t
                   isoin2IFwet(n,l)=isooutwet3
     *                             -isoconcIF(n,l)*qowet2(n)*t
                   isoin2fswet(n,l)=isooutwet4
     *                            -isoconcGWfs(n,l)*qowet2(n)*t
                   isoin2SWfswet(n,l)=isooutwet5
     *                            -isoconcSWfs(n,l)*qowet2(n)*t
                   isoin2IFfswet(n,l)=isooutwet6
     *                            -isoconcIFfs(n,l)*qowet2(n)*t

                   isoin2GW(n,l)=isoin2GW(n,l)+isoconcGW(n,l)
     *                                        *qowet2(n)*t
                   isoin2SW(n,l)=isoin2SW(n,l)+isoconcSW(n,l)
     *                                        *qowet2(n)*t
                   isoin2IF(n,l)=isoin2IF(n,l)+isoconcIF(n,l)
     *                                        *qowet2(n)*t
                   isoin2GWfs(n,l)=isoin2GWfs(n,l)+isoconcGWfs(n,l)
     *                                        *qowet2(n)*t
                   isoin2SWfs(n,l)=isoin2SWfs(n,l)+isoconcSWfs(n,l)
     *                                        *qowet2(n)*t
                   isoin2IFfs(n,l)=isoin2IFfs(n,l)+isoconcIFfs(n,l)
     *                                        *qowet2(n)*t

!      if(n.eq.67) pause 'routing Tracer 5 thru wetland @3, qowet2 -ve'

	           endif

	         else
!                SET INITIAL CONCENTRATION AS 1 kg/m^3
!                UNITS:  = [KG]
!                TRACE MASS OF GW ENTERING CHANNEL FROM QLZ
	           isoin2GW(n,l)=isoin2GW(n,l)+isoout2LZS(n,l)
                 isoin2SW(n,l)=isoin2SW(n,l)+1.0*sumq1(n)*t
                 isoin2IF(n,l)=isoin2IF(n,l)+1.0*sumqint(n)*t
	           isoin2GWfs(n,l)=isoin2GWfs(n,l)+isoout2LZSfs(n,l)
                 isoin2SWfs(n,l)=isoin2SWfs(n,l)+1.0*sumq1fs(n)*t
                 isoin2IFfs(n,l)=isoin2IFfs(n,l)+1.0*sumqintfs(n)*t

!      if(n.eq.67)pause 'routing Tracer 5 @4, no wetlands'


	         endif

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! WETLANDS AND QOWET2 IS -VE, SO REBALANCE THE WETLAND CONC'N WITH CHANNEL
!            INFLOW. THEN ROUTE WHATEVER'S LEFT IN THE CHANNEL.
             if(wetland_flag(n)==.true.
     *                       .and.qowet2(n).lt.0)then ! QOWET2 IS -ve
!              CHANNEL FEEDS THE WETLAND, THEREFORE MIX CHANNEL CONC'N
!              WITH WETLAND AND UPDATE WETLAND CONC'N FOR NEXT TIME STEP.
!              ROUTE TRACER FROM CHANNEL TO WETLAND:
               isowetold1=1.0e+25
               isowetold2=1.0e+25
               isowetold3=1.0e+25
               isowetold4=1.0e+25
               isowetold5=1.0e+25
               isowetold6=1.0e+25
!              ROUTE THE TRACER THROUGH THE WETLAND STORAGE:
               call tracewet(iz,jz,time,n,l,t)

!      if(n.eq.67) pause 're-routing Tracer 5 @5, qowet2 -ve'

             endif   ! QOWET IS -VE

!  * * * * * * * * * END WETLAND CODE HERE  * * * * * * * * * * * * * * * * 
! NO WETLANDS AND QOWET2 IS +VE (OR ALREADY REROUTED THRU WETLAND)
!            CHANNEL ROUTING OR WETLAND->CHANNEL ROUTING
!            ROUTE TRACER THROUGH CHANNEL - CONVERGENCE LOOP:

      
!              UPDATE THE TRACER MASSES IN CHANNEL
               isostore2GW(n,l)=isostore1GW(n,l)+(isoin1GW(n,l)+
     *           isoin2GW(n,l)-isoout1GW(n,l)-isoout2GW(n,l))/2.
     	         isostore2GW(n,l)=amax1(isostore2GW(n,l),0.00001)

               isostore2SW(n,l)=isostore1SW(n,l)+(isoin1SW(n,l)+
     *           isoin2SW(n,l)-isoout1SW(n,l)-isoout2SW(n,l))/2.
	         isostore2SW(n,l)=amax1(isostore2SW(n,l),0.00001)

               isostore2IF(n,l)=isostore1IF(n,l)+(isoin1IF(n,l)+
     *          isoin2IF(n,l)-isoout1IF(n,l)-isoout2IF(n,l))/2.
	         isostore2IF(n,l)=amax1(isostore2IF(n,l),0.00001)

               isostore2GWfs(n,l)=isostore1GWfs(n,l)+(isoin1GWfs(n,l)+
     *           isoin2GWfs(n,l)-isoout1GWfs(n,l)-isoout2GWfs(n,l))/2.
	         isostore2GWfs(n,l)=amax1(isostore2GWfs(n,l),0.00001)

               isostore2SWfs(n,l)=isostore1SWfs(n,l)+(isoin1SWfs(n,l)+
     *          isoin2SWfs(n,l)-isoout1SWfs(n,l)-isoout2SWfs(n,l))/2.
	         isostore2SWfs(n,l)=amax1(isostore2SWfs(n,l),0.00001)

               isostore2IFfs(n,l)=isostore1IFfs(n,l)+(isoin1IFfs(n,l)+
     *          isoin2IFfs(n,l)-isoout1IFfs(n,l)-isoout2IFfs(n,l))/2.
	         isostore2IFfs(n,l)=amax1(isostore2IFfs(n,l),0.00001)


!              CONVERGENCE LOOP:
               do ijk=1,50

!                UPDATE THE CONCENTRATIONS IN THE CHANNEL:
                 isoconcGW(n,l)=isostore2GW(n,l)/store2(n)
                 isoconcSW(n,l)=isostore2SW(n,l)/store2(n)
                 isoconcIF(n,l)=isostore2IF(n,l)/store2(n)
                 isoconcGWfs(n,l)=isostore2GWfs(n,l)/store2(n)
                 isoconcSWfs(n,l)=isostore2SWfs(n,l)/store2(n)
                 isoconcIFfs(n,l)=isostore2IFfs(n,l)/store2(n)

!                COMPUTE THE TRACER MASSES LEAVING THE GRID...
!                EVAP LOSS TO PRESERVE FLOW BALANCES.
!                THIS SHOULD BE BACK DATED ONE TIME STEP
                 isoout2GW(n,l)=isoconcGW(n,l)*qo2(n)*t*coeff(n)
     *                         -isoconcGW(n,l)*strloss(n)*t
     	           isoout2GW(n,l)=amax1(isoout2GW(n,l),0.00001)

                 isoout2SW(n,l)=isoconcSW(n,l)*qo2(n)*t*coeff(n)
     *                         -isoconcSW(n,l)*strloss(n)*t
     	           isoout2SW(n,l)=amax1(isoout2SW(n,l),0.00001)

                 isoout2IF(n,l)=isoconcIF(n,l)*qo2(n)*t*coeff(n)
     *                         -isoconcIF(n,l)*strloss(n)*t
     	           isoout2IF(n,l)=amax1(isoout2IF(n,l),0.00001)

                 isoout2GWfs(n,l)=isoconcGWfs(n,l)*qo2(n)*t*coeff(n)
     *                           -isoconcGWfs(n,l)*strloss(n)*t
     	           isoout2GWfs(n,l)=amax1(isoout2GWfs(n,l),0.00001)

                 isoout2SWfs(n,l)=isoconcSWfs(n,l)*qo2(n)*t*coeff(n)
     *                           -isoconcSWfs(n,l)*strloss(n)*t
     	           isoout2SWfs(n,l)=amax1(isoout2SWfs(n,l),0.00001)

                 isoout2IFfs(n,l)=isoconcIFfs(n,l)*qo2(n)*t*coeff(n)
     *                           -isoconcIFfs(n,l)*strloss(n)*t
     	           isoout2IFfs(n,l)=amax1(isoout2IFfs(n,l),0.00001)


!                TS: ADDED RELAXATION FOR ITERATING TO CODE TO HELP CONVERGENCE (03/11/05)
                 isowt=amax1(.5,real(ijk/51.))

                 isoout2GW(n,l)=(1.0-isowt)*isoout2GW(n,l)+isowt*oldiso1
	           oldiso1=isoout2GW(n,l)

                 isoout2SW(n,l)=(1.0-isowt)*isoout2SW(n,l)+isowt*oldiso2
	           oldiso2=isoout2SW(n,l)

                 isoout2IF(n,l)=(1.0-isowt)*isoout2IF(n,l)+isowt*oldiso3
	           oldiso3=isoout2IF(n,l)

                 isoout2GWfs(n,l)=(1.0-isowt)*isoout2GWfs(n,l)
     *                                                    +isowt*oldiso4
	           oldiso4=isoout2GWfs(n,l)

                 isoout2SWfs(n,l)=(1.0-isowt)*isoout2SWfs(n,l)
     *                                                    +isowt*oldiso5
	           oldiso5=isoout2SWfs(n,l)

                 isoout2IFfs(n,l)=(1.0-isowt)*isoout2IFfs(n,l)
     *                                                    +isowt*oldiso6
	           oldiso6=isoout2IFfs(n,l)


!                UPDATE THE TRACER MASSES IN CHANNEL - LAGGED BY 1 TIMESTEP
                 isostore2GW(n,l)=isostore1GW(n,l)+(isoin1GW(n,l)+
     *            isoin2GW(n,l)-isoout1GW(n,l)-isoout2GW(n,l))/2.
     	           isostore2GW(n,l)=amax1(isostore2GW(n,l),0.00001)

                 isostore2SW(n,l)=isostore1SW(n,l)+(isoin1SW(n,l)+
     *            isoin2SW(n,l)-isoout1SW(n,l)-isoout2SW(n,l))/2.
	           isostore2SW(n,l)=amax1(isostore2SW(n,l),0.00001)

                 isostore2IF(n,l)=isostore1IF(n,l)+(isoin1IF(n,l)+
     *            isoin2IF(n,l)-isoout1IF(n,l)-isoout2IF(n,l))/2.
	           isostore2IF(n,l)=amax1(isostore2IF(n,l),0.00001)

                 isostore2GWfs(n,l)=isostore1GWfs(n,l)+(isoin1GWfs(n,l)+
     *            isoin2GWfs(n,l)-isoout1GWfs(n,l)-isoout2GWfs(n,l))/2.
	           isostore2GWfs(n,l)=amax1(isostore2GWfs(n,l),0.00001)

                 isostore2SWfs(n,l)=isostore1SWfs(n,l)+(isoin1SWfs(n,l)+
     *            isoin2SWfs(n,l)-isoout1SWfs(n,l)-isoout2SWfs(n,l))/2.
	           isostore2SWfs(n,l)=amax1(isostore2SWfs(n,l),0.00001)

                 isostore2IFfs(n,l)=isostore1IFfs(n,l)+(isoin1IFfs(n,l)+
     *            isoin2IFfs(n,l)-isoout1IFfs(n,l)-isoout2IFfs(n,l))/2.
	           isostore2IFfs(n,l)=amax1(isostore2IFfs(n,l),0.00001)

!      if(n.eq.67) pause 'routing Tracer 5 through channel'


!                CHECK FOR CONVERGENCE - WE NEED THIS TO ITERATE ON ISOOUT2()
                 GWchk=abs(isoout2GW(n,l)-isooutold1)
                 SWchk=abs(isoout2SW(n,l)-isooutold2)
                 IFchk=abs(isoout2IF(n,l)-isooutold3)
                 GWfschk=abs(isoout2GWfs(n,l)-isooutold4)
                 SWfschk=abs(isoout2SWfs(n,l)-isooutold5)
                 IFfschk=abs(isoout2IFfs(n,l)-isooutold6)
                 if(GWchk.lt.0.00001.and.
     *              SWchk.lt.0.00001.and.
     *              IFchk.lt.0.00001.and.
     *              GWfschk.lt.0.00001.and.
     *              SWfschk.lt.0.00001.and.
     *              IFfschk.lt.0.00001) GOTO 50
                 isooutold1=isoout2GW(n,l)
                 isooutold2=isoout2SW(n,l)
                 isooutold3=isoout2IF(n,l)
                 isooutold4=isoout2GWfs(n,l)
                 isooutold5=isoout2SWfs(n,l)
                 isooutold6=isoout2IFfs(n,l)
               end do
   50          CONTINUE

! ** ERROR CHECK: TAS Jan 11/06 -- make new s/r for this
!      if(isoconcSW(n,l).gt.1)then
!	  print*,'isoconcSW',n,l
!        print*,isostore2SW(n,l),isostore1SW(n,l),store2(n)
!        print*,isoin1SW(n,l),isoin2SW(n,l)
!        print*,isoout1SW(n,l),isoout2SW(n,l)
!	  print*,isoconcSW(n,l)
!	  pause
!      endif
!	if(isoconcSWfs(n,l).gt.1)then
!	  print*,'isoconcSWfs',n,l
!        print*,isostore2SWfs(n,l),isostore1SWfs(n,l),store2(n)
!        print*,isoin1SWfs(n,l),isoin2SWfs(n,l)
!        print*,isoout1SWfs(n,l),isoout2SWfs(n,l)
!	  print*,isoconcSWfs(n,l)
!	  pause
!      endif
!      if(isoconcIF(n,l).gt.1)then
!	  print*,'isoconcIF',n,l
!        print*,isostore2IF(n,l),isostore1IF(n,l),store2(n)
!        print*,isoin1IF(n,l),isoin2IF(n,l)
!        print*,isoout1IF(n,l),isoout2IF(n,l)
!	  print*,isoconcIF(n,l)
!	  pause
!      endif
!      if(isoconcIFfs(n,l).gt.1)then
!	  print*,'isoconcIFfs',n,l
!        print*,isostore2IFfs(n,l),isostore1IFfs(n,l),store2(n)
!        print*,isoin1IFfs(n,l),isoin2IFfs(n,l)
!        print*,isoout1IFfs(n,l),isoout2IFfs(n,l)
!	  print*,isoconcIFfs(n,l)
!	  pause
!      endif
!      if(isoconcGW(n,l).gt.1)then
!	  print*,'isoconcGW',n,l
!        print*,isostore2GW(n,l),isostore1GW(n,l),store2(n)
!        print*,isoin1GW(n,l),isoin2GW(n,l)
!        print*,isoout1GW(n,l),isoout2GW(n,l)
!	  print*,isoconcGW(n,l)
!	  pause
!      endif
!      if(isoconcGWfs(n,l).gt.1)then
!	  print*,'isoconcGWfs',n,l
!        print*,isostore2GWfs(n,l),isostore1GWfs(n,l),store2(n)
!        print*,isoin1GWfs(n,l),isoin2GWfs(n,l)
!        print*,isoout1GWfs(n,l),isoout2GWfs(n,l)
!	  print*,isoconcGWfs(n,l)
!	  pause
!	endif

            else
              isostore2SW(n,l)=0.0
              isoout2SW(n,l)=0.0
              isostore2SWfs(n,l)=0.0
              isoout2SWfs(n,l)=0.0
              isostore2IF(n,l)=0.0
              isoout2IF(n,l)=0.0
              isostore2IFfs(n,l)=0.0
              isoout2IFfs(n,l)=0.0
              isostore2GW(n,l)=0.0
              isoout2GW(n,l)=0.0
              isostore2GWfs(n,l)=0.0
              isoout2GWfs(n,l)=0.0

!      if(n.eq.67) pause 'skipped routing Tracer 5'

            end if


!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!           ISOTOPES ARE CONSERVATIVE!
!           NOTE: CAN USE THIS TO MODEL GW RETARDATION OF FLOW
            decay = 0.0
	      if(lnxt.ne.0)then
!             UPDATE MASS BEING INPUT INTO NEXT GRID CELL IF THE NEXT
!             GRID IS NOT A "DUMMY GRID"
              isoin2GW(lll,lnxt)=(1-(decay/100))*(isoin2GW(lll,lnxt)+
     *                                        isoout2GW(n,l))
              isoin2SW(lll,lnxt)=(1-(decay/100))*(isoin2SW(lll,lnxt)+
     *                                       isoout2SW(n,l))
              isoin2IF(lll,lnxt)=(1-(decay/100))*(isoin2IF(lll,lnxt)+
     *                                        isoout2IF(n,l))
              isoin2GWfs(lll,lnxt)=(1-(decay/100))*(isoin2GWfs(lll,lnxt)
     *                                      +isoout2GWfs(n,l))
              isoin2SWfs(lll,lnxt)=(1-(decay/100))*(isoin2SWfs(lll,lnxt)
     *                                       +isoout2SWfs(n,l))
              isoin2IFfs(lll,lnxt)=(1-(decay/100))*(isoin2IFfs(lll,lnxt)
     *                                       +isoout2IFfs(n,l))

	      elseif(lnxt.ne.0 .and. res(n).gt.0) then
!             SET=0 AT LAKE OUTLETS B/C NO ROUTING THROUGH LAKES
              isostore2GW(n,l)=0.0
	        isoout2GW(n,l)=0.0
	        isoin2GW(lll,lnxt)=0.0
              isostore2GW(lll,lnxt)=0.0
              isostore2SW(n,l)=0.0
	        isoout2SW(n,l)=0.0
	        isoin2SW(lll,lnxt)=0.0
              isostore2SW(lll,lnxt)=0.0
              isostore2IF(n,l)=0.0
	        isoout2IF(n,l)=0.0
	        isoin2IF(lll,lnxt)=0.0
              isostore2IF(lll,lnxt)=0.0
              isostore2GWfs(n,l)=0.0
	        isoout2GWfs(n,l)=0.0
	        isoin2GWfs(lll,lnxt)=0.0
              isostore2GWfs(lll,lnxt)=0.0
              isostore2SWfs(n,l)=0.0
	        isoout2SWfs(n,l)=0.0
	        isoin2SWfs(lll,lnxt)=0.0
              isostore2SWfs(lll,lnxt)=0.0
              isostore2IFfs(n,l)=0.0
	        isoout2IFfs(n,l)=0.0
	        isoin2IFfs(lll,lnxt)=0.0
              isostore2IFfs(lll,lnxt)=0.0
            endif


!           ADD MASS OUTFLOW OF CURRENT GRID AS INFLOW TO LAKE (NEXT GRID) (cms)
!           BUT ONLY IF ON THE EDGE OF A LAKE (NOT IN A LAKE)
c            if(ireach(lll).gt.0.and. iz.ne.jz.and.dds_flag.ne.1) then
            if(ireach(lll).gt.0.and. iz.ne.jz) then

              isolakeGW(rbin)=isolakeGW(rbin)+isoout2GW(n,l)/t      


!	write(777,'(4i10,999(f12.3))')n,l,lll,rbin,isolakeGW(rbin,jjz),
!     *            isoout2GW(n,l)/t,isoin2GW(lll,lnxt)/t

	      endif


!            ERROR CHECK: ALL COMPONENTS MUST SUM TO TOTAL SIM'D Q
!            CAUSE SOME CONC'NS > 1.0 FROM ABOVE...
             isosumQ(n)=isoout2SW(n,l)+isoout2SWfs(n,l)+
     *                  isoout2IF(n,l)+isoout2IFfs(n,l)+
     *                  isoout2GW(n,l)+isoout2GWfs(n,l)
!	print*,isosumQ(n)


!           ACCUMULATE MASSES FOR MASS BALANCE CHECK AT END OF RUN:
            massin(n)=(isoin1GW(n,l)+isoin2GW(n,l))/2.
     *               +(isoin1GWfs(n,l)+isoin2GWfs(n,l))/2.
     *               +(isoin1SW(n,l)+isoin2SW(n,l))/2.
     *               +(isoin1SWfs(n,l)+isoin2SWfs(n,l))/2.
     *               +(isoin1IF(n,l)+isoin2IF(n,l))/2.
     *               +(isoin1IFfs(n,l)+isoin2IFfs(n,l))/2.
            massout(n)=(isoout1GW(n,l)+isoout2GW(n,l))/2.
     *                +(isoout1GWfs(n,l)+isoout2GWfs(n,l))/2.
     *                +(isoout1SW(n,l)+isoout2SW(n,l))/2.
     *                +(isoout1SWfs(n,l)+isoout2SWfs(n,l))/2.
     *                +(isoout1IF(n,l)+isoout2IF(n,l))/2.
     *                +(isoout1IFfs(n,l)+isoout2IFfs(n,l))/2.
            masstore(n)=(isostore2GW(n,l)-isostore1GW(n,l))
     *                 +(isostore2GWfs(n,l)-isostore1GWfs(n,l))
     *                 +(isostore2SW(n,l)-isostore1SW(n,l))
     *                 +(isostore2SWfs(n,l)-isostore1SWfs(n,l))
     *                 +(isostore2IF(n,l)-isostore1IF(n,l))
     *                 +(isostore2IFfs(n,l)-isostore1IFfs(n,l))
            ISOdelta(n)=massin(n)-massout(n)

 !          CHECK OVERALL MASS BALANCE:
 !          MASS IN - MASS OUT = MASS STORED (IF EVERYTHING WORKS!) 
            if(iz.ne.jz.and.n.eq.nnprint)then
              sqerr=(ISOdelta(n)-masstore(n))**2.
              write(91,8000)time,massin(n),massout(n),ISOdelta(n),
     *                      masstore(n),sqerr
	      endif

          endif  ! SLOPE>=0 AND L.NE.0

          else   ! IREACH(N).LE.0
!           IN A LAKE, DON'T RUN TRACER!
	      CONTINUE 
	    endif


        end do   ! GRID NO. LOOP

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!       WRITE OUT DATA TO TRACER.CSV FILE
!        if(jan.eq.3)then    ! Prints every 0.25 hours
        if(iz.ne.jz)then     ! Prints every hour
          do l=1,no
!           STORES THE GRID NO. AS FXN OF THE GAUGE(BSN) NO.
            nn(l)=s(iy(l),jx(l))
	    end do

!         WRITE TO THE FILE ONLY WHEN YOU ARE AT THE FLAGGED OUTPUT
!         GRID=NNPRINT. 
!         UNITS: isoout2IBN(n,ll)=[kg]/[s]=[kg/s]=[m^3/s]
!         NB: TRACER OUTFLOW BASED ON MASS TO ACCOUNT FOR EVAP. LOSS

c        if(no.gt.30.and.dds_flag.ne.1) then
        if(no.gt.30) then
          write(90,99887)time,(isoout2GW(nn(l),l)/t,l=1,no)

c        elseif(dds_flag.ne.1)then
          write(90,9005)time,
     *     (isoout2SW(nn(l),l)/t,!isoconcSW(nn(l),l),
     *      isoout2SWfs(nn(l),l)/t,!isoconcSWfs(nn(l),l),
     *      isoout2IF(nn(l),l)/t,!isoconcIF(nn(l),l),
     *      isoout2IFfs(nn(l),l)/t,!isoconcIFfs(nn(l),l),
     *      isoout2GW(nn(l),l)/t,!isoconcGW(nn(l),l),
     *      isoout2GWfs(nn(l),l)/t,l=1,no)!isoconcGWfs(nn(l),l),

	write(999,12345)time,(qo2(nn(l)),isosumQ(nn(l))/t,
     *                (qo2(nn(l))-isosumQ(nn(l))/t),l=1,no)

12345 FORMAT(f10.3,<no>(3(',',f15.5)))

        endif

c        if(wetflg.eq.'y'.and.dds_flag.ne.1) write(93,9101)time,
        if(wetland_flag(n)) write(93,9101)time,
     *    (l,nn(l),qlz(nn(l)),qswevp(nn(l)),isoout2wet(nn(l),l)/t,
     *     isoconcwet(nn(l),l),l=1,no)
        endif
!       filename(90)='..\simout\tracer.csv'
!       filename(91)='..\simout\tracerMB.csv'     ! added Oct.30/03 TS
!       filename(93)='..\simout\tracerWET.csv'    ! added Dec.01/03 TS
!       filename(94)='..\simout\tracerWETMB.csv'  ! added Dec.01/03 TS


	endif      ! SNOWMELT TRACER



! FORMATS:
 8000 FORMAT(f10.3,256(',',f15.5))
 9005 FORMAT(f10.3,<no>(12(',',f15.5)))
 9000 FORMAT(f10.3,<no>(2(',',I15),4(',',f15.5)))
 9101 FORMAT(f10.3,<no>(2(',',I15),4(',',f15.5)))
20400 FORMAT(f10.3,999(',',f15.5))
99887 FORMAT(f10.3,<no>((',',f15.5)))
99888 FORMAT(f10.3,<no>(2(',',f15.5)))


      RETURN
      END SUBROUTINE tracer5