      SUBROUTINE tracer3(iz,jz,time,t,jan,tdum)

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
! RAIN-ON-STREAM TRACER
! This subroutine is designed to route isotope tracers through the
! WATFLOOD model to track the amount of water coming from new rain.  
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

c      USE area1
c	USE area2
c      USE area3
c      USE area4
c      USE area5
c      USE area6
cc	USE area9
c	USE area10
cc      USE areaet
cc	USE areamelt
c      USE areatrc
c	USE areawet

      USE area_watflood
      implicit none



!*****NOTE: THIS PROGRAM ASSUMES (FOR NOW) ZERO BACKGROUND CONCENTRATIONS
!*****OF THE TRACERS.  FROM RUN-TO-RUN THOUGH, TRACER MASSES ARE CONSERVED.
!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      INTEGER :: n,lll,l,iz,jz,ii,jj,i,j,jan,ios,inxt,jnxt,lnxt
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
!     TRACER 3 => RAIN-ON-STREAM TRACER (QSTREAM)               * *
!                                                               * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      if(itrace.eq.3)then

!       INITIALIZE AND RESET VALUES FOR EACH TIME STEP
        do n=1,na
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)	    
	    if(l.ne.0)then
!           INITIALIZE TRACERS:
            isoin1P(n,l)=isoin2P(n,l)
            isoin2P(n,l)=0.0
            isoout1P(n,l)=isoout2P(n,l)
            isoout2P(n,l)=0.0
            isostore1P(n,l)=isostore2P(n,l)
            isostore2P(n,l)=0.0
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

          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)

!         RUN TRACER CODE ONLY IF NOT IN A LAKE
          if(ireach(n).le.0.0 .OR.
     *       ireach(n).gt.0.and.res(n).gt.0) then

!         WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
!         IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!         IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
          if(slope(n).gt.0.0.and.l.ne.0)then
            if(store2(n).ne.0.0)then
               isooutold=1.0e+25

!              SET INITIAL CONCENTRATION AS 1 kg/m^3
!              INCLUDES ONLY GW
               isoin2P(n,l)=isoin2P(n,l)+1.0*qstrm(n)*t

!              CONVERGENCE LOOP:
               do ijk=1,50

!                UPDATE THE TRACER MASS IN CHANNEL n
                 isostore2P(n,l)=isostore1P(n,l)+(isoin1P(n,l)+
     *                      isoin2P(n,l)-isoout1P(n,l)-isoout2P(n,l))/2.
	           isostore2P(n,l)=amax1(isostore2P(n,l),0.00001)

!                UPDATE THE CONCENTRATION
                 isoconcP(n,l)=isostore2P(n,l)/store2(n)

!                COMPUTE THE TRACER MASS LEAVING THE GRID...ACCOUNTING FOR
!                EVAP LOSS TO PRESERVE FLOW BALANCES.
!                THIS SHOULD BE BACK DATED ONE TIME STEP
                 isoout2P(n,l)=isoconcP(n,l)*qo2(n)*t
     *                      -isoconcP(n,l)*strloss(n)*t
	           isoout2P(n,l)=amax1(isoout2P(n,l),0.00001)

!                CHECK FOR CONVERGENCE - WE NEED THIS TO ITERATE ON ISOOUT2()
                 if(abs(isoout2P(n,l)-isooutold).lt.0.0001) 
     *              GOTO 30
                 isooutold=isoout2P(n,l)
               end do
   30          CONTINUE

            else
              isostore2P(n,l)=0.0
              isoout2P(n,l)=0.0
            end if

!            DECREASED BY THE AMOUNT OF DEPOSITION:            
!            wi2(lll)=(1-(sdep/100))*(wi2(lll)+wo2(n))

!           ISOTOPES ARE CONSERVATIVE!
            decay=0.0
!           UPDATE MASS BEING INPUT INTO NEXT GRID CELL IF THE NEXT
!           GRID IS NOT A "DUMMY GRID"
            if(lnxt.ne.0)
     *         isoin2P(lll,lnxt)=(1-(decay/100))*(isoin2P(lll,lnxt)+
     *                                      isoout2P(n,l))

!           ACCUMULATE MASSES FOR MASS BALANCE CHECK AT END OF RUN:
            massin(n)=(isoin1P(n,l)+isoin2P(n,l))/2.
            massout(n)=(isoout1P(n,l)+isoout2P(n,l))/2.
            masstore(n)=(isostore2P(n,l)-isostore1P(n,l))
            ISOdelta(n)=massin(n)-massout(n)


 !          CHECK OVERALL MASS BALANCE:
 !          MASS IN - MASS OUT = MASS STORED (IF EVERYTHING WORKS!) 
            if(iz.ne.jz.and.n.eq.nnprint)then
              sqerr=(ISOdelta(n)-masstore(n))**2.
              write(91,8000)time,massin(n),massout(n),ISOdelta(n),
     *                      masstore(n),sqerr
	      endif 

          endif  ! SLOPE>=0 AND L.NE.0

          else
!           IN A LAKE, DON'T RUN TRACER!
	      CONTINUE 
	    endif


        end do   ! GRID NO. LOOP

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
          write(90,9003)time,
     *     (l,nn(l),qo2(nn(l)),qstrm(nn(l)),strloss(nn(l)),
     *      isoout2P(nn(l),l)/t,isoconcP(nn(l),l),l=1,no)
        endif

	endif      ! RAIN-ON-STREAM TRACER



! FORMATS:
 8000 FORMAT(f10.3,256(',',f15.5))
 9000 FORMAT(f10.3,<no>(2(',',I15),4(',',f15.5)))
 9003 FORMAT(f10.3,<no>(2(',',I15),5(',',f15.5)))


      RETURN
      END SUBROUTINE tracer3
