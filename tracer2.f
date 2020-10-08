      SUBROUTINE tracer2(iz,jz,time,t,jan,tdum)

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
! LANDCOVER TRACER
! This subroutine is designed to route isotope tracers through the
! WATFLOOD model to track the amount of runoff coming from individual landclasses.  
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
c	USE area10
cc      USE areaet
c	USE areamelt
c      USE areatrc
c	USE areawet

      use area_watflood
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
!     TRACER 2 => LANDCOVER TRACER                              * *
!                                                               * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      if(itrace.eq.2)then    

!       INITIALIZE AND RESET VALUES FOR EACH TIME STEP
        do n=1,na
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)	    
	    if(l.ne.0)then
            isosum(n,l)=0.0
		  do ii=1,classcount
!             INITIALIZE TRACERS:
              isoin1LC(n,l,ii)=isoin2LC(n,l,ii)
              isoin2LC(n,l,ii)=0.0
              isoout1LC(n,l,ii)=isoout2LC(n,l,ii)
              isoout2LC(n,l,ii)=0.0
              isostore1LC(n,l,ii)=isostore2LC(n,l,ii)
              isostore2LC(n,l,ii)=0
	        isoin1BLZS(n,l,ii)=isoin2BLZS(n,l,ii)
              isoin2BLZS(n,l,ii)=0.0
	        isoout1BLZS(n,l,ii)=isoout2BLZS(n,l,ii)
	        isoout2BLZS(n,l,ii)=0.0
              isoBLZS1(n,l,ii)=isoBLZS2(n,l,ii)
	        isoBLZS2(n,l,ii)=0.0
	      end do
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


          do ii=1,classcount

!       ^^^^^^^^^^^^^ GW ROUTE THE TRACER ^^^^^^^^^^^^^^^^^
!           TO GET MASS OF TRACER ENTERING STREAM FROM QLZ
            if(qlz(n).ne.0.0.and.l.ne.0)then
              isooutold=1.0e+25

!             CALC MASS OF TRACER GOING INTO GW STORAGE, SET INITIAL
!             CONCENTRATION AS 1 kg/m^3 * m^s/3 * s = KG
              isoin2BLZS(n,l,ii)=isoin2BLZS(n,l,ii)+
     *                       1.0*drng(n,ii)*frac(n)*tdum*(1-sca(n,ii))*t
     *                      +1.0*drngfs(n,ii)*frac(n)*tdum*sca(n,ii)*t

!             CONVERGENCE LOOP:
              do ijk=1,50
!               UPDATE THE TRACER MASS IN GW STORAGE
                isoBLZS2(n,l,ii)=isoBLZS1(n,l,ii)+
     *            (isoin1BLZS(n,l,ii)+isoin2BLZS(n,l,ii)-
     *             isoout1BLZS(n,l,ii)-isoout2BLZS(n,l,ii))/2.
	          isoBLZS2(n,l,ii)=amax1(isoBLZS2(n,l,ii),0.00001)

!               UPDATE CONCENTRATION USING INSTANT MIXING!! 
!               GIVES UNITS KG/m^3
                isoconcBLZS(n,l,ii)=isoBLZS2(n,l,ii)/lzs(n)

!               COMPUTE THE TRACER MASS LEAVING THE LZS
!               THIS SHOULD BE BACK DATED ONE TIME STEP
!               KG/m^3 * m^3/s * s = KG
!               ASSUME: NO EVAP IN SOIL COLUMN!!
                isoout2BLZS(n,l,ii)=isoconcBLZS(n,l,ii)*qlz(n)*t

!               CHECK FOR CONVERGENCE - WE NEED THIS TO ITERATE ON ISOOUT2()
                if(abs(isoout2BLZS(n,l,ii)-isooutold).lt.0.01*isooutold) 
     *            GOTO 21
                isooutold=isoout2BLZS(n,l,ii)
              end do
   21         CONTINUE
            else
             isoBLZS2(n,l,ii)=0.0
             isoout2BLZS(n,l,ii)=0.0
            end if

!       ^^^^^^^^^^^^^ CHANNEL ROUTE THE TRACER ^^^^^^^^^^^^^^^^^
!           ROUTED QLZ CONTRIBUTION ADDED TO CHANNEL AND ROUTED AGAIN.
!           WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
!           IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!           IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
            if(slope(n).gt.0.0.and.l.ne.0)then
              if(store2(n).ne.0.0)then
                isooutold=1.0e+25
!               SET INITIAL CONCENTRATION AS 1 kg/m^3
!        * * * NOTE: THIS DOES NOT ACCOUNT FOR WETLANDS!!!
!               ADD TRACER TO SW FLOWS, PLUS ROUTED GW TRACER:
                isoin2LC(n,l,ii)=isoin2LC(n,l,ii)+1.0*(q1(n,ii)+
     *                           q1fs(n,ii)+qint(n,ii)+qintfs(n,ii))*t+
     *                           isoout2BLZS(n,l,ii)

!               CONVERGENCE LOOP:
                do ijk=1,50
!                 UPDATE THE TRACER MASS IN CHANNEL n
                  isostore2LC(n,l,ii)=isostore1LC(n,l,ii)+
     *              (isoin1LC(n,l,ii)+isoin2LC(n,l,ii)-
     *               isoout1LC(n,l,ii)-isoout2LC(n,l,ii))/2.
	            isostore2LC(n,l,ii)=amax1(isostore2LC(n,l,ii),0.00001)

!                 UPDATE THE CONCENTRATION
                  isoconcLC(n,l,ii)=isostore2LC(n,l,ii)/store2(n)

!                 COMPUTE THE TRACER MASS LEAVING THE GRID...
!                 EVAP LOSS TO PRESERVE FLOW BALANCES.
!                 THIS SHOULD BE BACK DATED ONE TIME STEP
                  isoout2LC(n,l,ii)=isoconcLC(n,l,ii)*qo2(n)*t
     *                             -isoconcLC(n,l,ii)*strloss(n)*t
	            isoout2LC(n,l,ii)=amax1(isoout2LC(n,l,ii),0.00001)

!                 CHECK FOR CONVERGENCE - WE NEED THIS TO ITERATE ON ISOOUT2()
                  if(abs(isoout2LC(n,l,ii)-isooutold).lt.0.01*isooutold) 
     *              GOTO 20
                  isooutold=isoout2LC(n,l,ii)
                end do
   20           CONTINUE

              else
               isostore2LC(n,l,ii)=0.0
               isoout2LC(n,l,ii)=0.0
              end if
!             ISOTOPES ARE CONSERVATIVE!
              decay = 0.0
!             UPDATE MASS BEING INPUT INTO NEXT GRID CELL IF THE NEXT 
!             GRID IS NOT A "DUMMY GRID"
              if(lnxt.ne.0)
     *           isoin2LC(lll,lnxt,ii)=(1-(decay/100))*
     *                       (isoin2LC(lll,lnxt,ii)+isoout2LC(n,l,ii))

!             ACCUMULATE MASSES FOR MASS BALANCE CHECK AT END OF RUN:
              massin(n)=(isoin1LC(n,l,ii)+isoin2LC(n,l,ii))/2.
              massout(n)=(isoout1LC(n,l,ii)+isoout2LC(n,l,ii))/2.
              masstore(n)=(isostore2LC(n,l,ii)-isostore1LC(n,l,ii))
              ISOdelta(n)=massin(n)-massout(n)

            endif    ! SLOPE>=0 AND L.NE.0
            isosum(n,l)=isosum(n,l)+isoout2LC(n,l,ii)


	    end do     ! classcount LOOP


!         CHECK OVERALL MASS BALANCE:
!         MASS IN - MASS OUT = MASS STORED (IF EVERYTHING WORKS!) 
          if(iz.ne.jz.and.n.eq.nnprint)then
            sqerr=(ISOdelta(n)-masstore(n))**2.
            write(91,8000)time,massin(n),massout(n),ISOdelta(n),
     *                    masstore(n),sqerr
	    endif

          else
!          IN A LAKE, DON'T RUN TRACER!
	     CONTINUE 
	    endif

        end do       ! GRID NO./BASIN NO. LOOP

!       WRITE OUT DATA TO TRACER.CSV FILE
        if(iz.ne.jz)then     ! Prints every hour
!        if(jan.eq.3)then    ! Prints every 0.25 hours
          do l=1,no
!           STORES THE GRID NUMBERS AS FXN OF GAUGE NO.
            nn(l)=s(iy(l),jx(l))
          end do
!         WRITE TO THE FILE ONLY WHEN YOU ARE AT THE FLAGGED OUTPUT
!         GRID=NNPRINT. USE ii TO CYCLE THRU LCOVERS CAUSE W/IN l DO-LOOP.
!         UNITS: isoout2LC(n,ii)=[kg]/[s]=[m^3/s]
!         NB: TRACER OUTFLOW BASED ON MASS TO ACCOUNT FOR EVAP. LOSS
          write(90,9002)time,(l,qo2(nn(l)),isosum(nn(l),l)/t,
     *     (isoout2LC(nn(l),l,ii)/t,isoconcLC(nn(l),l,ii),
     *                ii=1,classcount),l=1,no)

!     debug check added May 11/10  nk
	if(iopt.ge.1)then

!        added Mar. 17/06    nk  temporary fix to get glacier component.
!        outputs are grouped by subwatershed
         if(mod(iz,24).eq.0)then
           ii=6
!           write(92,9003)time,
!     *     ((qo2(nn(l))*(1.0-isoconcLC(nn(l),l,ii)),ii=1,classcount),l=1,no)
           write(92,9003)time,
     *     (qo2(nn(l))*(1.0-isoconcLC(nn(l),l,ii)),l=1,no)
         endif

      endif

	  endif

	endif      ! LANDCOVER TRACER



! FORMATS:
 8000 FORMAT(f10.3,256(',',f15.5))
 9000 FORMAT(f10.3,<no>(2(',',I15),4(',',f15.5)))
 9002 FORMAT(f10.3,<no>(2(',',I15),2(',',f15.5),
     *                <classcount>(2(',',f15.5))))
 9003 format(f10.0,<no*classcount>(',',f10.3))


      RETURN
      END SUBROUTINE tracer2