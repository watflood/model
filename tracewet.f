      SUBROUTINE tracewet(iz,jz,time,n,l,t)

!***********************************************************************
!    Copyright (C) 2016 by Tricia Stadnyk and Tegan Holmes  
        
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
! S/R TRACERWET
! This subroutine is designed to route isotope tracers through the WETLAND
! portion of the WATFLOOD model.  It will do so based on a mixed cell model
! (continuity), assuming instant mixing.
! Adapted from LFL's sediment transport subroutine.
!
! Created:  November 2003 - TS
!
! REV. 9.1.49  Nov.  23/03 - TS: Added wetlands to GW Tracer + Wetland Tracer
! REV.         Feb.  4/16  - TH: Replaced convergence loop with direct solution, matched evap back with route.f
!
!*****************************************************************************

   

      use area_watflood
      implicit none	
	
	REAL    :: t,time
	REAL*4  :: GWchk,GWfschk,SWchk,SWfschk,IFchk,IFfschk
	INTEGER :: iwet,n,l,lll,iz,jz


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	if(itrace.eq.1)then      ! TRACER FOR GLACIERS
!       3D ARRAYS

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	elseif(itrace.eq.2)then      ! TRACER FOR LCOVER
!       3D ARRAYS



! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	elseif(itrace.eq.4)then      ! TRACER W/ FLOWTYPES

!       UPDATE THE TRACER MASS IN WETLAND
        isowstore2(n,l)=(isowstore1(n,l)+(isoin1wet(n,l)+
     *               isoin2wet(n,l)-isoout1wet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
	  isowstore2(n,l)=amax1(isowstore2(n,l),0.00001)

        isowstore2SW(n,l)=(isowstore1SW(n,l)+(isoin1SWwet(n,l)+
     *          isoin2SWwet(n,l)-isoout1SWwet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
	  isowstore2SW(n,l)=amax1(isowstore2SW(n,l),0.00001)

        isowstore2IF(n,l)=(isowstore1IF(n,l)+(isoin1IFwet(n,l)+
     *          isoin2IFwet(n,l)-isoout1IFwet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
	  isowstore2IF(n,l)=amax1(isowstore2IF(n,l),0.00001)

!         UPDATE THE CONCENTRATION IN THE WETLAND
        isoconcwet(n,l)=isowstore2(n,l)/wstore2(n)
        isoconcSWwet(n,l)=isowstore2SW(n,l)/wstore2(n)
        isoconcIFwet(n,l)=isowstore2IF(n,l)/wstore2(n)

!         COMPUTE THE TRACER MASS LEAVING THE GRID...
!         EVAP LOSS TO PRESERVE FLOW BALANCES.
!         THIS SHOULD BE BACK DATED ONE TIME STEP
        isoout2wet(n,l)=isoconcwet(n,l)*qowet2(n)*t
     *                   +isoconcwet(n,l)*qswevp(n)*t
     	  isoout2wet(n,l)=amax1(isoout2wet(n,l),0.00001)

        isoout2SWwet(n,l)=isoconcSWwet(n,l)*qowet2(n)*t
     *                   +isoconcSWwet(n,l)*qswevp(n)*t
        isoout2SWwet(n,l)=amax1(isoout2SWwet(n,l),0.00001)

        isoout2IFwet(n,l)=isoconcIFwet(n,l)*qowet2(n)*t
     *                   +isoconcIFwet(n,l)*qswevp(n)*t
        isoout2IFwet(n,l)=amax1(isoout2IFwet(n,l),0.00001)

        isooutwet1=isoconcwet(n,l)*qowet2(n)*t
        isooutwet2=isoconcSWwet(n,l)*qowet2(n)*t
        isooutwet3=isoconcIFwet(n,l)*qowet2(n)*t

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	elseif(itrace.eq.5)then      ! TRACER W/ FLOWTYPES+MELT
!       2D ARRAYS
        
!       UPDATE THE TRACER MASS IN WETLAND
        isowstore2(n,l)=(isowstore1(n,l)+(isoin1wet(n,l)+
     *               isoin2wet(n,l)-isoout1wet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
!	  isowstore2(n,l)=amax1(isowstore2(n,l),0.00001)

        isowstore2SW(n,l)=(isowstore1SW(n,l)+(isoin1SWwet(n,l)+
     *          isoin2SWwet(n,l)-isoout1SWwet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
!	  isowstore2SW(n,l)=amax1(isowstore2SW(n,l),0.00001)

        isowstore2IF(n,l)=(isowstore1IF(n,l)+(isoin1IFwet(n,l)+
     *          isoin2IFwet(n,l)-isoout1IFwet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
!	  isowstore2IF(n,l)=amax1(isowstore2IF(n,l),0.00001)

        isowstore2fs(n,l)=(isowstore1fs(n,l)+(isoin1fswet(n,l)+
     *          isoin2fswet(n,l)-isoout1fswet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
!	  isowstore2fs(n,l)=amax1(isowstore2fs(n,l),0.00001)

        isowstore2SWfs(n,l)=(isowstore1SWfs(n,l)+(isoin1SWfswet(n,l)+
     *    isoin2SWfswet(n,l)-isoout1SWfswet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
!	  isowstore2SWfs(n,l)=amax1(isowstore2SWfs(n,l),0.00001)

        isowstore2IFfs(n,l)=(isowstore1IFfs(n,l)+(isoin1IFfswet(n,l)+
     *    isoin2IFfswet(n,l)-isoout1IFfswet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
!	  isowstore2IFfs(n,l)=amax1(isowstore2IFfs(n,l),0.00001)


!         UPDATE THE CONCENTRATION IN THE WETLAND
          isoconcwet(n,l)=isowstore2(n,l)/wstore2(n)
          isoconcSWwet(n,l)=isowstore2SW(n,l)/wstore2(n)
          isoconcIFwet(n,l)=isowstore2IF(n,l)/wstore2(n)
          isoconcfswet(n,l)=isowstore2fs(n,l)/wstore2(n)
          isoconcSWfswet(n,l)=isowstore2SWfs(n,l)/wstore2(n)
          isoconcIFfswet(n,l)=isowstore2IFfs(n,l)/wstore2(n)

!         COMPUTE THE TRACER MASS LEAVING THE GRID...
!         EVAP LOSS TO PRESERVE FLOW BALANCES.
!         THIS SHOULD BE BACK DATED ONE TIME STEP
          isoout2wet(n,l)=isoconcwet(n,l)*qowet2(n)*t
     *                   +isoconcwet(n,l)*qswevp(n)*t
!          isoout2wet(n,l)=amax1(isoout2wet(n,l),0.00001)

          isoout2SWwet(n,l)=isoconcSWwet(n,l)*qowet2(n)*t
     *                   +isoconcSWwet(n,l)*qswevp(n)*t
!     	    isoout2SWwet(n,l)=amax1(isoout2SWwet(n,l),0.00001)

          isoout2IFwet(n,l)=isoconcIFwet(n,l)*qowet2(n)*t
     *                   +isoconcIFwet(n,l)*qswevp(n)*t
!     	    isoout2IFwet(n,l)=amax1(isoout2IFwet(n,l),0.00001)

          isoout2fswet(n,l)=isoconcfswet(n,l)*qowet2(n)*t
     *                   +isoconcfswet(n,l)*qswevp(n)*t
!          isoout2fswet(n,l)=amax1(isoout2fswet(n,l),0.00001)

          isoout2SWfswet(n,l)=isoconcSWfswet(n,l)*qowet2(n)*t
     *                   +isoconcSWfswet(n,l)*qswevp(n)*t
!     	    isoout2SWfswet(n,l)=amax1(isoout2SWfswet(n,l),0.00001)

          isoout2IFfswet(n,l)=isoconcIFfswet(n,l)*qowet2(n)*t
     *                   +isoconcIFfswet(n,l)*qswevp(n)*t
!     	    isoout2IFfswet(n,l)=amax1(isoout2IFfswet(n,l),0.00001)

      
        isooutwet1=isoconcwet(n,l)*qowet2(n)*t
        isooutwet2=isoconcSWwet(n,l)*qowet2(n)*t
        isooutwet3=isoconcIFwet(n,l)*qowet2(n)*t
        isooutwet4=isoconcfswet(n,l)*qowet2(n)*t
        isooutwet5=isoconcSWfswet(n,l)*qowet2(n)*t
        isooutwet6=isoconcIFfswet(n,l)*qowet2(n)*t

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	elseif(itrace.eq.0. OR. itrace.eq.3 .OR. itrace.eq.100)then    
                                   ! TRACER FOR GW, RAIN, AND SUB-BASINS

!       UPDATE THE TRACER MASS IN WETLAND
        isowstore2(n,l)=(isowstore1(n,l)+(isoin1wet(n,l)+
     *              isoin2wet(n,l)-isoout1wet(n,l))/2.)
     *             /(1.0+(qowet2(n)+qswevp(n))*t
     *              /2.0/wstore2(n))
	  isowstore2(n,l)=amax1(isowstore2(n,l),0.00001)

!         UPDATE THE CONCENTRATION IN THE WETLAND
        isoconcwet(n,l)=isowstore2(n,l)/wstore2(n)

!         COMPUTE THE TRACER MASS LEAVING THE GRID...
!         EVAP LOSS TO PRESERVE FLOW BALANCES.
!         THIS SHOULD BE BACK DATED ONE TIME STEP
        isoout2wet(n,l)=isoconcwet(n,l)*qowet2(n)*t
     *                   -isoconcwet(n,l)*qswevp(n)*t
        isoout2wet(n,l)=amax1(isoout2wet(n,l),0.00001)

        isooutwet=isoconcwet(n,l)*qowet2(n)*t

      endif           ! GW TRACER
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *



! FORMATS:
 8000 FORMAT(g10.3,256(',',g15.5))    ! nk  changed Jan 5/11
c 8000 FORMAT(f10.3,256(',',f15.5))
10000 FORMAT('   time   ,',
     *<no>('         basin#,          grid#,',
     *'        qlz    ,       qswevp ,        ISOoutQ,',
     *'        ISOconc,'))
11000 FORMAT('    time  ,        IN     ,       OUT     ,',
     *'        DELTA  ,       STORE   ,       SSE     ,')



      RETURN
	END SUBROUTINE tracewet