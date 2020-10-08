      SUBROUTINE tracer4(iz,jz,time,t,jan,tdum,jjz,jjzold)

!***********************************************************************
!    Copyright (C) 2008 by Tricia Stadnyk, Tegan Holmes and Nicholas Kouwen 
        
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
! FLOW-TYPE TRACER
! This subroutine is designed to route isotope tracers through the
! WATFLOOD model to track the proportion of runoff from each different flow type.  
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
!     REV. 9.8.18  Apr.  26/11 - NK: Added in-basin check in tracer4
!     REV. 9.8.77  July  2013  - TH: Added lake routing, diversions & nudging
!     REV.         Feb.  4/16  - TH: Replaced convergence loop with direct solution
!
!*****************************************************************************


      use area_watflood
	implicit none

!     SEE tracerGW.f FOR MORE COMPLETE COMMENTS
!     FROM RUN-TO-RUN THOUGH, TRACER MASSES ARE CONSERVED.
!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      INTEGER :: n,m,lll,l,iz,jz,i,j,jan,inxt,jnxt,lnxt
	INTEGER :: rbin,resnum,jjz,jjzold,ii,jj
      REAL*4  :: t,time,GWchk,SWchk,IFchk
     	REAL*4  :: tdum,roe 
     	logical :: firstpass,trace_off
     	
     	data firstpass/.true./
      data trace_off/.false./

      if(firstpass)then
          do l=1,no
!     rev. 9.7.16  Jan.  05/11  - NK: Fixed init flows outside sub-basin
              if(inbsnflg(l).eq.1)then    ! added nk Jan. 05/11
!             STORES THE GRID NO. AS FXN OF THE GAUGE(BSN) NO.
                  nn(l)=s(iy(l),jx(l))
              else
!     rev. 9.8.18  Apr.  26/11  - NK: Added in-basin check in tracer4
                  nn(l)=-1
              endif
          end do
          do l=1,no
d             print*,l,iy(l),jx(l),nn(l)    !ok
              if(nn(l).le.0)then
                  trace_off=.true.
                  trcflg='n'
                  print*
                  print*,'WARNING:'
                  print*,'Tracers do not work when lakes are outside '
                  print*,'the model domain as when using a '
                  print*,'sub-watershed derived from a larger domain '
                  print*,'without making new reach numbers'
                  print*,'Reach # ',l,' is outside the watershed'
                  print*,'Grid #  ',s(iy(l),jx(l))
                  print*,'inbsnflg',inbsnflg
                  print*,'NOTE:'
                  Print*,'Program continues without doing tacers<<<<<<<'
                  print*
              endif
c              pause 'In tracer4'
c              return
          end do
	endif

      res_set: if(index==0)then
        do n=1,naa
          res(n)=0
        end do

        do n=1,naa
          lll=next(n)
          do l=1,noresv
            if(yyy(n).eq.ires(l).and.xxx(n).eq.jres(l)) res(n)=l
          end do
        end do
      endif res_set
      

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                                                               * *
!     TRACER 4 => FLOW-TYPE TRACER (SW, IF, GW)                 * *
!                                                               * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     INITIALIZE AND RESET VALUES FOR EACH TIME STEP
      initialize: do n=1,naa
        i=yyy(n)
        j=xxx(n)
        
        
        
        l=nbasin(i,j)
c        l=1
        
        
	  lll=next(n)
	  resnum=ireach(n)
	  rbin=ireach(lll)

        if(rbin.ne.0) then
          isolakeGW(rbin)=0.0
          isolakeSW(rbin)=0.0
          isolakeIF(rbin)=0.0
          isolakest(rbin)=0.0
	  endif

	  isosumQ(n)=0.0

	  if(l.ne.0)then
          isoin1GW(n,l)=isoin2GW(n,l)
          isoin2GW(n,l)=0.0
          isoout1GW(n,l)=isoout2GW(n,l)
          isoout2GW(n,l)=0.0
          isostore1GW(n,l)=isostore2GW(n,l)
          isostore2GW(n,l)=0.0
          isoin1SW(n,l)=isoin2SW(n,l)
          isoin2SW(n,l)=0.0
          isoout1SW(n,l)=isoout2SW(n,l)
          isoout2SW(n,l)=0.0
          isostore1SW(n,l)=isostore2SW(n,l)
          isostore2SW(n,l)=0.0
          isoin1IF(n,l)=isoin2IF(n,l)
          isoin2IF(n,l)=0.0
          isoout1IF(n,l)=isoout2IF(n,l)
          isoout2IF(n,l)=0.0
          isostore1IF(n,l)=isostore2IF(n,l)
          isostore2IF(n,l)=0.0
	    if(wetflg.eq.'y')then
            isoin1wet(n,l)=isoin2wet(n,l)
            isoin2wet(n,l)=0.0
            isoout1wet(n,l)=isoout2wet(n,l)
            isoout2wet(n,l)=0.0
            isowstore1(n,l)=isowstore2(n,l)
            isowstore2(n,l)=0.0
            isoin1SWwet(n,l)=isoin2SWwet(n,l)
            isoin2SWwet(n,l)=0.0
            isoout1SWwet(n,l)=isoout2SWwet(n,l)
            isoout2SWwet(n,l)=0.0
            isowstore1SW(n,l)=isowstore2SW(n,l)
            isowstore2SW(n,l)=0.0
            isoin1IFwet(n,l)=isoin2IFwet(n,l)
            isoin2IFwet(n,l)=0.0
            isoout1IFwet(n,l)=isoout2IFwet(n,l)
            isoout2IFwet(n,l)=0.0
            isowstore1IF(n,l)=isowstore2IF(n,l)
            isowstore2IF(n,l)=0.0
	    endif
	  endif
      end do initialize 

! trish  had to add the divertflg check

      if(divertflg.eq.'y')then
        divert: do m=1,nodivert
          n=gridgive(m)
          i=yyy(n)
          j=xxx(n)
          l=nbasin(i,j)
c          l=1
! NK  Jul. 22/2014         
!         It is possible that the give grid is outside the watershed
!         in which case nbasin(i,j) = 0   >>>  no good, so skip
!         it's not part of the tracer routing then anyways         
          if(l.gt.0)then
            isoin2GW(n,l)=isoin2GW(n,l)+0.75*qr(n)*t
            isoin2IF(n,l)=isoin2IF(n,l)+0.1*qr(n)*t
            isoin2SW(n,l)=isoin2SW(n,l)+0.05*qr(n)*t
          endif
        end do divert
      endif

!       ^^^^^^^^^^^  TRACER ROUTING  ^^^^^^^^^^^^^^^^^^^^^^^^^^
!       isoin    = isotope mass in kg
!       isoout   = isotope mass out in kg
!       isostore = isotope in storage in kg
!       isoconc  = isotope concentration in the unit in kg/m^3
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
      route: do n=1,naa
        lll=next(n)
        
        in_basin: if(lll>0)then       
        inxt=yyy(lll)
        jnxt=xxx(lll)
	  lnxt=nbasin(inxt,jnxt)
	  resnum=ireach(n)
	  rbin=ireach(lll)
        i=yyy(n)
        j=xxx(n)
          
          
          
          l=nbasin(i,j)
c          l=1
          
          
          

!	  WHEN IN LAKE, ADD OPEN WATER EVAPORATION FOR GRID 
!       TO TOTAL EVAPORATION FOR THE LAKE 
        lake_evap: if(resnum>0) then
          isolakest(resnum)=isolakest(resnum)+strloss(n)
        end if lake_evap
             
!       RUN TRACER CODE ONLY IF NOT IN A NON_OUTLET LAKE GRID
        if(resnum<=0.0 .OR.
     *     resnum>0.and.res(n)>0) then


!       WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASIN: SKIP ROUTE
!       IF L=0, THEN THE 1ST GAUGE IN THE BASIN IS NOT AT THE OUTLET, AND THERE
!       IS BASIN AREA BELOW.  B/C WE HAVE NO DATA FOR THIS AREA - IGNORE IT!
        if(slope(n)>0.0.and.l.ne.0)then
          if(store2(n)>0.0)then

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!              ROUTE TRACER THROUGH WETLAND FIRST (IF THERE ARE WETLANDS):
               if(wetland_flag(n)==.true.)then

!                SET INITIAL CONCENTRATION AS 1 kg/m^3]
!                UNITS:  = [KG]
!                TRACE MASS ENTERING WETLAND FROM FLOWPATH
                 isoin2wet(n,l)=isoin2wet(n,l)+1.0*qlz(n)*t
                 isoin2SWwet(n,l)=isoin2SWwet(n,l)
     *                           +1.0*(sumq1(n)+sumq1fs(n))*t
                 isoin2IFwet(n,l)=isoin2IFwet(n,l)
     *                           +1.0*(sumqint(n)+sumqintfs(n))*t


                 

                 if(qowet2(n).ge.0.0)then
!                  SPECIFY INITIAL MASS FOR CHANNEL FROM WETLAND:
                   call tracewet(iz,jz,time,n,l,t)
                   isoin2GW(n,l)=isoin2GW(n,l)+isooutwet1
                   isoin2SW(n,l)=isoin2SW(n,l)+isooutwet2
                   isoin2IF(n,l)=isoin2IF(n,l)+isooutwet3

	           else   ! qowet2(n) is -ve
!                  SPECIFY INITIAL MASS IN WETLAND (+MASS FROM CHANNEL)
!                  NOTE: -VE CAUSE QOWET2 IS -VE BUT MASS IS ADDED!

                   isoin2wet(n,l)=isoin2wet(n,l)  
     *                           -isoconcGW(n,l)*qowet2(n)*t
                   isoin2SWwet(n,l)=isoin2SWwet(n,l)
     *                             -isoconcSW(n,l)*qowet2(n)*t
                   isoin2IFwet(n,l)=isoin2IFwet(n,l)
     *                             -isoconcIF(n,l)*qowet2(n)*t

                   isoin2GW(n,l)=isoin2GW(n,l)
     *                          +isoconcGW(n,l)*qowet2(n)*t
                   isoin2SW(n,l)=isoin2SW(n,l)
     *                          +isoconcSW(n,l)*qowet2(n)*t
                   isoin2IF(n,l)=isoin2IF(n,l)
     *                          +isoconcIF(n,l)*qowet2(n)*t
                   call tracewet(iz,jz,time,n,l,t)
	           endif

	         else
!                SET INITIAL CONCENTRATION AS 1 kg/m^3
!                UNITS:  = [KG]
	           isoin2GW(n,l)=isoin2GW(n,l)+1.0*qlz(n)*t
                 isoin2SW(n,l)=isoin2SW(n,l)
     *                        +1.0*(sumq1(n)+sumq1fs(n))*t
                 isoin2IF(n,l)=isoin2IF(n,l)
     *                        +1.0*(sumqint(n)+sumqintfs(n))*t
     
!	           AT LAKE OUTLET, MAKE INFLOW EQUAL THE ACCUMULATED INFLOWS TO WHOLE LAKE
	           if(resnum>0.AND.res(n)>0)then
	               isoin2GW(n,l)=isolakeGW(resnum)
	               isoin2SW(n,l)=isolakeSW(resnum)
	               isoin2IF(n,l)=isolakeIF(resnum)
	           end if
	               
	        endif


! NO WETLANDS AND QOWET2 IS +VE (OR ALREADY REROUTED THRU WETLAND)
!            CHANNEL ROUTING OR WETLAND->CHANNEL ROUTING
!            ROUTE TRACER THROUGH CHANNEL - CONVERGENCE LOOP:

!            UPDATE THE TRACER MASSES IN CHANNEL     

!            COMPUTE THE TRACER MASSES LEAVING THE GRID...
!            EVAP LOSS ADDED TO PRESERVE FLOW BALANCES.
             if(resnum>0.AND.res(n)>0)then
               isostore2GW(n,l)=(isostore1GW(n,l)+(isoin1GW(n,l)+
     *               isoin2GW(n,l)-isoout1GW(n,l))/2.)
     *              /(1.0+(qo2(n)*coeff(n)+isolakest(resnum))*t
     *              /2.0/store2(n))
	         isostore2GW(n,l)=amax1(isostore2GW(n,l),0.00001)
	         isoconcGW(n,l)=isostore2GW(n,l)/store2(n)
               isoout2GW(n,l)=isoconcGW(n,l)*qo2(n)*t*coeff(n)
     *                     +isoconcGW(n,l)*isolakest(resnum)*t
     	         isoout2GW(n,l)=amax1(isoout2GW(n,l),0.00001)
     	         
     	         isostore2SW(n,l)=(isostore1SW(n,l)+(isoin1SW(n,l)+
     *               isoin2SW(n,l)-isoout1SW(n,l))/2.)
     *               /(1.0+(qo2(n)*coeff(n)+isolakest(resnum))*t
     *              /2.0/store2(n))
	         isostore2SW(n,l)=amax1(isostore2SW(n,l),0.00001)
	         isoconcSW(n,l)=isostore2SW(n,l)/store2(n)
     	         isoout2SW(n,l)=isoconcSW(n,l)*qo2(n)*t*coeff(n)
     *                     +isoconcSW(n,l)*isolakest(resnum)*t
     	         isoout2SW(n,l)=amax1(isoout2SW(n,l),0.00001)
     	         
     	         isostore2IF(n,l)=(isostore1IF(n,l)+(isoin1IF(n,l)+
     *               isoin2IF(n,l)-isoout1IF(n,l))/2.)
     *               /(1.0+(qo2(n)*coeff(n)+isolakest(resnum))*t
     *              /2.0/store2(n))
	         isostore2IF(n,l)=amax1(isostore2IF(n,l),0.00001)
	         isoconcIF(n,l)=isostore2IF(n,l)/store2(n)  
     	         isoout2IF(n,l)=isoconcIF(n,l)*qo2(n)*t*coeff(n)
     *                     +isoconcIF(n,l)*isolakest(resnum)*t
     	         isoout2IF(n,l)=amax1(isoout2IF(n,l),0.00001)
     	           
             else
               isostore2GW(n,l)=(isostore1GW(n,l)+(isoin1GW(n,l)+
     *               isoin2GW(n,l)-isoout1GW(n,l))/2.)
     *              /(1.0+(qo2(n)*coeff(n)+strloss(n))*t
     *              /2.0/store2(n))
	         isostore2GW(n,l)=amax1(isostore2GW(n,l),0.00001)
	         isoconcGW(n,l)=isostore2GW(n,l)/store2(n)
               isoout2GW(n,l)=isoconcGW(n,l)*qo2(n)*t*coeff(n)
     *                       +isoconcGW(n,l)*strloss(n)*t
       	       isoout2GW(n,l)=amax1(isoout2GW(n,l),0.00001)

               isostore2SW(n,l)=(isostore1SW(n,l)+(isoin1SW(n,l)+
     *               isoin2SW(n,l)-isoout1SW(n,l))/2.)
     *               /(1.0+(qo2(n)*coeff(n)+strloss(n))*t
     *              /2.0/store2(n))
	         isostore2SW(n,l)=amax1(isostore2SW(n,l),0.00001)
	         isoconcSW(n,l)=isostore2SW(n,l)/store2(n)
               isoout2SW(n,l)=isoconcSW(n,l)*qo2(n)*t*coeff(n)
     *                       +isoconcSW(n,l)*strloss(n)*t
     	         isoout2SW(n,l)=amax1(isoout2SW(n,l),0.00001)

               isostore2IF(n,l)=(isostore1IF(n,l)+(isoin1IF(n,l)+
     *               isoin2IF(n,l)-isoout1IF(n,l))/2.)
     *               /(1.0+(qo2(n)*coeff(n)+strloss(n))*t
     *              /2.0/store2(n))
	         isostore2IF(n,l)=amax1(isostore2IF(n,l),0.00001)
	         isoconcIF(n,l)=isostore2IF(n,l)/store2(n)
               isoout2IF(n,l)=isoconcIF(n,l)*qo2(n)*t*coeff(n)
     *                       +isoconcIF(n,l)*strloss(n)*t
     	         isoout2IF(n,l)=amax1(isoout2IF(n,l),0.00001)
     	         
     	       end if
     	              
            isosumQ(n)=isoout2GW(n,l)+isoout2IF(n,l)+isoout2SW(n,l)
   
          else   ! follow what's done in rerout
!           TS: added May 21/08 s.t tracer works when store2<0
            if(isoin2GW(n,l).lt.0.0)then
              isoout2GW(n,l)=isoin2GW(n,l)/2.0
            else
              isoout2GW(n,l)=0.0
	      endif
            isostore2GW(n,l)=isostore1GW(n,l)+(isoin1GW(n,l)+
     *            isoin2GW(n,l)-isoout1GW(n,l)-isoout2GW(n,l))/2.
	      isostore2GW(n,l)=amax1(isostore2GW(n,l),0.001)
            isostore2SW(n,l)=0.0
            isoout2SW(n,l)=0.0
            isostore2IF(n,l)=0.0
            isoout2IF(n,l)=0.0
          end if
!       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
!         ISOTOPES ARE CONSERVATIVE!
!         NOTE: CAN USE THIS TO MODEL GW RETARDATION OF FLOW

!         ADD MASS OUTFLOW OF CURRENT GRID AS INFLOW TO LAKE (NEXT GRID) (cms)
!         BUT ONLY IF ON THE EDGE OF A LAKE (NOT IN A LAKE)
          if(rbin>0) then
           if(nopt(l)==2.and.n==iflowgrid(l))then
            isolakeGW(rbin)=isolakeGW(rbin)
     *            +qhyd(l,iz+1)*t*isoconcGW(n,l)
            isolakeSW(rbin)=isolakeSW(rbin)
     *            +qhyd(l,iz+1)*t*isoconcSW(n,l)
            isolakeIF(rbin)=isolakeIF(rbin)
     *            +qhyd(l,iz+1)*t*isoconcIF(n,l)
           else
            isolakeGW(rbin)=isolakeGW(rbin)+isoout2GW(n,l)
            isolakeSW(rbin)=isolakeSW(rbin)+isoout2SW(n,l)
            isolakeIF(rbin)=isolakeIF(rbin)+isoout2IF(n,l)

	     endif
	    end if


          if(lnxt/=0) then
          
           if(nopt(l)==2.and.n==iflowgrid(l))then
           
              isoin2GW(lll,lnxt)=isoin2GW(lll,lnxt)
     *            +qhyd(l,iz+1)*t*isoconcGW(n,l)
              isoin2IF(lll,lnxt)=isoin2IF(lll,lnxt)
     *            +qhyd(l,iz+1)*t*isoconcIF(n,l)
              isoin2SW(lll,lnxt)=isoin2SW(lll,lnxt)
     *            +qhyd(l,iz+1)*t*isoconcSW(n,l)
     
           else
            if(resnum>0.AND.res(n)>0)then
            
             isoin2GW(lll,lnxt)=isoin2GW(lll,lnxt)
     *       +isoout2GW(n,l)-isoconcGW(n,l)*isolakest(resnum)*t
     
             isoin2IF(lll,lnxt)=isoin2IF(lll,lnxt)
     *       +isoout2IF(n,l)-isoconcIF(n,l)*isolakest(resnum)*t
     
             isoin2SW(lll,lnxt)=isoin2SW(lll,lnxt)
     *       +isoout2SW(n,l)-isoconcSW(n,l)*isolakest(resnum)*t
     
            else
             
             isoin2GW(lll,lnxt)=isoin2GW(lll,lnxt)
     *       +isoout2GW(n,l)-isoconcGW(n,l)*strloss(n)*t
     
             isoin2IF(lll,lnxt)=isoin2IF(lll,lnxt)
     *       +isoout2IF(n,l)-isoconcIF(n,l)*strloss(n)*t
     
             isoin2SW(lll,lnxt)=isoin2SW(lll,lnxt)
     *       +isoout2SW(n,l)-isoconcSW(n,l)*strloss(n)*t
     
            end if
           end if
          endif


!         ERROR CHECK: ALL COMPONENTS MUST SUM TO TOTAL SIM'D Q
!         CAUSE SOME CONC'NS > 1.0 FROM ABOVE...
          isosumQ(n)=isoout2SW(n,l)+isoout2IF(n,l)+isoout2GW(n,l)

 !        ACCUMULATE MASSES FOR MASS BALANCE CHECK AT END OF RUN:
          massin(n)=(isoin1GW(n,l)+isoin2GW(n,l))/2.
     *             +(isoin1SW(n,l)+isoin2SW(n,l))/2.
     *             +(isoin1IF(n,l)+isoin2IF(n,l))/2.
          massout(n)=(isoout1GW(n,l)+isoout2GW(n,l))/2.
     *              +(isoout1SW(n,l)+isoout2SW(n,l))/2.
     *              +(isoout1IF(n,l)+isoout2IF(n,l))/2.
          masstore(n)=(isostore2GW(n,l)-isostore1GW(n,l))
     *               +(isostore2SW(n,l)-isostore1SW(n,l))
     *               +(isostore2IF(n,l)-isostore1IF(n,l))
          ISOdelta(n)=massin(n)-massout(n)
 
 !        CHECK OVERALL MASS BALANCE:
 !        MASS IN - MASS OUT = MASS STORED (IF EVERYTHING WORKS!) 
          if(iz.ne.jz.and.n.eq.nnprint)then
            sqerr=(ISOdelta(n)-masstore(n))**2.
            write(91,8000)time,massin(n),massout(n),ISOdelta(n),
     *                    masstore(n),sqerr
	    endif
	    
        endif  ! SLOPE>=0 AND L.NE.0

        else   ! IREACH(N).gt.0 but res(n)=0
!         IN A LAKE, PASS UPSTREAM GROUNDWATER TO D/S GRID


! trish   rbin was 0 here so crashed when running debug  NK 20140626
! NK  Jul. 22/2014         
!         It is possible that the give grid is outside the watershed
!         in which case l = nbasin(i,j) = 0   >>>  no good, so skip
!         it's not part of the tracer routing then anyways         
            if(rbin.gt.0.and.l.gt.0)then
	        isoout2GW(n,l)=isolakeGW(rbin) 
	      endif

	  endif

	  tt(n)=rl(n)*chaxa(n)/qo2(n)
c	 if(n==595)write(230,*)isolakeGW(19),t,iz
c     * ,jz

        end if in_basin         !  lll<0
     
      end do route              ! GRID NO. LOOP
     

!     WRITE OUT DATA TO TRACER.CSV FILE
!      if(jan.eq.3)then    ! Prints every 0.25 hours
c      if(iz.ne.jz.and.dds_flag.ne.1)then     ! Prints every hour
      if(iz.ne.jz)then     ! Prints every hour
c        do l=1,no
c!     rev. 9.7.16  Jan.  05/11  - NK: Fixed init flows outside sub-basin
c          if(inbsnflg(l).eq.1)then    ! added nk Jan. 05/11
c!           STORES THE GRID NO. AS FXN OF THE GAUGE(BSN) NO.
c            nn(l)=s(iy(l),jx(l))
c          else
c!     rev. 9.8.18  Apr.  26/11  - NK: Added in-basin check in tracer4
c            nn(l)=-1
c          endif
c	  end do
!       WRITE TO THE FILE ONLY WHEN YOU ARE AT THE FLAGGED OUTPUT
!       GRID=NNPRINT. 
!       UNITS: isoout2IBN(n,ll)=[kg]/[s]=[kg/s]=[m^3/s]
!       NB: TRACER OUTFLOW BASED ON MASS TO ACCOUNT FOR EVAP. LOSS

!Trish:
! problems for partial watersheds   fix fix fix  nk
c      if(mod(int(totaltime),kt).eq.0.and.nn(l).gt.0)

      if(firstpass)then
          do l=1,no
              if(nn(l).le.0)trace_off=.true.
          end do
      endif


        if(mod(int(totaltime),kt).eq.0)then

            
         write(90,9014)totaltime,
     *     (qo2(nn(l))*isoconcSW(nn(l),l),                   !isoconcSW(nn(l),l),
     *      qo2(nn(l))*isoconcIF(nn(l),l),                   !isoconcIF(nn(l),l),
     *      qo2(nn(l))*isoconcGW(nn(l),l),l=1,no)            !isoconcGW(nn(l),l),l=1,no)

         endif

cc	  if(dds_flag.ne.1)
cd          write(999,12345)time,
cd    *          (qo2(nn(l)),isosumQ(nn(l))/t,
cd    *          (qo2(nn(l))-isosumQ(nn(l))/t),l=1,no)


      endif

      firstpass=.false.

! FORMATS:
 8000 FORMAT(f10.3,256(',',f15.5))
 9004 FORMAT(f10.3,<no>(2(',',I15),10(',',f15.5)))
 9014 FORMAT(f10.3,<no>(3(',',f15.5)))
12345 FORMAT(f10.3,<no>(3(',',f15.5)))


      RETURN
      END SUBROUTINE tracer4