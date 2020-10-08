      SUBROUTINE tracer(iz,jz,time,t,jan,tdum,irte)

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen, Tricia Stadnyk and Tegan Holmes
        
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
! PROGRAM TRACER
! This subroutine is designed to route isotope tracers through the
! WATFLOOD model.  It will do so based on a mixed cell model (continuity).
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

      INTEGER :: iz,jz,jan,n,i,j,l,ii,jjz,jjzold,rbin,irte
	REAL*4  :: time,tdum,t,ustar,Pe,Cr
      real*4  :: tt1,elat,elon1,elon2,lm,fr,hydd

      DATA jjzold/0/

      jjzold=jjz
      jjz=jz
	if(jjz.lt.1) jjz=1

	if(icount.eq.0)then
c        allocate(nn(no),vel(naa),disp(naa),coeff(naa),stat=iAll)
        allocate(nn(no),vel(na),disp(na),coeff(na),stat=iAll)
	  if(iAll.ne.0) STOP 'Error allocating arrays in Tracer.for'
      endif

!     CALCULATE THE DISPERSION COEFFICIENT FOR TRACER MIXING
!***TS:  to check -- slope(n) already sqrt'd
        do n=1,naa
c          if(slope(n).gt.0.0.and.hcha2(n).gt.0.001.and.ireach(n).le.0
c     *     .or.ireach(n).gt.0.and.res(n).eq.0)then     nk
          if(slope(n).gt.0.0.and.hcha2(n).gt.0.001.and.ireach(n).le.0
     *     .or.ireach(n).gt.0.and.res(n).eq.0.and.hcha2(n).gt.0.001)then
            
            vel(n)=qo2(n)/chaxa(n)
            disp(n)=vel(n)*al/2.
            tt1=al*chaxa(n)/qo2(n)
            hcha2(n)=amax1(hcha2(n),0.0)  ! nk

            ustar=sqrt(9.81*hcha2(n))*slope(n)   !slope(n) already sqrt'd
!            print*,vel(n),ustar,hcha2(n)
!            print*,disp(n),tt1

	      coeff(n)=2*vel(n)*t/(al) !(chawid(n)/vel(n))
!            print*,coeff(n)
!            PAUSE
	      if(coeff(n).gt.1.0) coeff(n)=1.0
            if(coeff(n).lt.0) coeff(n)=0.0
            
          else
	      coeff(n)=0.0
          endif

!         DETERMINE THE MIXING LENGTH
!         ELON = longitudinal dispersion (m^2/s)
!         ELAT = lateral dispersion (m^2/s)
          elon1=0.011*vel(n)**2*chawid(n)**2/hcha2(n)/ustar  ! std

          hcha2(n)=amax1(hcha2(n),0.0)   ! nk
          hydd=(2/3.*chawid(n)*hcha2(n))/   ! hydraulic depth for parabolic channel
     *         (chawid(n)*sqrt(hcha2(n)/chadep(n)))

	    fr=vel(n)/sqrt(9.81*hcha2(n))         ! Froude number

	    elon2=0.05937*qo2(n)/slope(n)**2/chawid(n)   ! Fr<0.5
          elat=0.6*hcha2(n)*ustar  !
	    lm=0.1*vel(n)*chawid(n)**2/elat
          Pe=al*vel(n)/elon1
	    Cr=vel(n)*t/al

d         if(n.eq.nnprint)write(120,1200)time,n,qo2(n),hcha2(n),
d    *    chadep(n),vel(n),hydd,fr,elon1,elon2,elat,lm,Pe,Cr

!	    if(dds_flag.ne.1) write(888,*)n,coeff(n),vel(n),disp(n),Pe,Cr
          coeff(n)=1.0  ! only for Gr2k!

        end do






! FOR NOW, HARD CODE THE TYPE OF TRACER BEING USED
! SET IN AREATRC SO THAT SPL CAN READ IT...
! 0 = SUB-GAUGE TRACER
! 1 = GLACIER MELT TRACER (ii=6)
! 2 = LANDCOVER TRACER
! 3 = RAIN-ON-STREAM TRACER AS FXN OF SUB-BASIN
! 4 = FLOW TYPE TRACER (SW+IF+GW) AS FXN OF SUB-BASIN
! 5 = SNOWMELT TRACER (SW+IF) AS FXN OF SUB-BASIN
!100= ORIGINAL GW TRACER (NK) AS FXN OF SUB-BASIN
!101= WETLAND FLOW TRACER (qowet2)
! **********************************************************************
      if(itrace.eq.0)then

	  print*,'Sub-basin tracer under development'
	  print*,'Select a different itrace number and re-run'
	  STOP 'Program terminated in trace.for @0'

        if(icount.eq.0)then
!         ALLOCATE THE SUB-BASIN TRACER ARRAYS:
          allocate(isoin1IBN(na),isoin2IBN(na),
     *    isoout1IBN(na),isoout2IBN(na),isoconcIBN(na),
     *    isostore1IBN(na),isostore2IBN(na),
     *    massin(na),massout(na),masstore(na),ISOdelta(na),stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating IBN Tracer arrays'

!         ALLOCATE WETLAND ARRAYS IF WETLAND OPTION TURNED ON
          if(wetflg.eq.'y')
     *     allocate(isoin1wet(na,no),isoin2wet(na,no),isoout1wet(na,no),
     *     isoout2wet(na,no),isoconcwet(na,no),isowstore1(na,no),
     *     isowstore2(na,no),wmassin(na),wmassout(na),wmasstore(na),
     *     wISOdelta(na),stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating Wetland Tracer arrays'

!         ALLOCATE LAKE GW TRACER
          allocate(isolakeGW(noresv),stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating Lake Tracer array'          

          do n=1,naa
	      i=yyy(n)
	      j=xxx(n)
	      l=nbasin(i,j)
!           INITIALIZE 'OUT' ARRAYS SO 'IN' ARRAYS SET TO ZERO @t=1:
            isoin2IBN(n)=0.0
            isoout2IBN(n)=0.0
            isostore2IBN(n)=0.0
!           INITIALIZE WETLAND ARRAYS
            if(wetland_flag(n))then
              isoin2wet(n,l)=0.0
              isoout2wet(n,l)=0.0
              isowstore2(n,l)=0.0
            endif
	    end do
        endif

!       WRITE OUTPUT FILE HEADERS
	  if(icount.eq.0)write(91,11000)
        if(icount.eq.0)write(92,50000)

	  call tracer0(iz,jz,time,t,jan,tdum)

! **********************************************************************
	elseif(itrace.eq.1)then
c	  print*,'Landcover tracer under development'
c	  print*,'Select a different itrace number and re-run'
c	  STOP 'Program terminated in trace.for @2'

        if(icount.eq.0)then 
!         ALLOCATE THE BC (BASIN-COVER) TRACER ARRAYS:
          allocate(isoin1LC(na,no,classcount-1),
     *    isoin2LC(na,no,classcount-1),
     *    isoout1LC(na,no,classcount-1),isoout2LC(na,no,classcount-1),
     *    isoconcLC(na,no,classcount-1),isoconcBLZS(na,no,classcount-1),
     *    isostore1LC(na,no,classcount-1),
     *    isostore2LC(na,no,classcount-1),
     *    isoin1BLZS(na,no,classcount-1),isoin2BLZS(na,no,classcount-1),
     *    isoout1BLZS(na,no,classcount-1),
     *    isoout2BLZS(na,no,classcount-1),
     *    isoBLZS1(na,no,classcount-1),isoBLZS2(na,no,classcount-1),
     *    isosum(na,no),
     *    massin(na),massout(na),masstore(na),ISOdelta(na),stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating B-LC Tracer arrays'
          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            l=nbasin(i,j)	    
	      if(l.ne.0)then
	        do ii=1,classcount-1
!              INITIALIZE 'OUT' ARRAYS SO 'IN' ARRAYS SET TO ZERO @t=1:
                isoin2LC(n,l,ii)=0.0
                isoout2LC(n,l,ii)=0.0
                isostore2LC(n,l,ii)=0.0
                isoin2BLZS(n,l,ii)=0.0
	          isoout2BLZS(n,l,ii)=0.0
                isoBLZS2(n,l,ii)=isoBLZS2(n,l,ii)*1.0*(drng(n,ii)+
     *                    drngfs(n,ii))/1000.*aclass(n,ii)*frac(n)*al**2
              end do
            endif
	    end do
        endif

!       WRITE OUTPUT FILE HEADERS
        if(icount.eq.0)write(91,11000)
        if(icount.eq.0)write(92,52000)

	  call tracer1(iz,jz,time,t,jan,tdum)


! **********************************************************************
	elseif(itrace.eq.2)then

c	  print*,'Landcover tracer under development'
c	  print*,'Select a different itrace number and re-run'
c	  STOP 'Program terminated in trace.for @2'

        if(icount.eq.0)then 
!         ALLOCATE THE BC (BASIN-COVER) TRACER ARRAYS:
          allocate(isoin1LC(na,no,classcount-1),
     *     isoin2LC(na,no,classcount-1),
     *    isoout1LC(na,no,classcount-1),isoout2LC(na,no,classcount-1),
     *    isoconcLC(na,no,classcount-1),isoconcBLZS(na,no,classcount-1),
     *    isostore1LC(na,no,classcount-1),
     *    isostore2LC(na,no,classcount-1),
     *    isoin1BLZS(na,no,classcount-1),isoin2BLZS(na,no,classcount-1),
     *    isoout1BLZS(na,no,classcount-1),
     *    isoout2BLZS(na,no,classcount-1),
     *    isoBLZS1(na,no,classcount-1),isoBLZS2(na,no,classcount-1),
     *    isosum(na,no),
     *    massin(na),massout(na),masstore(na),ISOdelta(na),stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating B-LC Tracer arrays'
          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            l=nbasin(i,j)	    
	      if(l.ne.0)then
	        do ii=1,classcount-1
!              INITIALIZE 'OUT' ARRAYS SO 'IN' ARRAYS SET TO ZERO @t=1:
                isoin2LC(n,l,ii)=0.0
                isoout2LC(n,l,ii)=0.0
                isostore2LC(n,l,ii)=0.0
                isoin2BLZS(n,l,ii)=0.0
	          isoout2BLZS(n,l,ii)=0.0
                isoBLZS2(n,l,ii)=isoBLZS2(n,l,ii)*1.0*(drng(n,ii)+
     *                    drngfs(n,ii))/1000.*aclass(n,ii)*frac(n)*al**2
              end do
            endif
	    end do
        endif

!       WRITE OUTPUT FILE HEADERS
        if(icount.eq.0)write(91,11000)
        if(icount.eq.0)write(92,52000)

	  call tracer2(iz,jz,time,t,jan,tdum)

! **********************************************************************
	elseif(itrace.eq.3)then

	  print*,'Rain tracer under development'
	  print*,'Select a different itrace number and re-run'
	  STOP 'Program terminated in trace.for @3'

        if(icount.eq.0)then 
!         ALLOCATE THE GW TRACER ARRAYS:
          allocate(isoin1P(na,no),isoin2P(na,no),isoout1P(na,no),
     *    isoout2P(na,no),isoconcP(na,no),isostore1P(na,no),
     *    isostore2P(na,no),massin(na),massout(na),masstore(na),
     *    ISOdelta(na),stat=iAll)
          if(iAll.ne.0) STOP 'Error allocating RAIN Tracer arrays'
          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            l=nbasin(i,j)	    
!           INITIALIZE 'OUT' ARRAYS SO 'IN' ARRAYS SET TO ZERO @t=1:
!           INITIALIZE STORAGE TO SOME BACKGROUND CONC'N OF TRACER:
	      if(l.ne.0)then
              isoin2P(n,l)=0.0
              isoout2P(n,l)=0.0
              isostore2P(n,l)=0.0
              isostore2P(n,l)=0.0
            endif
          end do
        endif

!       WRITE OUTPUT FILE HEADERS
        if(icount.eq.0)write(91,11000)
        if(icount.eq.0)write(92,53000)

	  call tracer3(iz,jz,time,t,jan,tdum)

! **********************************************************************
	elseif(itrace.eq.4)then
      
        if(icount.eq.0)then 

!         ALLOCATE THE FLOW-TYPE TRACER ARRAYS:
          allocate(isoin1GW(na,no),isoin2GW(na,no),isoout1GW(na,no),
     *    isoout2GW(na,no),isoconcGW(na,no),
     *    isostore1GW(na,no),isostore2GW(na,no),isoin1SW(na,no),
     *    isoin2SW(na,no),isoout1SW(na,no),isoout2SW(na,no),
     *    isoconcSW(na,no),isostore1SW(na,no),
     *    isostore2SW(na,no),isoin1IF(na,no),isoin2IF(na,no),
     *    isoout1IF(na,no),isoout2IF(na,no),
     *    isoconcIF(na,no),isostore1IF(na,no),isostore2IF(na,no),
     *    massin(na),massout(na),masstore(na),ISOdelta(na),
     *    isosumQ(na),tt(na),stat=iAll)
		if(iAll.ne.0) STOP 'Error allocating FLOW-TYPE Tracer arrays'

!         ALLOCATE WETLAND ARRAYS IF WETLAND OPTION TURNED ON
          if(wetflg.eq.'y')
     *     allocate(isoin1wet(na,no),isoin2wet(na,no),isoout1wet(na,no),
     *     isoout2wet(na,no),isoconcwet(na,no),isowstore1(na,no),
     *     isowstore2(na,no),isoin1SWwet(na,no),isoin2SWwet(na,no),
     *     isoout1SWwet(na,no),isoout2SWwet(na,no),isoconcSWwet(na,no),
     *     isowstore1SW(na,no),isowstore2SW(na,no),isoin1IFwet(na,no),
     *     isoin2IFwet(na,no),isoout1IFwet(na,no),isoout2IFwet(na,no),
     *     isoconcIFwet(na,no),isowstore1IF(na,no),isowstore2IF(na,no),
     *     wmassin(na),wmassout(na),wmasstore(na),wISOdelta(na),
     *     stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating Wetland Tracer arrays'

!         ALLOCATE LAKE TRACER ARRAYS
          allocate(isolakeGW(noresv),
     *    isolakest(noresv),isolakeSW(noresv),
     *    isolakeIF(noresv),stat=iAll)

	print*,'4 isolakeGW allocated with',noresv

	    if(iAll.ne.0) STOP 'Error allocating Lake Tracer array'          

!         INITIALIZE PRELIMINARY ITERATION VALUES
          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            l=nbasin(i,j)
!           INITIALIZE 'OUT' ARRAYS SO 'IN' ARRAYS SET TO ZERO @t=1
!           BASED ON INITIALIZATION LOOP BELOW WHERE 1'S = 2'S:
!           INITIALIZE STORAGE TO SOME BACKGROUND CONC'N OF TRACER:
	      if(l.ne.0)then
              isoconcGW(n,l)=0.0
              isoin2GW(n,l)=0.0
              isoout2GW(n,l)=0.0
              isostore2GW(n,l)=0.75*store2(n)
              isoconcSW(n,l)=0.0
              isoin2SW(n,l)=0.0
              isoout2SW(n,l)=0.0
	        isostore2SW(n,l)=0.0
              isoconcIF(n,l)=0.0
              isoin2IF(n,l)=0.0
              isoout2IF(n,l)=0.0
	        isostore2IF(n,l)=0.0
	        isosumQ(n)=0.0
!             INITIALIZE WETLAND ARRAYS
  	        if(wetland_flag(n))then
                isoconcwet(n,l)=0.0
                isoin2wet(n,l)=0.0
                isoout2wet(n,l)=0.0
                isowstore2(n,l)=0.75*wstore2(n)
                isoconcSWwet(n,l)=0.0
                isoin2SWwet(n,l)=0.0
                isoout2SWwet(n,l)=0.0
                isowstore2SW(n,l)=0.0
                isoconcIFwet(n,l)=0.0
                isoin2IFwet(n,l)=0.0
                isoout2IFwet(n,l)=0.0
                isowstore2IF(n,l)=0.0
              endif
            endif
          end do
        endif

!       WRITE OUTPUT FILE HEADERS        
        if(icount.eq.0)write(91,11000)
        if(icount.eq.0)write(90,57000)        
                
        call tracer4(iz,jz,time,t,jan,tdum,jjz,jjzold)

!       FOR ISOTOPE MODELLING
!       + rd_iso.for to read in isotope data

d       if(iopt.eq.2)print*,
d    *   'tracer.f - line 387, about to call isotopes.f'

!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
!     This call moved to sub
!	  if(frcflg.eq.'y') call isotopes(iz,jz,jjz,time,t,irte)

! **********************************************************************
	elseif(itrace.eq.5)then

        if(icount.eq.0)then 
!         ALLOCATE THE FLOW-TYPE TRACER ARRAYS:
          allocate(isoin1SW(na,no),isoin2SW(na,no),isoout1SW(na,no),
     *    isoout2SW(na,no),isoconcSW(na,no),isostore1SW(na,no),
     *    isostore2SW(na,no),isoin1SWfs(na,no),isoin2SWfs(na,no),
     *    isoout1SWfs(na,no),isoout2SWfs(na,no),isoconcSWfs(na,no),
     *    isostore1SWfs(na,no),isostore2SWfs(na,no),isoin1IF(na,no),
     *    isoin2IF(na,no),isoout1IF(na,no),isoout2IF(na,no),
     *    isoconcIF(na,no),isostore1IF(na,no),isostore2IF(na,no),
     *    isoin1IFfs(na,no),isoin2IFfs(na,no),isoout1IFfs(na,no),
     *    isoout2IFfs(na,no),isoconcIFfs(na,no),isostore1IFfs(na,no),
     *    isostore2IFfs(na,no),isoin1GW(na,no),isoin2GW(na,no),
     *    isoout1GW(na,no),isoout2GW(na,no),isoconcGW(na,no),
     *    isostore1GW(na,no),isostore2GW(na,no),isoin1GWfs(na,no),
     *    isoin2GWfs(na,no),isoout1GWfs(na,no),isoout2GWfs(na,no),
     *    isoconcGWfs(na,no),isostore1GWfs(na,no),isostore2GWfs(na,no),
     *    isoLZS1(na,no),isoLZS2(na,no),isoconcLZS(na,no),
     *    isoin1LZS(na,no),isoin2LZS(na,no),isoout1LZS(na,no),
     *    isoout2LZS(na,no),isoconcLZSfs(na,no),isoin1LZSfs(na,no),
     *    isoin2LZSfs(na,no),isoout1LZSfs(na,no),isoout2LZSfs(na,no),
     *    isoLZS1fs(na,no),isoLZS2fs(na,no),massin(na),massout(na),
     *    masstore(na),ISOdelta(na),isosumQ(na),isosumQfs(na),
     *    tt(na),stat=iAll)
          if(iAll.ne.0) STOP 'Error allocating SNOWMELT Tracer arrays'

!         ALLOCATE WETLAND ARRAYS IF WETLAND OPTION TURNED ON
          if(wetflg.eq.'y')
     *     allocate(isoin1wet(na,no),isoin2wet(na,no),isoout1wet(na,no),
     *     isoout2wet(na,no),isoconcwet(na,no),isowstore1(na,no),
     *     isowstore2(na,no),isoin1SWwet(na,no),isoin2SWwet(na,no),
     *     isoout1SWwet(na,no),isoout2SWwet(na,no),isoconcSWwet(na,no),
     *     isowstore1SW(na,no),isowstore2SW(na,no),isoin1IFwet(na,no),
     *     isoin2IFwet(na,no),isoout1IFwet(na,no),isoout2IFwet(na,no),
     *     isoconcIFwet(na,no),isowstore1IF(na,no),isowstore2IF(na,no),
     *     isoin1fswet(na,no),isoin2fswet(na,no),isoout1fswet(na,no),
     *     isoout2fswet(na,no),isoconcfswet(na,no),isowstore1fs(na,no),
     *     isowstore2fs(na,no),
     *     isoin1SWfswet(na,no),isoin2SWfswet(na,no),
     *     isoout1SWfswet(na,no),isoout2SWfswet(na,no),
     *     isoconcSWfswet(na,no),isowstore1SWfs(na,no),
     *     isowstore2SWfs(na,no),isoin1IFfswet(na,no),
     *     isoin2IFfswet(na,no),isoout1IFfswet(na,no),
     *     isoout2IFfswet(na,no),isoconcIFfswet(na,no),
     *     isowstore1IFfs(na,no),isowstore2IFfs(na,no),
     *     wmassin(na),wmassout(na),wmasstore(na),wISOdelta(na),
     *     stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating Wetland Tracer arrays'
          if(iopt.gt.1) print*,'allocations completed'

!         ALLOCATE LAKE GW TRACER
          allocate(isolakeGW(noresv),isolakest(noresv),
     *             stat=iAll)

	print*,'5 isolakeGW allocated with',noresv

	    if(iAll.ne.0) STOP 'Error allocating Lake Tracer array'          

          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            l=nbasin(i,j)
!           INITIALIZE 'OUT' ARRAYS SO 'IN' ARRAYS SET TO ZERO @t=1
!           BASED ON INITIALIZATION LOOP BELOW WHERE 1'S = 2'S:
!           INITIALIZE STORAGE TO SOME BACKGROUND CONC'N OF TRACER:
	      if(l.ne.0)then
	        isosumQ(n)=0.0
	        isosumQfs(n)=0.0
	        isoconcSW(n,l)=0.0
              isoin2SW(n,l)=0.0
              isoout2SW(n,l)=0.0
              isostore2SW(n,l)=0.0
	        isoconcSWfs(n,l)=0.0
              isoin2SWfs(n,l)=0.0
              isoout2SWfs(n,l)=0.0
              isostore2SWfs(n,l)=0.0
	        isoconcIF(n,l)=0.0
              isoin2IF(n,l)=0.0
              isoout2IF(n,l)=0.0
              isostore2IF(n,l)=0.0
	        isoconcIFfs(n,l)=0.0
              isoin2IFfs(n,l)=0.0
              isoout2IFfs(n,l)=0.0
              isostore2IFfs(n,l)=0.0
	        isoconcGW(n,l)=0.0
              isoin2GW(n,l)=0.0
              isoout2GW(n,l)=0.0
              isostore2GW(n,l)=0.0*store2(n)
	        isoconcGWfs(n,l)=0.0
              isoin2GWfs(n,l)=0.0
              isoout2GWfs(n,l)=0.0
              isostore2GWfs(n,l)=0.0
	        isoconcLZS(n,l)=0.0
              isoin2LZS(n,l)=0.0
              isoout2LZS(n,l)=0.0
              isoLZS2(n,l)=0.0
	        isoconcLZSfs(n,l)=0.0
              isoin2LZSfs(n,l)=0.0
              isoout2LZSfs(n,l)=0.0
              isoLZS2fs(n,l)=0.0

!             SUM ALL CONTRIBUTIONS FROM ALL LANDCLASSES TO GET TOTAL INITIAL STORAGE
!             FOR UZS AND LZS => CAN'T INITIALIZE CHANNEL/SURFACE STORAGES!!
!              do ii=1,classcount-1
!                isoLZS2(n,l)=isoLZS2(n,l)+1.0*drng(n,ii)*
!     *                    (1-sca(n,ii))/1000.*aclass(n,ii)*frac(n)*al**2
!                isoLZS2fs(n,l)=isoLZS2fs(n,l)+1.0*drngfs(n,ii)*
!     *                     sca(n,ii)/1000.*aclass(n,ii)*frac(n)*al**2
!              end do

!             INITIALIZE WETLAND ARRAYS
  	        if(wetland_flag(n))then
                isoconcwet(n,l)=0.0
                isoin2wet(n,l)=0.0
                isoout2wet(n,l)=0.0
                isowstore2(n,l)=0.0*wstore2(n)
                isoconcSWwet(n,l)=0.0
                isoin2SWwet(n,l)=0.0
                isoout2SWwet(n,l)=0.0
                isowstore2SW(n,l)=0.0
                isoconcIFwet(n,l)=0.0
                isoin2IFwet(n,l)=0.0
                isoout2IFwet(n,l)=0.0
                isowstore2IF(n,l)=0.0
                isoconcfswet(n,l)=0.0
                isoin2fswet(n,l)=0.0
                isoout2fswet(n,l)=0.0
                isowstore2fs(n,l)=0.0
                isoconcSWfswet(n,l)=0.0
                isoin2SWfswet(n,l)=0.0
                isoout2SWfswet(n,l)=0.0
                isowstore2SWfs(n,l)=0.0
                isoconcIFfswet(n,l)=0.0
                isoin2IFfswet(n,l)=0.0
                isoout2IFfswet(n,l)=0.0
                isowstore2IFfs(n,l)=0.0
              endif
            endif
          end do
        endif


!       WRITE OUTPUT FILE HEADERS
        if(icount.eq.0)write(91,11000)
        if(icount.eq.0)write(90,58000)


        if(iopt.gt.1) print*,'initializations completed'
	  call tracer5(iz,jz,time,t,jan,tdum,jjz,jjzold)

! **********************************************************************
	elseif(itrace.eq.100)then

        if(icount.eq.0)then

!         ALLOCATE THE FLOW-TYPE TRACER ARRAYS:
          allocate(isoin1GW(na,no),isoin2GW(na,no),isoout1GW(na,no),
     *     isoout2GW(na,no),isoconcGW(na,no),
     *     isostore1GW(na,no),isostore2GW(na,no),massin(na),massout(na),
     *     masstore(na),ISOdelta(na),tt(na),stat=iAll)
          if(iAll.ne.0) PAUSE 'Error allocating GW Tracer arrays'

!         ALLOCATE WETLAND ARRAYS IF WETLAND OPTION TURNED ON
          if(wetflg.eq.'y')
     *     allocate(isoin1wet(na,no),isoin2wet(na,no),isoout1wet(na,no),
     *     isoout2wet(na,no),isoconcwet(na,no),isowstore1(na,no),
     *     isowstore2(na,no),wmassin(na),wmassout(na),wmasstore(na),
     *     wISOdelta(na),stat=iAll)
	    if(iAll.ne.0) STOP 'Error allocating Wetland Tracer arrays'

!         ALLOCATE LAKE GW TRACER
          allocate(isolakeGW(noresv),
     *             isolakest(noresv),stat=iAll)

	print*,'100 isolakeGW allocated with',noresv

	    if(iAll.ne.0) STOP 'Error allocating Lake Tracer array'          

!         INITIALIZE PRELIMINARY ITERATION VALUES 
          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            l=nbasin(i,j)
!	      INITIALIZE 'OUT' ARRAYS SO 'IN' ARRAYS SET TO ZERO @t=1
!           BASED ON INITIALIZATION LOOP BELOW WHERE 1'S = 2'S.
!           INITIALIZE STORAGE TO SOME BACKGROUND CONC'N OF TRACER:
            if(l.ne.0)then
	        isoconcGW(n,l)=0.0
              isoin2GW(n,l)=0.0
              isoout2GW(n,l)=0.0
              isostore2GW(n,l)=0.75*store2(n)
!             INITIALIZE WETLAND ARRAYS
  	        if(wetland_flag(n))then
                isoconcwet(n,l)=0.0
                isoin2wet(n,l)=0.0
                isoout2wet(n,l)=0.0
                isowstore2(n,l)=0.75*wstore2(n)
              endif
	      endif
          end do
!        INITILIZE LAKE ARRAYS
          do n=1,noresv
            isolakeGW(n)=0.0
            isolakest(n)=0.0
          end do  
	  endif

!       WRITE OUTPUT FILE HEADERS
        if(icount.eq.0)write(90,56000)
	  if(icount.eq.0)write(91,11000)
	  if(icount.eq.0)write(94,11000)

!       output only GW variables
!        if(no.le.20.and.icount.eq.0)write(92,10000)
!        if(icount.eq.0)write(90,56000)

        call tracerGW(iz,jz,time,t,jan,tdum,jjz,jjzold)


! **********************************************************************
!	elseif(itrace.eq.101)then
!	  call tracer101(iz,jz,time,t,jan,tdum)


! **********************************************************************
	else
        print*,'Ivalid itrace code - exiting run'
	  print*,'Specify a different itrace number in PAR file and re-run'
	  STOP 'Program terminated in tracer.for @101'

	endif

! **********************************************************************



! FORMATS:
 1200 FORMAT(f8.2,',',i8,6(',',f16.4),3(',',f20.6),',',f16.4,
     *2(',',f20.6))
50000 FORMAT('   time   ,',
     *'         basin#,          grid#,        qo2    ,',
     *'        qr     ,        ISOoutQ,        ISOconc,')
52000 FORMAT('    time  ,',
     *<no>('         basin#,          grid#,       qo2     ,      ',
     *'sumQii   ,',<classcount-1>('        ISOoutQ,        ISOconc,')))
53000 FORMAT('   time   ,',
     *<no>('         basin#,          grid#,        qo2    ,',
     *'       qstrm   ,       strloss ,        ISOoutQ,',
     *'        ISOconc,'))
54000 FORMAT('   time   ,',
     *<no>('         basin#,          grid#,        tt     ,',
     *'        q02    ,        qlz    ,      ISOoutQGW,',
     *'      ISOconcGW,       q1      ,      ISOoutQSW,',
     *'      ISOconcSW,         qint  ,         ISOoutQIF,',
     *'      ISOconcIF,        ISOsumQ,'))
55000 FORMAT('      time   ,',
     *<no>('       qo2     ,        q1     ,',
     *'      ISOoutQSW,      ISOconcSW,      q1fs     ,',
     *'    ISOoutQSWfs,    ISOconcSWfs,     qint      ,',
     *'      ISOoutQIF,      ISOconcIF,     qintfs    ,',
     *'    ISOoutQIFfs,    ISOconcIFfs,     qdrng     ,',
     *'      ISOoutQGW,      ISOconcGW,     qdrngfs   ,',
     *'    ISOoutQGWfs,    ISOconcGWfs,'))
56000 FORMAT('      time   ,',
     *<no>('QGW  ,'))
57000 FORMAT('      time,',
     *<no>('QSW  ,','QIF  ,',
     *'QGW  ,'))
58000 FORMAT('      time,',
     *<no>('QSW  ,','QSWfs  ,','QIF  ,',
     *'QIFfs  ,','QGW  ,','QGWfs  ,'))
10000 FORMAT('   time   ,',
     *<no>('         basin#,          grid#,          qo2  ,',
     *'        qlz    ,       strloss ,        ISOoutQ,',
     *'        ISOconc,'))
10100 FORMAT('   time   ,',
     *<no>('         basin#,          grid#,          qlz  ,',
     *'        qowet  ,       qswevap ,        ISOoutQ,',
     *'        ISOconc,'))
11000 FORMAT('    time  ,        IN     ,       OUT     ,',
     *'        DELTA  ,       STORE   ,       SSE     ,')


      icount=1
      
      RETURN
      END SUBROUTINE tracer