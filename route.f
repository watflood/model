      SUBROUTINE route(div,thr,dtmin,jz,iz,time,date,tdum)

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen 
        
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
!
!  THE INFLOW INTO THE CHANNEL IN EACH ELEMENT IS THE SUM OF THE
!  OUTFLOW OF THE ELEMENTS WHICH DRAIN INTO IT FROM ABOVE AND THE
!  SURFACE RUNOFF AND SUBSURFACE FLOW FROM THE ELEMENT IN WHICH THE
!  OUTFLOW IS BEING CALCULATED                                 
!
!     REV. 8.23 - Mar.  25/96 -  fixed bug in route - keep qo2 for res
!     rev. 8.99mm Dec. 13/2001-     added check for <= 0 init res flow
!     rev. 8.99n  Dec. 31/2001-     fixed nat. res initial flow (JW)
!
!     REV. 9.00    Mar.  2000 - TS: CONVERTED TO FORTRAN 90
!     REV. 9.03    Nov.  2000 - TS: ADDED WATFLOOD SWAMP ROUTING 
!     rev  9.1.03  July  24/01  - added polinomial to reservoir routing
!     rev. 9.1.10  Jan.  29/02  - flow nudging added for nopt(l)=2
!     rev. 9.1.14  Mar.  24/02  - fixed wetland min time step and outflow
!     rev. 9.1.16  Apr.  03/02  - Added wetland conditional to select river w/wo wetland
!     rev. 9.1.31  Nov.  13/02  - Fixed the wetland Q to account for wetland area
!     rev. 9.1.33  Dec.  05/02  - Fixed instability in wetland flow    
!     rev. 9.1.38  Mar.  31/03  - revised str header and routing dt selectable
!     rev. 9.1.39  Apr.  06/03  - Fixed wetland routing when channel is dry
!     rev. 9.2.11  Sep.  11/05  - NK: added Manning's n  r1n & r2n
!     rev. 9.2.13  Sep.  28/05  - NK: added freeze and break up to route
!     rev. 9.2.23  Nov.  22/05  - NK: Fixed res(n)=0 bug in route 
!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!     rev. 9.2.39  May.  09/06  - NK: t added to route & rerout arg listrules
      
!     rev. 9.2.43  Jun.  21/06  - NK: fixed spikes in route
!     rev. 9.3.04  Oct.  24/06  - NK: routing parameters dim to na in rte
!     rev. 9.3.10  Jan.  29/07  - NK: routing pars changed to gridded values
!                                 eg: lzf(ii) -> lzf(n) 
!     rev. 9.5     Sep.  07/07  - NK: changed wetland/channel routing 
!     rev. 9.5.01  Oct.  15/07  - NK: added wetland continuity check
!     rev. 9.5.02  Oct.  21/07  - NK: set init qdwpr=0.0 in route
!     rev. 9.5.04  Dec.  27/07  - NK: fixed bug in wetland routing
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
!     rev. 9.5.21  Mar.  06/08  - NK: fixed dtmin for first time step each event
!     rev. 9.5.25  Mar.  20/08  - NK: fixed lake initiation - moved code route -> flowinit
!     rev. 9.5.39  Oct.  15/08  - NK: fixed bug in reservoit routing
!     rev. 9.5.46  Dec.  23/08  - NK: trying to fix problem with -ve storage. Changed conditional to .lt.
!     rev. 9.5.47  Dec.  26/08  - NK: add flwinitflg to warn about initial flows
!     rev. 9.5.54  Feb.  11/09  - NK: undid rev. 9.2.28
!     rev. 9.5.58  Apr.  16/09  - NK: added nudgeflg for forcing gauge flows
!     rev. 9.5.83  Feb.  17/10  - NK: non_basin exclusion for dds_flag=1
!     rev. 9.8.10  Dec.  06/11  - NK: Added message for FP overflow in route
!     rev. 9.8.12  Dec.  07/11  - NK: removed 30 char limit on find filetype 
!     rev. 9.8.29  Oct.  15/12  - NK: added wetland_flag to speed up route.f
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals
!     rev. 9.8.52  Mar.  20/13  - NK: deleted a pause for dds runs in route
!                                 qwdr = point withdrawals
!                                 taken from store2() in river routing
!                                 taken from lake inflows in lake routing
!
!     changes made to include c&g model stuff  nk  April. 26/07
!
!***********************************************************************

      use area_watflood
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      INTEGER :: istate(400,400),rbin,lsta,nnn1,jz,n,lll,l,iz,jjz,
     *             i,j,ii,ic,jj,ios,itracker,unt,ln,n_dt_min,
     *             month_last,ic_max
      integer :: iallocate
      REAL    :: old,oldwet,convert,kfactor,wth,old_at,aaa,tdum
      REAL(4) :: time,newstore,try1,try2,try3,div,thr,at,dtmin,
     *             wt,atemp,ax,xa,at1,dt_min_n,
     *             hold,whold,strlossvol,qtest
      integer :: frame_no1
!        REAL(4), DIMENSION(:) :: att(10000)
      character*1 :: flowsta,firstpass,flwinitflg
      character*80 :: junk
!        REAL    :: qdlast(6)
      character(14) :: date
      logical  :: printmsg
!     rev. 10.2.35 Oct.  08/18  - NK: Moved logical def. to area_watflood
c      logical, dimension(:), allocatable  :: printwarning
      
!     needs to fixed for dynamic mem
      DATA month_last/0/
      DATA istate /160000*0/
      DATA firstpass/'y'/
      data printmsg/.true./
     
!     rev. 9.5.25  Mar.  20/08  - NK: fixed lake initiation - moved code route -> flowinit
!     rev. 9.5.25  Mar.  20/08  - NK: fixed lake initiation - moved code route -> flowinit
      index=2

      if(firstpass.eq.'y')then   !section added Apr. 28/06  nk
!       check that qdwpr memory has been allocated
!       `firstpass` is passed through all routing routines          
        do n=1,naa
          if(ireach(n).gt.0.or.res(n).gt.0)then
            if(.not.allocated(qdwpr))then
              print*
              print*,'Memory not allocated for qdwpr'
              print*,'No of reservoirs in the .rel file is 0'
              print*,'but reaches in the map & shed files have been'
              print*,'defined. This is a problem.'
              print*,'Please either set all reach values = 0'
              print*,'or code in the reservoir locations in the rel'
              print*,'files'
              stop 'Program aborted in route @ 279'
            endif      
          endif
        end do

        allocate(printwarning(na),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *     'Error with allocation of printwarning in route @ 114'
        do n=1,naa
          printwarning(n)=.true.
          netoutflow(n)=0.0
        end do
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
        if(iopt99)then
          do n=1,naa
            qqlow(n)=10000.0
          end do
        endif
d        print*,'printwarning set to "T"'
        
        
! moved to spl  Jul. 08/13
c!     rev. 9.8.29  Oct.  15/12  - NK: added wetland_flag to speed up route.f
c!       wetland_flag(n) replaced multiple logical checks by just 1 in each
c!       time & grid loop - so will save a lot of time
c        allocate(wetland_flag(naa),stat=iAllocate)
c        if(iAllocate.ne.0) STOP
c     *     'Error with allocation of wetland_flag in route @ 112'
c        do n=1,naa
c          if(aclass(n,classcount-2).gt.0.0
c     *            .and.grid_area(n).gt.10.0   ! added Oct. 15/12 NK
c     *            .and.wetflg.eq.'y'
c     *            .and.theta(ibn(n)).gt.0.00001
c     *            .and.glacier_flag(n).ne.'y')then
c            wetland_flag(n)=.true.
c          else
c            wetland_flag(n)=.false.
c          endif
c	  end do  
c	  
	  flwinitflg='n'
	  old_at=dtmin
	  frame_no1=0
        qlow=0.000000
        
	endif

!     CALCULATIONS START IN THE HIGHEST ELEMENT IN THE WATERSHED 
!     AND PROCEED TO THE LOWEST.

!#ifdef TEST
      if(iopt.ge.3)then
        write(55,6002)iz,jz
        write(53,6002)iz,jz
        write(55,6003)
      endif
!#endif
 
!     ROUTING LOOP:

      jjz=jz
      if(jjz.lt.1) jjz=1
      dt_min_n=1.0e32
      
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

d      if(iopt.eq.2)print*, ' in route before the 1-naa loop'

c!DEC$ NOPARALLEL

      do n=1,naa
!       REV. 8.23 - Mar.  25/96 - FIXED BUG IN ROUTE, KEEP QO2 FOR RES

        store1(n)=store2(n)
        wstore1(n)=wstore2(n)
        hwet1(n)=hwet2(n)
        hcha1(n)=hcha2(n)
        qo1(n)=qo2(n)
        if(res(n).eq.0)qo2(n)=0.0
        qowet1(n)=qowet2(n)
        qi1(n)=qi2(n)
c        qi2(n)=1.0e-10
        qi2(n)=0.0
c        netoutflow(n)=0.0
        
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
        strlossvol=strloss(n)*div*2.0

!     REV. 10.1.24 Jan.  30/16  - NK: Added qUS1 & qUS2 for watbal
        qUS1(n)=qUS2(n)
c        qUS1(n)=0.0
        qUS2(n)=0.0
     
        qiwet1(n)=qiwet2(n)
        qiwet2(n)=1.0e-10

!     rev. 9.5.02  Oct.  21/07  - NK: set init qdwpr=0.0 in route
        if(ireach(n).gt.0)then
          rbin=ireach(n)
!          grid is part of a reservoir or lake 
          qdwpr(rbin,jjz)=0.0
        endif
      end do
      
      do n=1,naa
d        if(iopt.eq.2.and.n.eq.1)print*,'In route, passed 101'

        i=yyy(n)
        j=xxx(n)
        l=nhyd(i,j)

c        if(l.eq.0.and.numa.gt.0)then            ! NON-BASIN EXCLUSION
!     rev. 9.5.83  Feb.  17/10  - NK: non_basin exclusion for dds_flag=1
c         if(l.eq.0.and.dds_flag.eq.1)then            ! NON-BASIN EXCLUSION
c        if(l.eq.0)then            ! NON-BASIN EXCLUSION


!           DON'T DO ANYTHING - WE'RE NOT IN A GRID THAT COUNTS FOR OPT

c        else                    ! (l.eq.0.and.numa.gt.0)

!         WHEN THE SLOPE <= 0.0 THE ELEMENT IS NOT IN THE BASINfat=

!         AND THE ROUTING SEQUENCE IS SKIPPED

          if(slope(n).gt.0.0)then
!           REV. 7.2 Sept. 19/94 - ADDED IREACH() FOR DWOPER INPUT  

! * * * * * * * * * * LAKE or RESERVOIR ROUTING * * * * * * * * * * * * * * * * 
            if(ireach(n).gt.0.or.res(n).gt.0)then   ! res(n).ne.0 outlet
! * * * * * * * * * * RESERVOIR ROUTING * * * * * * * * * * * * * * * * *
! * * * * * * * * * * RESERVOIR ROUTING * * * * * * * * * * * * * * * * *
! * * * * * * * * * * RESERVOIR ROUTING * * * * * * * * * * * * * * * * *
! * * * * * * * * * * RESERVOIR ROUTING * * * * * * * * * * * * * * * * *
! * * * * * * * * * * RESERVOIR ROUTING * * * * * * * * * * * * * * * * *
!             IF WE ARE ROUTING EXTERNALLY WITH DWOPER OR DOING IT IN A 
!             CONTROLLED RESERVOIR DOWNSTREAM IF THERE IS ONE
!             THE FLOWS ARE ACCUMULATED IN A REACH-BIN FOR DWOPER 
!             SEVERAL ELEMENTS CAN CONTRIBUTE BUT NONE IS ROUTED TO 
!             DOWNSTREAM IN ed :

!             ADD UPSTREAM CONTRIBUTIONS AND LOCAL INFLOW
                
              lll=next(n)
              rbin=ireach(n)
c     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso model
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals
	        resvstore1(n)=store2(n)
              
              if(res(n).eq.0)then
!               grid is part of a reservoir or lake but is not an outlet 
c                if(trcflg.ne.'y')then
                  qdwpr(rbin,jjz)=qdwpr(rbin,jjz)+qi2(n)+qr(n)
     *                     +qstream(n)-strloss(n)-qwdr(n,month_now)   !+qlz(n): TS - Mar12/08: included in qr(n)

              else     ! res(n).ne.0     = outlet

!               THERE IS A DAM IN THIS GRID AND THE 
!               RESERVOIR ROUTING SUBROUTINE REROUT IS CALLED.

!               * * * FOR RESERVOIR ROUTING:* * *
d                if(iopt.ge.3)write(55,6004)n,res(n),index

!               IT'S ASSUMED THAT THE DAM IS LOCATED IN A DWOPER
!               REACH NUMBER

!     rev. 9.1.43  Jun.  01/03  - Fixed the qdwpr.txt function - re: last grid in lake
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipla & irrigation withdrawals
!               corrected this June 1/03
!               previously, the last grid was not added to the qdwpr.txt file
!               this next line copied from above
                if(rbin.gt.0)then     ! rbin=reach(n)
!                 there may be a dam but the grid may not have been 
!                 designated as part of a lake
!                 i.e. rbin - the lake number
!                 qdwoper collects all the inflow to the reach=lake=reservoir      
                  qdwpr(rbin,jjz)=qdwpr(rbin,jjz)+qi2(n)+qr(n)
     *                         +qstream(n)-strloss(n)-qwdr(n,month_now)  !+qlz(n): TS - Mar12/08: included in qr(n)
!                 qi2(n) is used in rerouta as the lake inflow for routing
                  qi2(n)=qdwpr(rbin,jjz)
                endif

!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!               note: for external routing b1 has to be -ve
!               so when b1 = 0, releases will be used.
c                if(routeflg.eq.'q'.and.b1(res(n)).lt.-0.001)then
                if(routeflg.eq.'q'.and.b1(res(n)).eq.0.00)then
!                 routing will be done by external routing             
!                 no routing in rerout
!                 this way any lakes where b1 <> 0 will be routed in rerout
!                 either by a rule or release table
                  qo2(n)=0.0   
                else
	          
!     rev. 9.9.65  Apr.  `3/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
!     rev. 10.2.13 Jan.  31/18  - NK: Re-wrote rules.f to mimic stop log operations -> rules_sl.f
                  if(ruleflg.and.resruleflg(res(n)))then
                    if(ruletype(1:7).eq.'StopLog')then  
!                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      call rules_sl(n,div,thr,res(n),
     *                        jz,at,dtmin,date,time,firstpass)
!                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    elseif(ruletype(1:11).eq.'TargetLevel')then
!                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                      call rules_tl(n,div,thr,res(n),
     *                        jz,at,dtmin,date,time,firstpass)
!                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    else
                      print*,'Error'
                      print*,'StopLog or TargetLevel key words'
                      print*,'not found in the rules.ts5 file'
                      print*,'Keyword found =',ruletype 
                      stop 'Program aborted in route @ 321'
                    endif
                  else
      
!                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    call rerout(n,div,thr,res(n),
     *                        jz,at,dtmin,date,time,firstpass)
!                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
                  endif   ! rules
                  
!     rev. 9.1.67  Oct.  21/04  - NK; added unit 80 for lake_stor & lake_outflow
!                  May   02/05  - TS: revised
!                 for flowsd.csv
                  ln=res(n)
!     rev. 10.1.60 Jan.  03/17  - NK: Fixed conditional in route
                  if(jz.ge.1.and.res(n).gt.0)then      ! added conditional Rev. 9.2.43  nk
                    lake_stor(ln,jz)=store2(n)
                    lake_outflow(ln,jz)=qo2(n)
                    if(n.eq.7713)write(334,*)jz,ln,totaltime,qo2(n)
                    net_lake_outflow(ln,jz)=qo2(n)
!     rev. 10.1.96 Sep   11/17  - NK: Added variable lake depth calculation lake_elv()-LKinvert()
                    LKdepth(ln)=lake_elv(ln,jz)-LKinvert(ln)
                    del_stor(ln,jz)=(qi2(n)-qo2(n))*div    !fixed feb. 15/08 -nk- 
                    
c                    if(ln.gt.1)net_lake_outflow(ln,jz)=
c     *                        qo2(n)-lake_outflow(ln-1,jz)
                    
!                    div added May 9/06 nk
                  endif    ! jz
                endif    ! routeflg

                  
                  
              endif    ! end reservoir routing
                  
! * * * * * * * * * *  TS - WETLAND ROUTING  * * * * * * * * * * * * * * * * 
! * * * * * * * * * *  TS - WETLAND ROUTING  * * * * * * * * * * * * * * * * 

!     rev. 9.8.29  Oct.  15/12  - NK: added wetland_flag to speed up route.f
!     rev. 10.2.08 Nov.  04/17  - NK: New rt_channel & rt-_wetland subroutines      
            elseif(wetland_flag(n))then    ! see conditionals above!!!
                
              call rt_wetland(div,thr,dtmin,jz,iz,time,date,tdum,n) 
              lll=next(n)

            elseif(pondflg(n))then

! * * * * * * * * * * * * * * * POND ROUTING * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * POND ROUTING * * * * * * * * * * * * *
                
c              if(.not.netCDFflg)then
c                  IF(iopt.ge.1.and.n.eq.nnprint)THEN
c                      write(554,5550)id,time,at/3600.,qi1(n),qi2(n),
c     *                   qo1(n),qo2(n),store1(n),store2(n),
c     *                   qr(n),qstream(n),strloss(n),ice_fctr(n)
c                  endif
c              endif
                !     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing      
!     rev. 10.2.08 Nov.  04/17  - NK: New rt_channel & rt-_wetland subroutines      
!             For grids with large water areas but not included in lakes or reservoirs, 
!             use pond routing   
              call rt_pond(div,thr,dtmin,jz,iz,time,date,tdum,n) 
              lll=next(n)
c              if(.not.netCDFflg)then
c                  IF(iopt.ge.1.and.n.eq.nnprint)THEN
c                      write(554,5550)id,time,at/3600.,qi1(n),qi2(n),
c     *                   qo1(n),qo2(n),store1(n),store2(n),
c     *                   qr(n),qstream(n),strloss(n),ice_fctr(n)
c                  endif
c              endif
            else     

! * * * * * * * * * * * * * * * CHANNEL ROUTING * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * CHANNEL ROUTING * * * * * * * * * * * * *

                
c        IF(iopt.ge.1.and.n.eq.nnprint)THEN
c          write(555,5550)id,time,at/3600.,qi1(n),qi2(n),
c     *                   qo1(n),qo2(n),store1(n),store2(n),
c     *                   qr(n),qstream(n),strloss(n),ice_fctr(n)
c        endif
                
                
                
                
                
!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing      
!             For grids with large water areas but not included in lakes or reservoirs, 
!             use pond routing   
              call rt_channel(div,thr,dtmin,jz,iz,time,date,tdum,n) 
              lll=next(n)
                  
            endif                         ! END OF ROUTING

!     rev. 9.1.10  Jan.  29/02  - flow nudging added for nopt(l)=2
!           NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING 
!           NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING 
!           NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING NUDGING 
!           for nudgeflg = 2
!           keywords:  nudge  nudgeflg  
            flowsta=' '
            do l=1,no
!             check to see if this grid has a flow station
              if(iflowgrid(l).eq.n.and.nopt(l).eq.2)then
!               we are at a flow station that is to be used for nudging
!               also,nudge only if there is observed flow 
!     REV. 10.1.28 Apr.  26/16  - NK: Fixed first day of output for master_inflows file 
c                if(jz.ge.kt.and.qhyd(l,jz).ge.0.0)then
                if(qhyd(l,jz).ge.0.0)then
                  flowsta='y'
                  lsta=l
                endif
              endif
            end do

!     rev. 9.7.16  Jan.  05/11  - NK: Fixed init flows outside sub-basin
            if(lll.gt.0)then   ! added Jan/11 nk
              if(flowsta.eq.'y')then  
!               we are in a grid with a flow station
!               used instead of the computed flow     
!               ADD observed FLOW TO DOWNSTREAM ELEMENT
                qi2(lll)=qi2(lll)+qhyd(lsta,jz)
!       if(lsta.eq.21)print*,'xxyyzz',jz,lll,qhyd(lsta,jz),qi2(lll)
!               Note: qo2(n) is NOT changed so the plots will show the computed 
!                     hydrograph. Also, error is based on computed flow.                
              ELSE
!               ADD computed FLOW TO DOWNSTREAM ELEMENT
!               but only if it is not a lake
               if(nbsflg.eq.'n')then
!                  Normal routing = add outflow to downstream grid                    
                   qi2(lll)=qi2(lll)+qo2(n)
               else
!                 nbsflg = y or 0 so don't route downstream if at lake outlet 
!                 nbsflg = 0 no routing AND no data written to NBS output file
                  if(res(n).eq.0)then  ! no reservoir outlet in this grid
!                   nbsflg = y means that the NBS file will be written
!                   nbsflg = 0 means the outflows are set to 0.0 but 
!                   flows will not be written      
                    qi2(lll)=qi2(lll)+qo2(n)
                  endif
                endif
              endif
            endif

!#ifdef TEST
            if(iopt.ge.3)then
!             this creates huge files
              at1=at/3600.0
              write(55,6037,iostat=ios)
     *             n,yyy(n),xxx(n),ic,at1,qr(n),qi1(n),qi2(n),
     *        qo1(n),qo2(n),qi2(lll),store1(n),store2(n),cap(n),lll
!              if(ios.ne.0)then
!                print*,'n,lll,na/',n,lll,na
!                stop
!              endif
            endif
!#endif
          endif                            ! SLOPE IF
c        endif           ! NON-BASIN EXCLUSION l.eq.0.and.dds_flag.eq.1 

!#ifdef TEST
d      if(iopt.eq.2.and.n.eq.naa)print*,'In route, passed 901'
!#endif

!       take out
!       print*,time,at,dtmin
        att(n)=at/3600.
!     added for water balance calculations

        if(frac(n).gt.0.0.and.aclass(n,ii_water).gt.0.0)then
           uzs(n,ii_water)=uzs(n,ii_water)
     *          +(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
     *          /frac(n)/aclass(n,ii_water) ! took out convert 14/04/08 nk
!     *          /convert/frac(n)/aclass(n,ii_water) 
        endif

        if(firstpass.eq.'y')then
!         write the header in rte.txt
          IF(iopt.ge.1.and.n.eq.nnprint)then
          write(55,*)'Routing state variables for debug grid # ',nnprint
          write(55,*)'Note time interval may change - i.e. not regular'
          write(55,5551)
5551      format('   id   time     at(hrs)        qi1         qi2'
     *         '         qo1         qo2       store1      store2'
     *         '          qr       qstream     strloss ice_fctr')
          endif
	  endif

!             added Dec. 12/00 nk.
        IF(iopt.ge.1.and.n.eq.nnprint)THEN
          write(55,5550)id,time,at/3600.,qi1(n),qi2(n),
     *                   qo1(n),qo2(n),store1(n),store2(n),
     *                   qr(n),qstream(n),strloss(n),ice_fctr(n)
5550      format(i5,f10.3,11e12.3)
        endif
        
!     REV. 10.1.24 Jan.  30/16  - NK: Added qUS1 & qUS2 for watbal
c        qUS2(next(n))=qUS2(next(n))+qo2(n)
        
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
        if(iopt99)then
          if(qo2(n).lt.0.000)then
            if(.not.warning(1))then
                 write(63,*)'Reporting flows lower than ',qlow,' = qlow'
                  write(63,*)'mode:  c= channel, w=wetland routing'
                   write(63,*)
     *            'mode     grid#          id       time       outflow'      
            endif
            write(63,*)'C',n,id,time,qo2(n)
            if(qo2(n).lt.qlow)then
              qlow=qo2(n)  ! global low flow for warning
              lowtime=time
              lowid=id
              lowgrid=n
              warning(1)=.true.
            endif
          endif
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
          if(ireach(n).eq.0)then
!           for rivers only - not lakes or reservoirs                  
            qqlow(n)=amin1(qqlow(n),qo2(n))
          endif
        endif
        
      end do    ! n=1,naa
        
! END NAA LOOP   END NAA LOOP   END NAA LOOP   END NAA LOOP   END NAA LOOP   END NAA LOOP   END NAA LOOP        
        
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      if(iopt99)then
        if(time.ge.mhtot)then
!         File replaced at the end of each event              
          open(unit=99,file='debug\qlow.xyz',status='unknown')
          do n=1,naa
              i=yyy(n)
              j=xxx(n)
              write(99,*)xorigin+float(j-1)*xdelta+xdelta/2.0,
     *               yorigin+float(i-1)*ydelta+ydelta/2.0,qqlow(n)
          end do
          close(unit=99,status='keep')
          write(98,*)'Info: debug\qlow.xyz written & closed'
        endif
      endif
!####################################################################################      
        
!     rev. 9.9.73  Aug.  31/15  - NK: Finshed rules s/r - ready for beta testing
!     nrules is the number of reservoir with water level target rules

c      if(ruleflg)write(987,98700)totaltime,
c     *             (lower_range(jul_day_now,i),
c     *              upper_range(jul_day_now,i),
c     *              lake_elv(resvNo(i),jz),i=1,nrules)
c98700 format(f10.0,<nrules*3>(',',f10.3)) 

!#ifdef TEST
d      if(iopt.eq.2)print*, ' in route after the 1-naa loop'
!#endif

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!#ifdef TEST
d      if(iopt.ge.2.and.iopt.le.10) write(55,6006)dtmin
d      if(iopt.ge.3.and.iopt.le.10) write(53,1001)jz
!#endif




!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
c        write(56,*)'time= ',totaltime
c        if(iz.ne.jz)then
c            do i=1,ycount
cc              write(56,56000)(bankfull(i,j),j=1,xcount)
c56000         format(<xcount>f8.0)              
c            end do
c        endif
        
        
        
        

d      if(iopt.eq.2)print*,'finished writing istate'

!     INFORMATION FOR WATBAL.FOR
!     change in volume of the channel or lake
!     REV. 10.1.24 Jan.  30/16  - NK: Added qUS1 & qUS2 for watbal
!     qUS1 and qUS2 are the channel inflow from upstream grid(s)
      do n=1,naa
        if(ireach(n).gt.0.and.res(n).eq.0)then
!         in a lake or reservoir        
          netoutflow(n)=0.0
        elseif(ireach(n).gt.0.and.res(n).ne.0)then
!         in the outlet grid of a lake or reservoir
          netoutflow(n)=netoutflow(n)
c     *	     +(qo1(n)+qo2(n)-qUS1(n)-qUS2(n))/tdum/2.0
     *	                +(qo2(n)-qUS2(n))/tdum   !checked out
        else
!         qUS1 & qUS2 is just the inflow from upstream        
          netoutflow(n)=netoutflow(n)
c     *	     +(qo1(n)+qo2(n)-qUS1(n)-qUS2(n))/tdum/2.0
     *	                +(qo2(n)-qUS2(n))/tdum   !checked out
        endif
        end do
        
!     rev. 9.8.80  Aug   09/13  - NK: Added withdraw.r2c output file in route.f
!     write the monthly withdrawals
      if(month_now.ne.month_last.and.iopt99)then
        xcount_temp=xcount
        ycount_temp=ycount
        do i=1,ycount
          do j=1,xcount
            outarray(i,j)=0.0
          end do
        end do
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          outarray(i,j)=qwdr(n,month_now)
        end do  
        frame_no1=frame_no1+1
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(23,23,ni*12,0,frame_no1,0,11)    
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  month_last=month_now  

      endif

!     rev. 9.5.47  Dec.  26/08  - NK: add flwinitflg to warn about initial flows
!     rev. 9.8.10  Dec.  06/11  - NK: Added message for FP overflow in route
      if(firstpass.eq.'y'.and.flwinitflg.eq.'y')then
	  print*,'WARNING:'
	  print*,'It looks like your initial flows are too low '
	  print*,'or none were given in the first event.'
	  Print*,'Please enter some appropriate value in the first'
	  print*,'line of the str file for each station.'
	  print*,'If you started the run in mid-winter, '
        print*,'that could be your problem - recorded flows tend to'
	  print*,'very small in cold climates.'
        flwinitflg='n'
c	  pause 'hit enter to continue but please fix 1st str.r2c file'
	endif

      old_at=at
      firstpass='n'

      RETURN

!#ifdef TEST
d      if(iopt.eq.2)print*, ' assigned netflow'
!#endif

! FORMATS


 1000 format(50i3)
 1001 format(3i5)
 6000 format(60i2)
 6001 format(' ',3i5,3f10.2)
 6002 format(' route:iz,jz/',2i5)
 6004 format(' gone to rerout - n,res(n),index/',3i5)
 6005 format(' n,res(n),jz,qo2(n)/',3i5,2f10.2)
 6006 format(' dtmin =',f10.2)
 6007 format(' ','reservoir locations wrt 1-naa')
 6037 format(4i5,7f9.3,3f15.0,i5)
 6003 format('    n    i    j   ic     at       qr      qi1',
     *'     qi2      qo1      qo2  qi2[lll]   store1    store2    cap   
     *lll')
 9801 format(a80)
 9802 format(i5,10g12.3)


      END SUBROUTINE route

