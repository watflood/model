      subroutine rt_wetland(div,thr,dtmin,jz,iz,time,date,tdum,n)
     
!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen and Tricia Stadnyk
        
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
     
!     rev. 10.2.08 Nov.  04/17  - NK: New rt_channel & rt-_wetland subroutines      
!     rev. 10.2.71 Nov.  18/19  - NK Bug fixes in wetland & reservoir routing NaNtest

      use area_watflood
      implicit none

      SAVE          ! SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT

      INTEGER :: istate(400,400),rbin,lsta,nnn1,jz,n,lll,l,iz,jjz,&
                  i,j,ii,ic,jj,ios,itracker,unt,ln,n_dt_min,&
                  month_last,ic_max
      integer :: iallocate
      REAL    :: oldwet,convert,kfactor,wth,old_at,aaa,tdum
      REAL(4) :: time,newstore,try1,try2,try3,div,thr,at,dtmin,&
                  wt,atemp,ax,xa,at1,dt_min_n,&
                  hold,whold,strlossvol,qtest
      integer :: frame_no1
!        REAL(4), DIMENSION(:) :: att(10000)
      character*1 :: flowsta,firstpass,flwinitflg
      character*80 :: junk
!        REAL    :: qdlast(6)
      character(14) :: date
      logical  :: printmsg,NaNtest
      
!     needs to fixed for dynamic mem
      DATA month_last/0/
      DATA istate /160000*0/
      DATA firstpass/'y'/
      data printmsg/.true./

      itracker=-1

!     rev. 9.1.16  Apr.  03/02  - Added wetland conditional to select river w/wo wetland
!     rev. 9.1.31  Nov.  13/02  - Fixed the wetland Q to account for wetland area
!     rev. 9.1.33  Dec.  05/02  - Fixed instability in wetland flow    
!     rev. 9.5     Sep.  07/07  - NK: changed wetland/channel routing 

!             WETLAND+CHANNEL ROUTING:
!             dacheck is just there to be able to do wetlands in headwater watersheds
!             and not do it in the larger rivers. This should be replaced by a
!             rivertype with -ve wetland parameters.
!             CHANNEL ROUTING PORTION OF CODE:
!             ADD UPSTREAM CONTRIBUTIONS AND LOCAL INFLOW:
      
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
              if(strlossvol.gt.store1(n)+(qi1(n)+qi2(n))*div)then
                    strloss(n)=0.0
              endif

              qi2(n)=qi2(n)+qstream(n)-strloss(n)
              qin(n)=qi2(n)
      
              qiwet2(n)=qr(n)+qswrain(n)-qswevp(n)+qlz(n)    !-qstream(n)+strloss(n)
     
              qiwetsum(n)=qiwetsum(n)+(qiwet1(n)+qiwet2(n))*div
              lll=next(n)
!c              old=qo1(n)
              oldwet=qowet1(n) 
              hold=1.0e+25
              whold=1.0e+25
              ic_max=1
              
!     rev. 9.9.53  Jan.  18/15  - NK: Prevent mode switch during iteration in wetland routing
!             The type of flow is no longer allowed to change
!             while iterating to prenent instabilities
!             The type of flow is based on conditions for the last time step.
              if(store1(n).le.0.0)then
!               NO OUTFLOW - CHANNEL IS EMPTY
!               outflow from the wetland is added to the channel.
                itracker=0
              endif
              hcha1(n)=store1(n)/chaarea(n)
!             determine what condition exists and keep it                   
              if(hcha1(n).le.chadep(n))then
!               CHANNEL FLOW ONLY
                if(wstore1(n).gt.wcap(n))then
!                 WETLAND IS FULL - OVERFLOWS INTO CHANNEL:
!                 WETLAND -> CHANNEL FLOW ONLY
                  itracker=1
                else
!                 WETLAND & CHANNEL DEPTH BELOW BANKFULL
!                 2 WAY  <---> FLOW
                  itracker=2
                endif    ! (wstore2(n).gt.wcap(n))
              else         ! hcha2(n) > chadep(n)     over > 0.0
!               CHANNEL DEPTH ABOVE BANKFUL
                if(wstore1(n).gt.wcap(n))then
!                 WETLAND & CHANNEL ABOVE BANKFULL
!                 2 WAY  <---> FLOW
                  itracker=3
                else       !  wstore < wcap
!                 CHANNEL ABOVE BANKFULL & WETLAND BELOW
!                 CHANNET --> WETLAND FLOW ONLY
                  itracker=4
                endif    
              endif     
              do ic=1,20   ! convergence loop
                if(hwet2(n).le.0.0)hwet2(n)=qlz(n)*div*2/wetarea(n)/theta(n)
              
!     rev. 9.9.48  Jan.  06/14  - NK: Added wetland cond. function for o/b flow
!               UP TO 20 ITERATIONS ARE ALLOWED
!     rev. 9.8.12  Dec.  07/11  - NK: removed 30 char limit on find filetype 
!               changed because both the channel routing & the wetland routing
!               to converge.
!c                if(abs(whold-wstore2(n)).gt.0.00001*whold.or.
!c     *                    abs(hold-store2(n)).gt.0.0003*hold)then
                if(abs(whold-wstore2(n)).gt.0.00001*whold.or.abs(hold-store2(n)).gt.0.00001*hold)then
!                 THIS ITERATES TO 3% OR ALLOWS UP TO 20 ITERATIONS
!c                  if(store2(n).le.0.0)then
!     rev. 9.9.53  Jan.  18/15  - NK: Prevent mode switch during iteration in wetland routing
                  if(itracker.eq.0)then
!                   NO OUTFLOW - CHANNEL IS EMPTY
!c                    itracker=0
                    flowxa(n)=0.000
                    hcha2(n)=0.000
                    qo2(n)=0.000
!     rev. 9.1.38  Apr.  06/03  - Fixed wetland routing when channel is dry
!                   Added this section to calculate wetland outflow
!                   even if the channel is empty.
!                   Also, made the convergence check to 1 mm in wetheight.
!     rev. 9.5.04  Dec.  27/07  - NK: fixed bug in wetland routing
                    hwet2(n)=wstore2(n)/wetarea(n)/theta(n)
                    qowet2(n)=kcond(n)*(hwet2(n)**2)*astep
!                   using astep makes it independent of grid size
!                   assumes eq. calculates flow/km
!c                  else    ! (store2(n).le.0.0)
                    hcha2(n)=store2(n)/chaarea(n)
!c                    if(hcha2(n).le.chadep(n))then
                  elseif(itracker.eq.1)then
!                     CHANNEL FLOW ONLY
!c                      itracker=1
!c                      if(wstore2(n).gt.wcap(n))then
!                       WETLAND IS FULL - OVERFLOWS INTO CHANNEL:
                        hwet2(n)=chadep(n)+(wstore2(n)-wcap(n))/wetarea(n)
!     rev. 9.9.61  Mar.  06/15  - NK: In route: restored hcha2(n)=store2(n)/chaarea(n)
                        hcha2(n)=store2(n)/chaarea(n)
                        if(abs(hwet2(n)-hcha2(n)).gt.0.0001)then
                        qowet2(n)=kcond(n)/2*(hwet2(n)**2-hcha2(n)**2)*astep
                        else
                          qowet2(n)=0.0000
                        endif
                 elseif(itracker.eq.2)then
!c                        itracker=2
!                       WETLAND IS NOT FULL - NOTHING OVERFLOWS:
!     rev. 9.9.61  Mar.  06/15  - NK: In route: restored hcha2(n)=store2(n)/chaarea(n)
                        hcha2(n)=store2(n)/chaarea(n)
                        hwet2(n)=wstore2(n)/wetarea(n)/theta(n)
!c                        if(abs(hwet2(n)-hcha2(n)).gt.0.0001.and.
!c     *                                         hwet2(n).gt.0.3)then
                        if(abs(hwet2(n)-hcha2(n)).gt.0.0001)then
                          qowet2(n)=kcond(n)/2*(hwet2(n)**2-hcha2(n)**2)*astep
                        else
                          qowet2(n)=0.0000
                        endif
!c                      endif    ! (wstore2(n).gt.wcap(n))
!c                    else         ! hcha2(n) > chadep(n)     over > 0.0
!                     CHANNEL + FLOOD PLAIN FLOW
!c                      if(wstore2(n).gt.wcap(n))then
                  elseif(itracker.eq.3)then
!c                        itracker=3
!                       WETLAND IS FULL - OVERLAND FLOW
                        hwet2(n)=chadep(n)+(wstore2(n)-wcap(n))/wetarea(n)
!     rev. 9.9.61  Mar.  06/15  - NK: In route: restored hcha2(n)=store2(n)/chaarea(n)
                        hcha2(n)=store2(n)/chaarea(n)
!                       this means that if water deprth on the wetland
!                       is 1 m the conductivity is 26 times kcond
                        if(abs(hwet2(n)-hcha2(n)).gt.0.0001)then
                          qowet2(n)=kcond(n)/2*(hwet2(n)**2-hcha2(n)**2)*astep
                        else
                          qowet2(n)=0.0000
                        endif
                  else       !  wstore < wcap
!c                        itracker=4
!                       CHANNEL OVERFLOWS INTO WETLANd:
                        hwet2(n)=wstore2(n)/wetarea(n)/theta(n)
!     rev. 9.9.61  Mar.  06/15  - NK: In route: restored hcha2(n)=store2(n)/chaarea(n)
                         hcha2(n)=store2(n)/chaarea(n)
                        if(abs(hwet2(n)-hcha2(n)).gt.0.0001)then
                          qowet2(n)=kcond(n)/2*(hwet2(n)**2-hcha2(n)**2)*astep
                        else
                          qowet2(n)=0.0000
                        endif
!c                    endif     ! (hcha2(n).le.chadep(n))
                  endif     ! (store2(n).le.0.0)
              
!     rev. 10.2.71 Nov.  18/19  - NK Bug fixes in wetland & reservoir routing
                  NaNtest=ISNAN(qowet2(n))
                  if(NaNtest)qowet2(n)=0.0
                  if(NaNtest)write(98,*)'Warning: wetland NaN in grid ',n,'at t= ',totaltime

                  flowxa(n)=store2(n)/rl(n)
                  qo2(n)=flowxa(n)**1.67*slope(n)/chawid(n)**0.667/r2n(n)*ice_fctr(n)
                  NaNtest=ISNAN(qo2(n))
                  if(NaNtest)qo2(n)=0.0
                  if(NaNtest)write(98,*)'Warning: flow NaN in grid ',n,'at t= ',totaltime

!     rev. 9.8.12  Dec.  07/11  - NK: removed 30 char limit on find filetype 
                  wt=amax1(.1,float(ic)/21.)
                  wt=amin1(1.00,wt)
                    
!c                  endif     ! (store2(n).le.0.0)
                  qi2(n)=qin(n)+qowet2(n)     ! moved from above
!c                  old=qo2(n)
!c                  qold(n)=qo2(n)                  !  not used really
                  oldwet=qowet2(n)
                  hold=store2(n)
                  whold=wstore2(n)
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals

                  store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n)-2.0*qwdr(n,month_now))*div
                  
!     rev. 9.9.52  Jan.  14/15  - NK: Fixed bug for channel store < 0 for withdrawals
                  if(store2(n).le.0.0)then
!                   cut off withdrawals                  
                    store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
                    if(store2(n).le.0.0)then
!                     if store2 still = 0 then reduce calculated outflow
!                     but this should not happen as it never did before 
!                     withdrawals were added to the code
                        
                      store2(n)=0.0
                      qo2(n)=store1(n)/div+qi1(n)+qi2(n)-qo1(n)         
                      if(qo2(n).lt.0.000.and.iopt99)write(63,*)'W',n,id,time,qo2(n)
                    endif
                  endif
!                            qwdr multiplied by 2 to offset /2 in div     
                  wstore2(n)=wstore1(n)+(qiwet1(n)+qiwet2(n)-qowet1(n)-qowet2(n))*div
                  store2(n)=(1.0-wt)*store2(n)+wt*hold
                  wstore2(n)=(1.0-wt)*wstore2(n)+wt*whold
                  if(wstore2(n).lt.0.0)wstore2(n)=0.0
                  satxa(n)=wstore2(n)/rl(n)/theta(n)
                else
!                 CONVERGENCE TO 3%
                  GO TO 26
                endif   ! (abs(whold-wstore2(n)).gt.0.00001*whold.or.
!     *                    abs(hold-store2(n)).gt.0.0003*hold)
              end do    ! ic=1,20
              ic_max=ic
              
26            continue
              if(dds_flag.eq.0)then            ! changed Nov. 27/14  nk
                if(store2(n).le.0.0)then
                  write(98,*)'Warning: In route @ line 486 msg # 1'
                  write(98,*)'Warning: River class',ibn(n), 'Time =',jz,'row',yyy(n),'col',xxx(n)
                  write(98,*)'Warning: store2(',n,' )-ve / needs work in route @ 603'
                  write(98,*)'Warning:', store1(n),qi1(n),qi2(n),qo1(n)
                  if(qwdr(n,month_now).gt.0.0)qwdr(n,month_now)=0.0
                  write(98,*)'Warning: withdrawals set = zero'
                  qo2(n)=store1(n)/div+qi1(n)+qi2(n)-qo1(n)
                  store2(n)=0.0
                  dtmin=a6
                endif      ! (store2(n).le.0.0)
              endif        ! (dds_flag.eq.0)            
              
      if(qo2(n).lt.0.000)write(98,*)'Warning: W',n,id,time,qo2(n)

! >  >  >     MAYBE NEXT LINE HAS TO BE CHECKED OUT
!             WHY IS IT 0 ANYWAYS ??????
!c              if(qo2(n).le.0.0)qo2(n)=0.0000

              qowetsum(n)=qowetsum(n)+(qowet1(n)+qowet2(n))*div

              if(qowet2(n).gt.1.0e+10)then
!     rev. 9.9.01  Dec.  12/13  - NK: Added `pintwarning' in route added
                if(printwarning(n))then
                  write(98,*)'WARNING: '
                  write(98,*)'WARNING: In route @ line 506 msg # 2'
                  write(98,*)'WARNING: Likely fp overflow: reduce kcond for'
                  write(98,*)'WARNING:  class',ibn(n),'kcond=',kcond(n)
                  write(98,*)'WARNING:  grid #',n,' row/col:',yyy(n),xxx(n)
                  if(ic_max.ge.20)then
                    write(98,*)'WARNING: OR most likely:'
                    write(98,*)'WARNING:  min time step A6 too large - try A6/2'
                  endif
                  write(98,*)'WARNING: Message for this grid is printed only once'
                if(writeflg(51))write(51,*)n,ibn(n),kcond(n),itrace,hwet2(n),hcha2(n),qowet2(n),qlz(n)
!c                writeflg='y'
                printwarning(n)=.false.
              endif     ! (printwarning(n))
              if(printmsg)then
                write(98,*)'WARNING: There are ways to correct this problem'
                write(98,*)'WARNING: Possible problem(s): '
                write(98,*)'WARNING: No initial flows for downstream stations'
                write(98,*)'WARNING:  init flow =',qda(n)
                write(98,*)'WARNING: If that is not it,'
                write(98,*)'WARNING:  Check initial lake levels and'
                write(98,*)'WARNING:  reservoir releases. If all ok,'
                write(98,*)'WARNING: You can increase the width/depth ratio and/or '
                write(98,*)'WARNING:  get rid of the wetland in this grid and/or'
                write(98,*)'WARNING:  you can reduce the conductivity and/or' 
                write(98,*)'WARNING:  increase the width/depth ratio in this '
                write(98,*)'WARNING:  and other grids by making a special river '
                write(98,*)'WARNING:  class for problem grids'
                write(98,*)'WARNING: After fixing one grid, others may show up'
                write(98,*)'WARNING:  you may find these grids by animating the '
                write(98,*)'WARNING:  grid outflow in GreenKenue in the log scale'
                write(98,*)'WARNING:  and noting the grids where the flow'
                write(98,*)'WARNING:  disappears'
                write(98,*)'WARNING: OR'
                write(98,*)'WARNING:  maybe min time step A6 too small - try A6/2'
                write(98,*)'WARNING: This message appears only once'
!     rev. 9.8.52  Mar.  20/13  - NK: deleted a pause for dds runs in route
!                 to skip the pause for dds & sensitivity runs	            
!     rev. 9.8.55  Mar.  20/13  - NK: fixed pause for dds runs in route see 9.8.52
!                 had it backwards  nk
!c	            if(abs(dds_flag).eq.0)pause 'Hit enter to continue the run'
                printmsg=.false.
              endif     ! (printmsg)

              endif       ! (qowet2(n).gt.1.0e+10)

              if(n.eq.naa.and.writeflg(51))write(51,*)'         n,          ibn,       kcond,        itrace,      hwet2,       hcha2,    qowet2       qlz(n)'

!             CALCULATE THE VELOCITY THROUGH EACH SQUARE
!             CALCULATE TRAVEL TIME FOR MAXIMUM VELOCITY.

!     rev. 9.1.14  Mar.  24/02  - fixed wetland min time step and outflow
          if(qo2(n).gt.0.000001.and.qowet2(n).gt.0.000001)then
              at=amin1(store2(n)/qo2(n),abs(wstore2(n)/qowet2(n)))

!             SELECT MIN TRAVEL TIME FOR THE TIME STEP CALC
!              >>>dtmin=amin1(at,dtmin,a6)
!             took out the a6 Mar. 26/01 nk
!             don't know why it was put in. Was not in spl8
!              dtmin=amin1(at,dtmin)
!             >>put it back in Oct. 24/01 because it crashed on tabaco creek
!             o2 and store2 do some weird things.

!c              dtmin=amin1(at,dtmin,a6)

!     rev. 9.9.49  Jan.  06/14  - NK: Added courantflg
              if(at.lt.a6)then
                courantflg=.false.
              endif

              dtmin=amin1(at,dtmin)
              dtmin=amax1(dtmin,a6)   ! dtmin > a6 no matter what
          endif
              
              

!             DTMIN IS THE TIME REQUIRED TO COMPLETELY DRAIN
!             THE FASTEST EMPTYING ELEMENT
!             CALCULATE THE CHANNEL STATE FOR GRAPHICAL OUTPUT:
              i=yyy(n)
              j=xxx(n)
              atemp=qo2(n)/(0.4*bnkfll(n))+1.0

!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
              bankfull(i,j)=qo2(n)/bnkfll(n)*100.0
              bankfull(i,j)=amax1(0.1,bankfull(i,j))
              bankfull(i,j)=amin1(10000.,bankfull(i,j))

!             TO PREVENT INTEGER UNDERFLOW OR OVEFLOW:  
              atemp=amax1(atemp,1.0)
              atemp=amin1(atemp,99.0)
              istate(i,j)=int(atemp)
!#ifdef TEST
              if(iopt.ge.4) write(55,1002)i,j,n,istate(i,j),bnkfll(n),qo2(n)
 1002 format(' i,j,n,istate,bnkfull,qo2/',4i5,2f10.3)
              if(iopt.ge.1)then
!               if(n.eq.nnprint.and.iz.ne.jz)
                if(n.eq.nnprint)then
!c                  if(iz.ne.jz)
                    call write_wetland(itracker,ic,n,totaltime,68)
                endif
                endif
!#endif

      return
      
      end subroutine rt_wetland
!     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
              
              
      subroutine write_wetland(itracker,ic,n,time,unt)

      use area_watflood
      implicit none

        INTEGER :: ii,ic,ios,itracker,unt,n,iz,jz
        REAL(4) :: time
        logical :: firstpass
        
        data firstpass/.true./
        
        if(firstpass)then
          if(writeflg(unt))write(unt,55998)'itracker,ic,totaltime,&
         store1(n),store2(n),qi1(n),qi2(n),qo1(n),qo2(n),&
         wstore1(n),wstore2(n),&
         hcha1(n),hcha2(n),hwet1(n),hwet2(n),&
         flowxa(n),chaxa(n),wetxa(n),satxa(n),&
         qiwet1(n),qiwet2(n),qowet1(n),qowet2(n),&
         qstream(n),qswrain(n),strloss(n),&
         qlz(n),qswevp(n),&
         kcond(n),theta(n),chaxa(n)/chawid(n),over(n)*rl(n)'
55998     format(a340)     
          
        endif

!       write(68 ..... for tracking purposes
        if(writeflg(unt))write(unt,55999)itracker,ic,totaltime,&
         store1(n),store2(n),qi1(n),qi2(n),qo1(n),qo2(n),&
         wstore1(n),wstore2(n),&
         hcha1(n),hcha2(n),hwet1(n),hwet2(n),&
         flowxa(n),chaxa(n),wetxa(n),satxa(n),&
         qiwet1(n),qiwet2(n),qowet1(n),qowet2(n),&
         qstream(n),qswrain(n),strloss(n),&
         qlz(n),qswevp(n),&
         kcond(n),theta(n),chaxa(n)/chawid(n),over(n)*rl(n)

55999 format(2(i5,','),f10.3,',',16(f12.3,','),9(f12.5,','),2(f5.2,','),5x,2(f12.3','))
  
      firstpass=.false.
  
      end subroutine write_wetland
