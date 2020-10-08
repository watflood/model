      subroutine rt_channel(div,thr,dtmin,jz,iz,time,date,tdum,n)
      
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
     
!     rev. 10.2.08 Nov.  04/17  - NK: New rt_channel & rt-_wetland subroutines      
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
      
!             CHANNEL ROUTING:
!             ADD UPSTREAM CONTRIBUTIONS AND LOCAL INFLOW:
!              if(aclass(n,classcount-2).gt.0.0)

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
      logical  :: printmsg
      
!     needs to fixed for dynamic mem
      DATA month_last/0/
      DATA istate /160000*0/
      DATA firstpass/'y'/
      data printmsg/.true./
      
      
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
                if(strlossvol*2.0.gt.store1(n)+(qi1(n)+qi2(n))*div)then
                    write(98,98000)'Warning:',time,'strloss ',strloss(n),' set = 0.0 for grid ',n
98000               format(a10,i10,a10,f10.3,a20,i10)                    
                    strloss(n)=0.0
                    warning(2)=.true.
                endif
!####################################################################################      

!     rev. 9.9.59  Feb.  18/15  - NK: In route: strloss option fracflg y/n
!     rev. 10.1.62 Jan.  08/17  - NK: Checkup on strloss effect on low flows
      qi2(n)=qi2(n)+qr(n)+qstream(n)-strloss(n)
        
      lll=next(n)
      qold(n)=qo1(n)
      hold=1.0e+25
      do ic=1,20
!       UP TO 20 ITERATIONS ARE ALLOWED
        if(abs(hold-store2(n)).gt.0.0003*hold)then
!         THIS ITERATES TO 3% OR ALLOWS UP TO 15 ITERATIONS
          if(store2(n).le.0.0)then
!           NO OUTFLOW - CHANNEL IS EMPTY
            ax=0.0
            qo2(n)=0.001
          else     ! (store2(n).le.0.0)
            over(n)=(store2(n)-cap(n))/rl(n)
            if(over(n).le.0.0)then
!             CHANNEL FLOW ONLY
              ax=store2(n)/rl(n)

!     rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
!             pool-riffle storage added
!             pool = % of bankfull area attributed to dead storage
            if(ax/chaxa(n).gt.pool(n))then
                qo2(n)=(ax-pool(n)*chaxa(n))**1.67*slope(n)/chawid(n)**0.667/r2n(n)*ice_fctr(n)
              else
                qo2(n)=0.0
            endif
            hcha2(n)=ax/chawid(n)
              
            else     ! (over(n).le.0.0)

!             CHANNEL + FLOOD PLAIN FLOW
!     rev. 9.2.43  Jun.  21/06  - NK: fixed spikes in route
              ax=cap(n)/rl(n)   ! added Jun 21/06 nk

!     rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
!     rev. 9.8.60  May   14/13  - NK: fixed ice factor for whole x-section
!             0.17 factor is based on 100:1 fp w/d ratio
!             flood plain width/depth assumes as 100
!             use quadratic equation to solve for fp. depth
              hwet2(n)=(-1.0+sqrt(1.+400.0*over(n)))/200.0
!             hcha2(n) is the bankfull depth here
              hcha2(n)=chaxa(n)/chawid(n)
!             xa is the total main channel cross section area
              xa=(hwet2(n)+hcha2(n))*chawid(n)
              if(over(n)-hwet2(n)*chawid(n).gt.0.0)then
                qo2(n)=&
                   (xa**1.67*slope(n)/chawid(n)**0.667/r2n(n)&
                   +(over(n)-hwet2(n)*chawid(n))**1.33&
                            *slope(n)*0.17/r1n(n))*ice_fctr(n)
              else
                qo2(n)=xa**1.67*slope(n)/chawid(n)**0.667/r2n(n)*ice_fctr(n)
              endif     ! (over(n)-hwet2(n)*chawid(n).gt.0.0)
            endif      ! (over(n).le.0.0)
!     rev. 9.8.12  Dec.  07/11  - NK: removed 30 char limit on find filetype 
!c            wt=amax1(.5,float(ic)/21.)
            wt=amax1(0.5,0.5+float(ic)/41.)
            wt=amin1(1.0,wt)
            qo2(n)=(1.0-wt)*qo2(n)+wt*qold(n)
            qold(n)=qo2(n)
            endif        ! (store2(n).le.0.0)
            hold=store2(n)
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals
            store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n)-2.0*qwdr(n,month_now))*div
!                     qwdr multiplied by 2 to offset /2 in div     

!     rev. 9.9.52  Jan.  14/15  - NK: Fixed bug for channel store < 0 for withdrawals
          if(store2(n).le.0.0)then
!           cut off withdrawals                  
            store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
            if(store2(n).le.0.0)then
!             if store2 still = 0 then reduce calculated outflow
!             but this should not happen as it never did before 
!             withdrawals were added to the code
              store2(n)=0.0
              qo2(n)=store1(n)/div+qi1(n)+qi2(n)-qo1(n)         
            endif
          endif
         else
!          convergence to 3%
           GO TO 16
         endif
      end do    ! ic=1,20
       
  16  continue
       
!     rev. 9.5.46  Dec.  23/08  - NK: trying to fix problem with -ve storage. Changed conditional to .lt.
!     rev. 9.6.05  Apr.  06/10  - NK: added store_error_flag for -ve storage grids
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      if(store2(n).lt.0.0.or.qo2(n).lt.0.0)then
          store_error_flag(n)='true'
          qo2(n)=store1(n)/div+qi1(n)+qi2(n)-qo1(n)
          store2(n)=0.0
          dtmin=a6
!     rev. 9.5.47  Dec.  26/08  - NK: add flwinitflg to warn about initial flows
!         if we are early in the first event,possibly
!         the initial flows are too low.
!         This can happen if run is started in mid-winter
!         and the lowest downstream flow station is not at
!          outlet.
!         Might be fixed by increasing initial flow in the str file
  
          flwinitflg='y'
      endif     ! (store2(n).lt.0.0)
      
      if(qo2(n).gt.0.000001)then
!         CALCULATE THE VELOCITY THROUGH EACH SQUARE
!         CALCULATE TRAVEL TIME FOR MAXIMUM VELOCITY.
          at=store2(n)/qo2(n)
!         SELECT MIN TRAVEL TIME FOR THE TIME STEP CALC
          dtmin=amin1(at,dtmin)
          dtmin=amax1(dtmin,a6)   ! dtmin > a6 no matter what
!     rev. 9.9.49  Jan.  06/14  - NK: Added courantflg
          if(at.lt.a6)then
              courantflg=.false.
          endif
     
!         DTMIN IS THE TIME REQUIRED TO COMPLETELY DRAIN
!         THE FASTEST EMPTYING ELEMENT

!         CALCULATE THE CHANNEL STATE FOR GRAPHICAL OUTPUT:
          i=yyy(n)
          j=xxx(n)
          atemp=qo2(n)/(0.4*bnkfll(n))+1.0

!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
          bankfull(i,j)=qo2(n)/bnkfll(n)*100.0
          bankfull(i,j)=amax1(0.1,bankfull(i,j))
          bankfull(i,j)=amin1(10000.,bankfull(i,j))

!         TO PREVENT INTEGER UNDERFLOW OR OVEFLOW:  
          atemp=amax1(atemp,1.0)
          atemp=amin1(atemp,99.0)

          istate(i,j)=int(atemp)
!         istate(i,j)=min(istate(i,j),99)
!         istate(i,j)=max(istate(i,j),1)
      endif
              
      if(iopt.ge.4)write(55,1002)i,j,n,istate(i,j),bnkfll(n),qo2(n)
 1002 format(' i,j,n,istate,bnkfull,qo2/',4i5,2f10.3)

      return
      
      end subroutine rt_channel
              