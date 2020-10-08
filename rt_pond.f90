      subroutine rt_pond(div,thr,dtmin,jz,iz,time,date,tdum,n)      

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
     
!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing      
!     Written Nov. 2/17 NK   

      use area_watflood
      implicit none

      SAVE          ! SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT

      integer       :: l,ic,n,jz,iz,i,j
      real          :: hold,wt,div,thr,dtmin,time,tdum,at,amin1
      real          :: strlossvol,old
      character*14  :: date

!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
               if(strlossvol*100.0.gt.store1(n)+(qi1(n)+qi2(n))*div)then
                    write(98,*)'Warning:',time,'strloss ',strloss(n),' set = 0.0 for grid ',n
                    strloss(n)=0.0
                    warning(2)=.true.
                endif
!####################################################################################      

!     rev. 9.9.59  Feb.  18/15  - NK: In route: strloss option fracflg y/n
!     rev. 10.1.62 Jan.  08/17  - NK: Checkup on strloss effect on low flows
!       TS preferred way:             
        qi2(n)=qi2(n)+qr(n)+qstream(n)-strloss(n)
        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
      
        old=qo1(n)
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
        hold=0.1e+26
        ic=0
!c         using a power function        
          if(store2(n).gt.0.0)then
            do while(abs(hold-store2(n)).gt.0.003*hold.and.ic.lt.20)
                ic=ic+1
                if(store2(n).gt.0.0)then
!                 have to do at least 3 iterations
                  qo2(n)=rlake(n)*store2(n)**1.75
!c                  qo2(n)=qo2(n)*lake_ice_factor(l,month_now)
                  wt=amax1(0.1,float(ic)/21.0)
                  qo2(n)=(1.0-wt)*qo2(n)+wt*old
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
!c                  qo2(n)=(1.0-wt)*qo2(n)+wt*qold(n)
                  old=qo2(n)
!c                  qold(n)=qo2(n)
                  hold=store2(n)
                  store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
                endif
            end do
          else     !(b3(l).gt.0.0)
            qo2(n)=0.0
            store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
!           store2(n) is allowed to go -ve     
          endif

!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
          i=yyy(n)
          j=xxx(n)
          bankfull(i,j)=-1.0
          if(qo2(n).gt.0.000001)then
!             CALCULATE THE VELOCITY THROUGH EACH SQUARE
!             CALCULATE TRAVEL TIME FOR MAXIMUM VELOCITY.
              at=store2(n)/qo2(n)
!             SELECT MIN TRAVEL TIME FOR THE TIME STEP CALC
              dtmin=amin1(at,dtmin)
              dtmin=amax1(dtmin,a6)   ! dtmin > a6 no matter what
              if(at.lt.a6)then
                courantflg=.false.
              endif
!             DTMIN IS THE TIME REQUIRED TO COMPLETELY DRAIN
!             THE FASTEST EMPTYING ELEMENT

          endif
          
      return
        
      end subroutine rt_pond