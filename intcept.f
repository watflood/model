        SUBROUTINE intcept(mon)

!***********************************************************************
!    Copyright (C) 1996 Todd Neff and Nicholas Kouwen 
        
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
! PROGRAM BY: Todd A.M. Neff,  MARCH 1996

! Subprogram intcept.f for the Watflood program runof5.f
! calculates ineterception v(n,ii);interception evaporation intev(n,ii);
! throughfall r(n,ii).  Based on Linsley's equation

!     REV. 8.6   - Oct.  28/97  -  replaced todd's incept routine 
!     REV. 9.00  - Mar.   2000  -  TS: CONVERTED TO FORTRAN 90
!     rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
!     rev. 9.9.26  Sep.  18/14  - NK: Added zero class bypass in intcept.f


!  vo           - interception storage during the last time step
!  ssumr        - sum of rain so far - used for plotting only
!  pint(n,ii)   -  the sum of intercepted rain
!  v(n,ii)      -  water in interception storage
!  intev(n,ii)  -  interception evaporation this time step
!  intevt(n,ii) -  the sum of intev
!  flint(ii)    -  an interception flag om the param file

!***********************************************************************

      USE area_watflood

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      implicit none
      SAVE

c      REAL*4 :: pint(17500,99),deficit(99)
!      moved to area_watflood & allocated

      integer  :: n,i,j,ii,mon


d     if(iopt.eq.2)print*,'in intcept @ 35'

      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        do ii=1,classcount
!     rev. 9.9.26  Sep.  18/14  - NK: Added zero class bypass in intcept.f
        if(aclass(n,ii).gt.0.0)then   ! skip if class area = 0.0

!       NOTE: classcount IS THE WATER CLASS AND IS DONE BELOW

!!!!!!! FIX THIS BOTH IF'S DO THE SAME THING !!!!!!!

         if(flint(ii).gt.0.9999)then
!           FLINT IS SET IN THE PARAM FILE
            if(h(mon,ii).gt.0.0)then
!              INT EVAP IS INCREASED BY FPET FOR ALL VEGETATION   
!              BECAUSE OF THE INCREASED ROUGHNESS AND ADVECTION
!              FIRST CHECK THAT PRECIP OCCURS
               if(p(i,j).gt.0.0)then
!                 DURING PRECIP EVENT THERE IS NO EVAPORATION FROM 
!                 THE LEAVES; WATER IS ADDED TO INTERCEPTION STORAGE 
!                 TO THE MAX INTERCEPTION STORAGE AND THE REMAINDER 
!                 REACHES THE GROUND...WE CAN DO AN EITHER-OR

!				This is different from what Todd has inhis thesis 
!                 which states that during precip the int. evaporation 
!                 occurs at the pet rate. But nk figures that pet=0.0
!                 during rain.

                  deficit(ii)=h(mon,ii)-v(n,ii)
                  if(p(i,j).ge.deficit(ii))then
!                    PRECIP IS GREATER THAN REMAINING INTERCEPTION 
!                    DEFICIT AND INTERCEPTION IS MAXED OUT
!                    THRUFALL IS PRECIP - INTERCEPTED WATER 
                     pint(n,ii)=pint(n,ii)+deficit(ii)
!                    CALCULATE THRUFALL   
                     r(n,ii)=p(i,j)-deficit(ii)
                     v(n,ii)=v(n,ii)+deficit(ii)
                     intev(n,ii)=0.0
                  else
!                    PRECIP IS LESS THAN INTERCEPTION DEFICIT 
!                    AND ALL PRECIP IS ADDED TO INTERCEPTION STORAGE
                     pint(n,ii)=pint(n,ii)+p(i,j)
                     r(n,ii)=0.0
                     v(n,ii)=v(n,ii)+p(i,j)
                     intev(n,ii)=0.0
                  endif

               else

!                 ELSE, NO RAIN OCCURS - NO WATER REACHES THE GROUND
!                 THIS IS THE SAME FOR BARE & SNOW COVERED AREAS
                  r(n,ii)=0.0

!                 IF THERE IS WATER IN INTERCEPTION STORAGE THEN 
!                 EVAPORATE IT AT FPET TIMES THE POTENTIAL RATE UNTIL
!                 THERE IS NONE.  RECALCULATE LINSLEY CONSTANT BASED ON
!                 CURRENT AMOUNT OF STORAGE.

                  if(v(n,ii).gt.0.0)then
!                    THERE IS WATER TO EVAP AT FPET TIMES THE POT RATE
                     if(fpet(ii)*pet(n,ii).gt.v(n,ii))then
!                       ALL WATER IN INT STORAGE WILL BE EVAPORATED 
                        intev(n,ii)=v(n,ii)
                        v(n,ii)=0.0
!                       PET IS REDUCED BY THIS AMOUNT FOR NEXT CALC
!                       FOR CALCULATING EVAP FROM SOIL IN AET.FOR
                     else
!                       NOT ALL WATER IN INT STORAGE WILL BE EVAPORATED
                        intev(n,ii)=fpet(ii)*pet(n,ii)
!                       AND INTERCEPTION STORAGE IS REDUCED BY
                        v(n,ii)=v(n,ii)-intev(n,ii)
                     endif
                  else
!                    THERE IS NO WATER IN INTERCEPTION STORAGE
                     intev(n,ii)=0.
                  endif
               endif

               intevt(n,ii)=intevt(n,ii)+intev(n,ii)

               sum_et(n,ii)=sum_et(n,ii)+intev(n,ii)

!              LOSS IS FROM BARE & SNOW COVERED AREA ALIKE
               eloss(n)=eloss(n)+intev(n,ii)*aclass(n,ii)
               if(pint(n,ii).lt.0.)     pint(n,ii)=0.
               if(v(n,ii).lt.0.)        v(n,ii)=0.
               if(v(n,ii).gt.h(mon,ii)) v(n,ii)=h(mon,ii)
               if(intev(n,ii).lt.0.)    intev(n,ii)=0.
               vo(n,ii)=v(n,ii)
               ssumr(n,ii)=ssumr(n,ii)+r(n,ii)
           endif

          else
!           IF FLINT EQ 0:
!           NO INTERCEPTION IS CALCULATED (eg. WATER, GLACIERS)
            intev(n,ii)=0.
            intevt(n,ii)=intevt(n,ii)+intev(n,ii)
            r(n,ii)=p(i,j)

!           REV. I ADDED THIS HERE - NK: JAN 29/97
            ssumr(n,ii)=ssumr(n,ii)+r(n,ii)
          endif

c          if(iopt.ge.1)then
c	      if(n.eq.nnprint.and.ii.eq.iiprint)then
c	       write(64,6401)n,ii,p(i,j),v(n,ii),deficit(ii),pint(n,ii),
c     *         r(n,ii),intev(n,ii),intevt(n,ii)
c6401         format(2i5,7e12.4)
c	      endif
c          endif

        endif  !  aclass()>0
        end do

!     rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
        r(n,classcount)=p(i,j)   !impervious area
        intev(n,classcount)=0.0
        v(n,classcount)=0.0

      end do
d     if(iopt.eq.2)print*,'in intcept @ 155'

      RETURN

      END SUBROUTINE intcept


