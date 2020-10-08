!***********************************************************************
      subroutine wqroute(t)
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
! This subroutine routes the water quality variables: sediments and
! nutrients.  It will do so based on a mixed cell model (continuity)
! Modified: Nov/98 - LFL

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      integer :: n,lll,l,nnnext
      real*4  ::    wi,wold,t,woldn,woldp,wi2keep,weight

! initialize and reset values for each time step
      do 10 n=1,naa
! sediment        
        wi1(n)=wi2(n)
        wi2(n)=0.0
        wo1(n)=wo2(n)
!        wo2(n)=0.0  ! (try without this line nk)
        ss1(n)=ss2(n)
        ss2(n)=0.0
! nutrients (n-nitrogen ; p-phosphorus)
        wi1n(n)=wi2n(n)
        wi2n(n)=0.0
        wo1n(n)=wo2n(n)
        wo2n(n)=0.0
        ss1n(n)=ss2n(n)
        ss2n(n)=0.0
        wi1p(n)=wi2p(n)
        wi2p(n)=0.0
        wo1p(n)=wo2p(n)
        wo2p(n)=0.0
        ss1p(n)=ss2p(n)
        ss2p(n)=0.0
   10 continue    

	
      if(index.eq.0)then
        do 20 n=1,naa

   20   res(n)=0
        do 30 n=1,naa
          lll=next(n)
          do 30 l=1,noresv
            if(yyy(n).eq.ires(l).and.xxx(n).eq.jres(l)) res(n)=l
   30  continue
      endif


!*****   route the sediment   ****************

      do 40 n=1,naa
        lll=next(n)
!       when the slope <= 0.0 the element is not in the basin: skip route

        if(slope(n).gt.0.0)then
          if (store2(n).ne.0.0)then
!            wold=1.0e+25

            wold=wo2(n)

!           yrot is the ss conc. in kg/m^3  9nk0
!           wi is the sediment mass inflow in kg this dt  (nk)      


!            Mar. 23. 2004- Ek changed these lines

!            wi2(n)=wi2(n)+yrot(n)*qr(n)*t

            wi2(n)=wi1(n)+yrot(n)*(q1(n,3)+qint(n,3))*t

!      if(n.eq.130)then
!      print*,n,yrot(n),wi2(n),q1(n,3),qint(n,3)
!      endif  
  
	      do 50 ijk=1,50
              weight=amax1(0.5,real(ijk)/50.0)
!             ss2 is the suspended solids in kg   (nk 21/03/04   )
              ss2(n)=ss1(n)+(wi1(n)+wi2(n)-wo1(n)-wo2(n))/2

!             yfinal is the sediment concentration in kg/m^3  (nk)
              yfinal(n)=ss2(n)/store2(n)

!             wo2=sediment mass leaving during dt
              wo2(n)=weight*qo2(n)*yfinal(n)*t+(1.0-weight)*wo2(n)
              if(abs(wo2(n)-wold).lt.0.03*wold)goto 60
              wold=wo2(n)
  50        continue
  
!      if(n.eq.218.and.wi2(n).gt.0.0)then
!      print*,n,wi2(n),wi1(n),wi2(lll)
!      print*,wold,q1(n,3),qint(n,3)
!      endif  
  
    
  60        continue
          else
            ss2(n)=0.0
            wo2(n)=0.0
          end if
! decreased by the amount of deposition:            


!      REMOVED WITH IMPLEMENTATION OF WETLAND DECAY ROUTINE!!!!!!!

!          wi2(lll)=(1.0-(sdep/100.))*(wi1(lll)+wo2(n))
          wi2(lll)=(wi1(lll)+wo2(n))

      if(wi1(lll).lt.0.0.or.wi2(lll).lt.0.0.or.
     *                       wo2(n).lt.0.0.or.wo2(n).lt.0.0)
     *      write(98,*)n,next(n),wi1(lll),wi2(lll),wo1(n),wo2(n)
      if(next(n).eq.nnprint)
     *         write(98,*)n,next(n),wi1(lll),wi2(lll),wo1(n),wo2(n)
        endif
   40 continue


!*****   route the nutrients   ****************

      do 70 n=1,naa
        lll=next(n)
!       when the slope <= 0.0 the element is not in the basin: skip route
        if(slope(n).gt.0.0)then
          if (store2(n) .ne. 0.0) then

!            woldn=1.0e+25




!            Mar. 23. 2004- Ek changed these lines
            woldn=wo2n(n)
!            wi2n(n)=wi2n(n)+cronrot(n)*qr(n)*t
                     
            
!*****************************            
            wi2n(n)=wi1n(n)+cronrot(n)*(q1(n,3)+qint(n,3))*t
!*************************************


            do 80 ijk=1,50
              weight=amax1(0.5,real(ijk)/50.0)
              ss2n(n)=ss1n(n)+(wi1n(n)+wi2n(n)-wo1n(n)-wo2n(n))/2
              nfinal(n)=ss2n(n)/store2(n)
              wo2n(n)=weight*qo2(n)*nfinal(n)*t+(1.0-weight)*wo2n(n)
              if(abs(wo2n(n)-woldn).lt.0.03*woldn)goto 100
              woldn=wo2n(n)
  80        continue
 100  continue
          else
            ss2n(n)=0.0
            wo2n(n)=0.0
          end if

!      REMOVED WITH IMPLEMENTATION OF WETLAND DECAY ROUTINE!!!!!!!
!         decreased by the amount of decay:            
!          wi2n(lll)=(1-(ndec/100.))*(wi2n(lll)+wo2n(n))

          wi2n(lll)=(wi2n(lll)+wo2n(n))

        endif
   70 continue

      do 170 n=1,naa
        lll=next(n)
!       when the slope <= 0.0 the element is not in the basin: skip route
        if(slope(n).gt.0.0)then
          if (store2(n) .ne. 0.0) then

!            Mar. 23. 2004- Ek changed these lines

            woldp=wo2p(n)
!            wi2p(n)=wi2p(n)+croprot(n)*qr(n)*t

            wi2p(n)=wi1p(n)+croprot(n)*(q1(n,3)+qint(n,3))*t

            do 180 ijk=1,50
              weight=amax1(0.5,real(ijk)/50.0)
              ss2p(n)=ss1p(n)+(wi1p(n)+wi2p(n)-wo1p(n)-wo2p(n))/2.
              pfinal(n)=ss2p(n)/store2(n)
              wo2p(n)=weight*qo2(n)*pfinal(n)*t+(1.0-weight)*woldp
              if(abs(wo2p(n)-woldp).lt.0.03*woldp)goto 200
              woldp=wo2p(n)
 180        continue
 200  continue
          else
            ss2p(n)=0.0
            wo2p(n)=0.0
          end if

!      REMOVED WITH IMPLEMENTATION OF WETLAND DECAY ROUTINE!!!!!!!
!         decreased by the amount of decay:            
!          wi2p(lll)=(1-(pdec/100.))*(wi2p(lll)+wo2p(n))
          wi2p(lll)=(wi2p(lll)+wo2p(n))

        endif
  170 continue
 
      return
      end
