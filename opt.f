      SUBROUTINE opt(*,*)

!***********************************************************************
!    Copyright (C) 1975  Monroe (NOAA) 
        
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
!   MAY 25,1975    NOAA TECH MEMO NWS HYDRO-12

!***********************************************************************

      use area_watflood
	implicit none

!     REPLACED AREAIZY (TO SAVE VARIABLES) WITH SAVE()
      SAVE 

      integer  :: lc,ittn,izy,ncoun,icoun,ldelt,nsave,ifirs,ll,iflg,i,ii
      real*4   :: ys,yx1,yy1,cc1,cc2,cc0,cd

!     THESE THINGS NO LONGER ARE SET TO 0 EACH TIME OPT IS CALLED BEING
!     LEFT TO THEIR OWN DEVICES IN UNIX -- NEED TO BE ZEROED
      DATA lc,ittn,izy,ncoun,icoun,ldelt,nsave,ys,yx1,yy1,ifirs,ll/12*0/

      if(nstart.gt.0) GO TO 2

!     INTIALIZATION ROUTINE 
!     check for paramaters to within constraints
      iflg=0
      do 1 i=1,numa
         les(i)=0
         ba(i)=a(i)
         b(i)=a(i)
         iclosl(i)=0
         iclosh(i)=0
!        nper:  if = 0 delta = absolute value
!               if = 1 delta = a multiplier
         if(nper.le.0)then
            odelta(i)=ddelta(i)
         else
            odelta(i)=abs(ddelta(i)*a(i))
         endif

         if(checkh(i)-checkl(i).le.3.0*odelta(i))then
!           THE CONSTRAINTS ARE TO CLOSE TOGETHER
            write(6,6510)i
            write(52,6510)i
            iflg=1
         endif

         cc1=a(i)-1.01*odelta(i)
         if(cc1.le.checkl(i))then
            write(6,6500)i,a(i),cc1,ddelta(i),checkl(i),checkh(i)
            write(52,6500)i,a(i),cc1,ddelta(i),checkl(i),checkh(i)
            iflg=-1
         endif

         cc2=a(i)+1.01*odelta(i)
         if(cc2.ge.checkh(i))then
            write(6,6500)i,a(i),ddelta(i),checkl(i),checkh(i)
            write(52,6500)i,a(i),ddelta(i),checkl(i),checkh(i)
            iflg=-1
         endif
    1 CONTINUE

      if(iflg.ne.0) STOP 'please change init values or constraints'
       
      lc=0
      ittn=1
      izy=0
      nnn=0
      ncoun=1
      icoun=0
      ifirs=0
      ldelt=0
      nstart=1
      nsave=0
!      write(52,6221)
      write(52,6003)(i,i=1,numa)
      write(52,6005)ncoun,nnn,ys,(a(i),i=1,numa)
    2 ys=optim
      nnn=nnn+1
      if(nnn.gt.maxn) RETURN 1
      if(ifirs.eq.1) GO TO 4
      yx1=optim
      yy1=yx1
      ifirs=1
    4 write(52,6005)ncoun,nnn,ys,(a(i),i=1,numa)
   44 if(les(ittn).eq.1) GO TO 14
      if(izy.gt.0) GO TO 8
      if(ys.le.yy1)then
         nsave=1
         yx1=ys
         yy1=ys
         write(52,6003)(i,i=1,numa)
      endif
    6 izy=izy+1
      ittn=izy
     
	if(les(izy).eq.1) GO TO 107
  108 ll=0

!     LOCAL EXCURSION ROUTINE

!     LOCAL EXCURSION WITH + DDELTA(I) FIRST

      a(izy)=a(izy)+odelta(izy)
      nsign(izy)=0
      if(iclosh(izy).eq.0) GO TO 7
      ll=ll+1
      GO TO 88
    7 ll=ll+1

      RETURN 2

    8 if(yx1.gt.ys) GO TO 11
   88 GO TO (9,10,12),ll
    9 a(izy)=a(izy)-2.0*odelta(izy)
      nsign(izy)=1
      if(iclosl(izy).eq.1) GO TO 10
      GO TO 7
   10 a(izy)=a(izy)+odelta(izy)
      nsign(izy)=0
      GO TO 12
   11 yx1=ys
   12 if(izy.lt.numa) GO TO 6
      ittn=1
      izy=0
      if(yy1.eq.yx1) GO TO 25
      yy1=yx1

      GO TO 210

!     LOCAL EXCURSION WITH - ODELTA(I) FIRST

   14 if(izy.gt.0) GO TO 16
      if(ys.le.yy1)then
         nsave=1
         yx1=ys
         yy1=ys
      endif
      write(52,6003)(i,i=1,numa)
  106 izy=izy+1
      ittn=izy
      if(les(izy).eq.0) GO TO 108
  107 ll=0
! ONCE AGAIN WHY IS THIS LINE USING DDELTA?
! CHANGED TO ODELTA BY FRANK S DEC/2000
c      a(izy)=a(izy)-ddelta(izy)
      a(izy)=a(izy)-odelta(izy)
      nsign(izy)=1
      if(iclosl(izy).eq.0) GO TO 15
      ll=ll+1
      GO TO 166
   15 ll=ll+1

      RETURN 2

   16 if(yx1.gt.ys) GO TO 19
  166 GO TO (17,18,20),ll
   17 a(izy)=a(izy)+2.0*odelta(izy)
      nsign(izy)=0
      if(iclosh(izy).eq.1) GO TO 18
      GO TO 15
   18 a(izy)=a(izy)-odelta(izy)
      nsign(izy)=1
      GO TO 20
   19 yx1=ys
   20 if(izy.lt.numa) GO TO 106
      ittn=1
      izy=0

      if(yy1.eq.yx1) GO TO 25
      yy1=yx1
!         WATCH THE JUMP INTO HERE!         
  210     if(nper.ne.0)then
             do 21 i=1,numa
                odelta(i)=abs(ddelta(i)*a(i))
   21        CONTINUE
          endif
          
          lc=0
          nsave=0
          write(52,6005) ncoun,nnn,yy1,(a(i),i=1,numa)
          write(52,6220)
          if(numa.gt.1)write(6,6220)
          ncoun=ncoun+1
         
!              PATTERN MOVE ROUTINE
         
          do 24 i=1,numa
             les(i)=nsign(i)
             ba(i)=a(i)
             a(i)=2.0*a(i)-b(i)
!              CHECK UPPER AND LOWER CONSTRAINTS
           
             cc0=a(i)-1.01*odelta(i)
             cd=a(i)+1.01*odelta(i)

             if(cc0.le.checkl(i))then
                iclosl(i)=1
                a(i)=ba(i)
             else
                iclosl(i)=0
              endif

             if(cd.ge.checkh(i))then
                iclosh(i)=1
                a(i)=ba(i)
             else 
                iclosh(i)=0
             endif

             b(i)=ba(i)
   24     CONTINUE
          
          RETURN 2

   25 lc=lc+1

!          DESTROY PRESENT PATTERN

      if(lc-1)225,26,28

  225 RETURN 1

   26 if(nsave.eq.1) GO TO 260

      do 27 i=1,numa
         a(i)=ba(i)
   27 CONTINUE
      icoun=icoun+1
      GO TO 30

   28 if(ldelt.ge.kc) RETURN 1

!          HALVE ODELTA(I) (RESOLUTION)

  260 nsave=0
      do 29 i=1,numa
         ddelta(i)=ddelta(i)*0.5
         odelta(i)=odelta(i)*0.5
! frank was wrong. Checked with NOAA tech memo 12 NK Feb. 05/10
!     THIS IS WRONG, I CHANGED IT TO WHAT IT SHOULD BE
!     FRANK S DEC 2000
c         odelta(i)=ddelta(i)*0.5
   29 CONTINUE
      ldelt=ldelt+1
   30 write(52,6031)icoun,ldelt
      if(numa.gt.1) write(6,6031)icoun,ldelt
      write(52,6005) ncoun,nnn,yy1,(a(i),i=1,numa)

      GO TO 44
  
! FORMATS

 6003 format(' trial run     criteria  ',50('      a(',i2,')'))
 6005 format(i6,i5,e16.7,50e11.4)
 6220 format(21x,'pattern move')
 6221 format(21x,'initial values of the coefficients')
 6031 format(20x,'pattern=',i4,' resolution=',i5)
 6500 format(' init a(',i2,') to close to its constraint',5e10.3)
 6510 format(' the constraints on a(',i3,') are too close together')          

      RETURN

      END SUBROUTINE opt

