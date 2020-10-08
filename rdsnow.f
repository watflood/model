      SUBROUTINE rdsnow(time,jan,new)

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
! SUBROUTINE TO READ DATA REQUIRED FOR SNOW ROUTINES

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90

!***********************************************************************

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

        CHARACTER(1) :: new
      integer  :: jan,n,i,j,it31a,it31b,nhrs,mz,idt33
      real*4   :: tlast,time,tst31a,tst31b

!-------------------------------------------------------------------
!     FIRST TIME STEP STUFF  
!--------------------------------------------------------------------
!     JAN IS <= 1 ON FIRST TIME THROUGH ONLY
!     JAN IS SET TO 2 AFTER THE FIRST EXECUTION OF SUBROUTINE RUNOF5
   

!--------------------------------------------------------------------  
!     JD'S OLD STUFF:

      if(jan.eq.1)then



	go to 999





!        OPEN THE DISTRIBUTED SNOW GAUGE AND TEMPERATURE DATA FILES
  
!        FLN(16)=yymmdd.tmn
!        FLN(17)=yymmdd.tmx
!        FLN(18)=yymmdd.dsn

!         open(unit=46,file=fln(16),status='unknown',err=999)
!         open(unit=47,file=fln(17),status='unknown',err=999)
!         open(unit=48,file=fln(18),status='unknown',err=999)

!        REWIND FILES --  NEEDED TO REINITIALIZE FOR OPTIMIZATION   
 
         tlast=-10
        
!        DAILY SNOWFALL TIMESTEP AND CONVERSION FACTOR
         tlst31=-9999
!         read(48,5006,err=999)idt31,conv31
!         if(numa.eq.0) write(57,5001)idt31,conv31

!        DAILY MINIMUM TEMPERATURE TIMESTEP AND CONVERSION FACTOR 
         tlst32=-9999
!         read(46,5006)idt32,conv32
!         if(numa.eq.0) write(57,5003)idt32,conv32

!        DAILY MAXIMUM TEMPERATURE TIMESTEP AND CONVERSION FACTOR 
         tlst33=-9999
!         read(47,5006)idt33,conv33
!         if(numa.eq.0) write(57,5008)idt33,conv33

!-------------------------------------------------------------------
!        INITIALIZE DAILY MIN TEMPERATURE DATA - ALREADY DISTRIBUTED

         tlst32=time
         if(numa.eq.0) write(57,5004)time+1

!         read(46,*) 
!         do i=it,ib,-1
!            read(46,1000)(tmn(i,j),j=2,xcount-1)
!            if(numa.eq.0) write(57,1000)(tmn(i,j),j=2,xcount-1)
!         end do

!        PUT INTO VECTOR FORMAT
         do n=1,naa
            i=yyy(n)
            j=xxx(n)
            tmin1(n)=conv32*tmn(i,j)
         end do

!------------------------------------------------------------------
!        DAILY MIN TEMPERATURE DATA FOR SECOND DAY
!          - NEEDED FOR DIST'D TEMPERATURE 2ND HALF OF DAY 1

         tlst32=time
         if(numa.eq.0) write(57,5004)time+1+float(idt32)

!         read(46,*) 
!         do i=it,ib,-1
!            read(46,1000)(tmn(i,j),j=2,xcount-1)
!            if(numa.eq.0) write(57,1000)(tmn(i,j),j=2,xcount-1)
!         end do

!        PUT INTO VECTOR FORMAT
         do n=1,naa
            i=yyy(n)
            j=xxx(n)
            tmin2(n)=conv32*tmn(i,j)
         end do

      endif

!     * * * * *  END OF INITIALIZATION  * * * * * 

!  -----------------------------------------------------------------
!     READ MET DATA 
!                   - idt31 = snowfall data time interv.(hrs)
!                   - idt32 = min. temperature time interv.(hrs)
!                   - idt33 = max. temperature time interv.(hrs)   
!  -----------------------------------------------------------------
!     DAILY SNOWFALL DATA - ALREADY DISTRIBUTED SPATIALLY
!     FOR EACH TIME STEP. STILL MUST SPLIT BETWEEN CLASSES
 
!     READ SNOWFALL DATA IF NEXT TIME INCREMENT IS REACHED

      if(idt31.le.0)then
         print*, 'idt31 in rdsnow .le. 0.0'
         print*, 'idt31 assumes as 1'
         print*
         idt31=1
      endif
      tst31a=time/idt31
      tst31b=tlst31/idt31
      it31a=int(tst31a)
      it31b=int(tst31b)
 
      if(int(time/idt31).ne.int(tlst31/idt31))then
         tlst31=time
         if(numa.eq.0) write(57,5002)time+1

!        CHECK THAT TIME HASN'T RUN OUT
         nhrs=nr
         mz=int(time)+1
         if(mz.gt.nhrs)then
            do n=1,naa
               dsnow(n)=0.0
            end do
         else

!------------------------------------------------------------------
!           CHECK IF ANY SNOW FALL DURING TIME INTERVAL
!           IF YES THEN READ AND ASSIGN TO DSN(i,j)
!           IF NO THEN DSN(i,j) = 0 FOR ALL i,j's

!            read(48,5005) k

!            if(k.gt.0) then
!               do i=it,ib,-1
!                  read(48,1000)(dsn(i,j),j=2,xcount-1)
!                  if(numa.eq.0) write(57,1000)(dsn(i,j),j=2,xcount-1)
!               end do
!            else
               do i=it,ib,-1
                  do j=2,xcount-1
                     dsn(i,j)=0.0
                  end do
                  if(numa.eq.0) write(57,1000)(dsn(i,j),j=2,xcount-1)
               end do
!            endif

!           PUT INTO VECTOR FORMAT
            do n=1,naa
               i=yyy(n)
               j=xxx(n)
               dsnow(n)=conv31*dsn(i,j)
            end do

         endif
 
      endif
 
!-------------------------------------------------------------------
!     DAILY MIN TEMPERATURE DATA - ALREADY DISTRIBUTED
!     NB: DAILY MIN TEMPERATURE IS READ 1 TIME STEP AHEAD OF 
!         THE SIMULATION IN ORDER THAT TEMPERATURE CALCULATIONS
!         CAN BE CARRIED OUT FOR THE SECOND HALF OF THE DAY - SEE
!         CALC'S IN TEMPER.FOR 

!     READ TEMPERATURE DATA IF NEXT TIME INCREMENT IS REACHED

      if(int(time/idt32).ne.int(tlst32/idt32))then
         tlst32=time
         if(numa.eq.0) write(57,5010)time+1+float(idt32)
!        IS TIME+1+IDT32 BECAUSE WE ARE READING MIN TEMPS 1 STEP AHEAD

!        ASSIGN TMIN2 TO TMIN1 
         do n=1,naa
            tmin1(n)=tmin2(n)
         end do

!        CHECK THAT TIME HASN'T RUN OUT
!        NB: TMIN IS BEING READ 1 TIME STEP AHEAD OF THE OTHER DATA
!            - LEAVE IN OLD VALUES FOR TMIN 
         nhrs=nr
         mz=int(time+1)+idt32

         if(mz.gt.nhrs)then
            do n=1,naa
               tmin2(n)=tmin2(n)
            end do
         else
!            read(46,*) 
!            do i=it,ib,-1
!               read(46,1000)(tmn(i,j),j=2,xcount-1)
!               if(numa.eq.0) write(57,1000)(tmn(i,j),j=2,xcount-1)
!            end do

!           PUT INTO VECTOR FORMAT
            do n=1,naa
               i=yyy(n)
               j=xxx(n)
               tmin2(n)=conv32*tmn(i,j)
            end do
         endif
      
      endif

!-------------------------------------------------------------------
!     DAILY MAX TEMPERATURE DATA - ALREADY DISTRIBUTED

!     READ TEMPERATURE DATA IF NEXT TIME INCREMENT IS REACHED

      if(int(time/idt33).ne.int(tlst33/idt33))then
         tlst33=time
         if(numa.eq.0) write(57,5009)time+1

!        CHECK THAT TIME HASN'T RUN OUT
         nhrs=nr
         mz=int(time+1)
         if(mz.gt.nhrs) then
            do n=1,naa
!              USE LAST DAYS MAX TEMPERATURE	  
               tmax(n)=tmax(n)
            end do
        else
!            read(47,*) 
!            do i=it,ib,-1
!               read(47,1000)(tmx(i,j),j=2,xcount-1)
!               if(numa.eq.0) write(57,1000)(tmx(i,j),j=2,xcount-1)
!            end do

!           PUT INTO VECTOR FORMAT
            do n=1,naa
               i=yyy(n)
               j=xxx(n)
               tmax(n)=conv33*tmx(i,j)
            end do
        endif
        
      
      endif
	
  999 new='t'
      idt31=1
      idt32=1
      idt33=1
!      write(*,5999)new

! FORMATS:

 1000 format(16f10.1)
 5001 format(1x,'snow data time step=',i4,'  conv factor=',f6.2)
 5002 format(1x,'snowpack input data: time=',f5.0)
 5003 format(1x,'tmin data time step=',i4,'  conv factor=',f6.2)
 5004 format(1x,'temperature input data: time=',f5.0)
 5005 format(6x,i5)
 5006 format(i5,f6.1)
 5008 format(1x,'tmax data time step=',i4,'  conv factor=',f6.2)
 5009 format(1x,'max. temperature input data: time=',f5.0)
 5010 format(1x,'min.temperature input data: time=',f5.0)
 5999 format(' new = ',a1)


      RETURN

      END SUBROUTINE rdsnow


