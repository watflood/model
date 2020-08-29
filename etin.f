      SUBROUTINE etin(mon,jan,ju)
      
!***********************************************************************
!    Copyright (C) 1996 by Todd Neff  
        
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

! Organizes the calculation of PET and calculates the value of the 
! temperature constraint (fpet2) on PET 

!     REV. 8.25  - May.   22/97  -  changed calc for pet in etin 
!     REV. 8.99g - Feb.   7/2000 -  added ttoinit to init evaporation
!     REV. 9.00  - Mar.     2000 -  TS: CONVERTED TO FORTRAN 90 
!     rev. 9.2.40  Jun.  09/06  - NK: added tto(),ttomin(),ttomax() to resume
!     rev. 9.9.68  Apr.  29/15  - NK: Fixed tto reset with resume

!***********************************************************************

      USE area_watflood
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      integer  :: ju,jan,n,ii,mon
      real*4   :: ddh
      character*1 :: firstpass   ! intermediate variable

      data firstpass/'y'/

!        DATA ittoflg/0/    moved to spl9

!     LOOP INITIALIZES THE DEGREE DAYS TTO(N)
!     see also soilinit
c      if(firstpass.eq.'y')then
      if(jan.eq.1)then   
        do n=1,naa   
          if(resumflg.ne.'y')then
            tto(n)=tton            ! tton set in the par file
	      if(ju.ge.90.and.ju.le.270.and.tton.lt.1.0
     *                   .and.n.eq.1)then
	        print*,'Run is started during the growing season and'
	        print*,'tton (in the par file) should be given some value'
              print*,'See manual sect. 2.4.2'
              print*
	        pause 'Hit enter to continue (in etin.for @ 43)'
            endif 
            ttomin(n)=1.0e+32
            ttomax(n)=-1.0e+32
!           SET DEFAULT VALUE FOR FPET2
!           WARNING WILL BE USED IF NOT ALTERED LATER  
            fpet2(n)=1.0
          endif
        end do
      endif

      if(ju.eq.1)then
!       reset the values on Jan.1 each year
        do n=1,naa   
!     rev. 9.9.68  Apr.  29/15  - NK: Fixed tto reset with resume
c          if(resumflg.ne.'y')then   ! take it out!!!! it still needs to ber reset
            tto(n)=0.0
            ttomin(n)=1.0e+32
            ttomax(n)=-1.0e+32

!           SET DEFAULT VALUE FOR FPET2
!           WARNING WILL BE USED IF NOT ALTERED LATER  
            fpet2(n)=1.0
c          endif
        end do
      endif
      
!     rev. 9.7.28  Jun.  14/11  - NK: Add degree_day for lake_ice_factor  dd_ice
      if(ju.eq.274)then
        do n=1,naa
          dd_ice(n)=0.0
        end do
      end if
      if(ju.eq.90)then
        do n=1,naa
          dd_thaw(n)=0.0
        end do
      end if


!     SETS THE FLAG TO RESTART THE ACCUMULATION OF DEGREE-DAYS
!     FROM THE START OF MARCH (MON=3)

!     MAXIMUM VALUE OF FPET2(n)=1.0

!     FLGTMP1 SHUTS OFF THE ACCUMULATION OF DEGREE-DAYS ONCE
!     THE CONSTRAINT REACHES A VALUE OF 1.  FPET2(n)=1

!     ACCUMULATES DEGREE DAYS (TTO(n)).  CALCULATES THE CONSTRAINT
!     FPET2(n) ON THE PET IN EACH ELEMENT

      do n=1,naa
         ddh=tempv(n)/24.
         tto(n)=tto(n)+ddh     !! ACUMMULATION OF DEGREE-DAYS
         dd_ice(n)=dd_ice(n)-ddh
         dd_thaw(n)=dd_thaw(n)+ddh
!           assumes it gets here hourly

         ttomin(n)=amin1(ttomin(n),tto(n))
         ttomax(n)=amax1(ttomin(n),tto(n))

!        CHECK FOR DIV BY 0      

!      fix fix   ttomin thing needs to be checked out
!                problems with a summer start -evap ramps up 
!                with degree day from start of run
!                degree days have accumulated before the start.

!         if(abs(tempa3-ttomin(n)).gt.0.1e-20)then
!          fpet2(n)=(tto(n)-ttomin(n))/(tempa3-ttomin(n))
!         endif

         fpet2(n)=(tto(n)-ttomin(n))/tempa3
!           seems ok. checked out nk  Jul.6/06
         fpet2(n)=amax1(0.02,fpet2(n))
         fpet2(n)=amin1(1.00,fpet2(n))
!        ITTOFLG WILL BE 0 IF TTOMIN NEVER FALLS BELOW TEMPA3
!        TEMPA3 SHOULD BE GIVEN A HIGHER VALUE (IN THE PAR FILE)
!        THIS LINE WAS COMMENTED OUT, I THINK IT SHOULD BE THERE FRANK S DEC/2000
         if(tempa3.gt.ttomin(n)) ittoflg=1
      end do

!     IF FLGEVP2=1 USE TABLE PET
!     IF FLGEVP2=2 USE HARGREAVES PET
!     IF FLGEVP2=3 USE PREISTLEY-TAYLOR PET
!     NOTE: PET = 0 IF TEMP < -5 DEG. C

      if(flgevp2.eq.1.)then

!        print*,' in etin calculating pan evap'
!        USE PAN EVAPORATION AS THE POT EVAP
         do ii=1,classcount
            do n=1,naa
               if(tempv(n).ge.-5.0)then
!                 non-water class
                  pet(n,ii)=evap(ii,mon)
                  fpet2(n)=1.0
               else
                  pet(n,ii)=0.
               endif
            end do
         end do
      elseif (flgevp2.eq.2.) then
!	   print*,' in etin calling etharg'
         call etharg(mon,ju)
!          water class use table for cold temperatures
      elseif (flgevp2.eq.4.) then
!	   print*,' in etin calling etharg'
         call etharg_beta(mon,ju)
!          water class use table for cold temperatures
      elseif (flgevp2.eq.3.) then
         call etpriest(mon)
         do n=1,naa
            do ii=1,classcount
!              TN: APRIL, 1997 
!              REV. 8.25 - May.  22/97 - CHANGED CALC FOR PET IN ETIN 
               pet(n,ii)=alb(ii)/albe*pet(n,ii)
            end do
         end do
      endif

      firstpass='n'

      RETURN

      END SUBROUTINE etin

