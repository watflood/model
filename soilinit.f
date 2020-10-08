      SUBROUTINE soilinit

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
!     REV. 8.31  - June   3/97 -  added initial uzs values in evap.par
!     REV. 8.96  - Apr.  26/99 -  lower zone function related to nbsn
!     REV. 8.99f - Jan. 7/2000 -  changed uzs calcs re: shari's data
!
!     REV. 9.00  - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90
!     REV. 9.03    Nov.  2000  -  TS: ADDED WATFLOOD SWAMP ROUTING 
!     rev. 9.2.42  Jun.  20/06  - NK: water class included in the water balance
!
!***********************************************************************

      use area_watflood
	implicit none

      integer       :: nflag,jj,ios,i,j,n,ii
      real (4)      :: ff1                    !,xxx1
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

d      if(iopt.eq.2)print*,' checkpoint 1 in soilinit'

	do i=1,ycount
         do j=1,xcount
            sn1(i,j)=0.0
         end do
      end do


!     Initialize these - unknown values so set = 0.0
      do n=1,naa
	  eloss(n)=0.0
        do ii=1,classcount
          wcl(n,ii)=0.0
	    v(n,ii)=0.0
        end do
      end do

      write(98,*)'Info: Free water in snowpack set = 0.0'
	write(98,*)'Info: Interceptions storage  set = 0.0'

      sqlz=0.0
      slzinflw=0.0

d      if(iopt.eq.2)print*,' checkpoint 2 in soilinit'


      do n=1,naa
            i=yyy(n)
            j=xxx(n)
            d2(n)=0.000001
            nflag=1
!           A VALUE OF NFLAG > 0 INDICATES THERE IS SNOW ON THE GROUND
    

!           REV. 7.5 SEPERATE SNOW COVERED AND BARE GROUND
!           SET INITIAL VALUES FOR EACH SOIL & EACH ELEMENT:
            do ii=1,classcount
              if(nclass(ii).ne.'water     ')then ! added Jul.11/06 nk
!               aclass(n,ii)=amax1(aclass(n,ii),0.0)
!               aclass(n,ii)=amin1(aclass(n,ii),1.0)
                v(n,ii)=0.0
!               the next two were added Jun. 3/02  nk
	          ev(n,ii)=0.0
	          evt(n,ii)=0.0
	          intevt(n,ii)=0.0
                d1(n,ii)=0.000001
                d1fs(n,ii)=0.000001
!     rev. 8.99f - Jan. 7/2000-     changed uzs calcs re: shari's data
                uzs(n,ii)=retn(ii)*api(n,ii)/spore(ii)
                uzsfs(n,ii)=retn(ii)*api(n,ii)/spore(ii)
	          uzs(n,ii)=amin1(retn(ii),uzs(n,ii))
	          uzsfs(n,ii)=amin1(retn(ii),uzsfs(n,ii))
                sumf(n,ii)=0.0
                sumffs(n,ii)=0.0
              else
!     rev. 9.2.42  Jun.  20/06  - NK: water class included in the water balance
!               this line added Sep. 13/06 nk
	          evt(n,ii)=0.0
	          intevt(n,ii)=0.0
                uzs(n,ii)=0.0       !used for the water class water balance only
                uzsfs(n,ii)=0.0     !used for the water class water balance only
                d1(n,ii)=0.0        !just so rff files don't have junk
                d1fs(n,ii)=0.0      !just so rff files don't have junk
	        endif

!             SET INITIAL UNSATURATED ZONE EFFECTIVE POROSITY:
!             THE EFFECTIVE POROSITY DEPENDS ON SOIL MOISTURE
!             WHICH IS ENTERED FOR EACH GRID
!             effpor is for the intermediate zone so same for bare & sca
              effpor(n,ii)=spore(ii)-api(n,ii)
              effpor(n,ii)=amax1(0.0001,effpor(n,ii))
              effpor(n,ii)=amin1(spore(ii),effpor(n,ii))
              if(api(n,ii).lt.0.0)then
!                write(*,6027)n,ssmc(i,j),effpor(n)
                api(n,ii)=0.0
              endif
            end do
      end do

d      if(iopt.eq.2)print*,' checkpoint 3 in soilinit'

d      if(iopt.eq.2)print*,' checkpoint 4 in soilinit'


!     These valuse may be reset in etin.for
!     but are initialized here to have values in the soil_init.r2c file
!     if this file is written only at start of spl run nk Nov. 13/06
      do n=1,naa   
        tto(n)=tton            ! tton set in the par file
        ttomin(n)=1.0e+32
        ttomax(n)=-1.0e+32
        fpet2(n)=1.0
      end do

d      if(iopt.eq.2)print*,' checkpoint 4a in soilinit'


!     * * * * END OF THE INITIALIZATION (JAN) * * * *

! FORMATS

 6008 format(/' Error: DA <= 0.0 in el. n=',i5/)
 6100 format
     *('0','    n    i    j        da       lzs     qbase in runof5'/)
 6102 format(' ',3i5,3e10.3)

      if(iopt.ge.1)then
        print*
        print*,'Soil initialization completed'
        print*
      endif  


      RETURN

      END SUBROUTINE soilinit

