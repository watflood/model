      SUBROUTINE etharg(mon,ju)
      
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

! THIS SUBROUTINE CALCULATES HARGREAVES PET FOR EACH ELEMENT
! HARGREAVES EQUATION FROM A HANDBOOK OF HYDROLOGY
! - NOTE THAT THERE ARE SOME CHANGES - SEE ORIGINAL PAPER
!   REFERENCED IN MASc THESIS T.Neff 1996

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90

! akt - reduction due to excessive humidity (>54%)

!***********************************************************************

      USE area_watflood
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      real*4 :: akt1,petn,dlta,dr,sindlta,cosdlta,tandlta,ws0,so0
      integer  :: mon,ju,n,ii

      if(hu(mon).lt.54)then
         akt = .125
      else
         akt = 0.035*(100.-hu(mon))**(1./3.)
      endif

c      n=1   
      call solar(ju,so0,naa/2)
!     since for the old etharg the solar angle was constant
!     in read_par_parser.f all solar sinlat,coslat & tanlat are set equal
!     for flgevp2=2.0
!     
      do n=1,naa
        if(tempv(n).ge.-5.)then
!         rh(n)=amin1(100.00,rh(n))
!         akt = 0.035*(100.00-rh(n))**.333
          petn=0.0075*so0*sqrt((diff(mon)*1.8))*(tempv(n)*1.8+32.)
     *                                     *akt/24.
        else
          petn=0.
        endif
        do ii=1,classcount
          pet(n,ii)=petn
        end do
      end do

c      write(590,59000)pet(nnprint,3),so0,diff(mon),tempv(nnprint),akt
c59000 format(99f12.3)      

      RETURN

      END SUBROUTINE etharg

!***********************************************************************
      SUBROUTINE solar(ju,so0,n)
!***********************************************************************

! THIS SUBROUTINE CALCULATES THE EXTRATERRESTRIAL SOLAR RADIATION
! TAKEN FROM THE HANDBOOK OF HYDROLOGY - CHAPTER 4

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90

! ju is the Julian day
! phi is the latitude (+ve in the N. hemisphere)  << not used?
! ws is the sunset hour angle (radians)
! dr is the relative distance between the earth and the sun
! so is the solar radiation in equivalent water evaporation (mm/d)
!    so used to be in area_watflood but was renamed to so0
! dlta is the solar declination (radians)
! pi = re: circles

!***********************************************************************

      USE area_watflood
      implicit none

      real*4 :: akt1,petn,dlta,dr,sindlta,cosdlta,tandlta,ws0,so0
      integer  :: mon,ju,n,ii

      dlta=0.4093*sin(0.0172142*float(ju)-1.405)
      dr=1.0+0.033*cos(0.0172142*float(ju))
      sindlta=sin(dlta)
      cosdlta=cos(dlta)
      tandlta=sindlta/cosdlta
      ws0=acos(-tanlat(n)*tandlta)
      so0=15.392*dr*(ws0*sinlat(n)*sindlta+coslat(n)*cosdlta*sin(ws0))
!     The functions were subsripted and renamed to accomodata variable solar angles      

c	  write(811,9091)n,dlta,dr,sindlta,cosdlta,tandlta,
c     *             sinlat(n),coslat(n),tanlat(n),ws0,so0
c	  write(811,9091)n,tandlta,tanlat(n),-tandlta*tanlat(n),ws0
c9091    format(i6,12f12.4)

      RETURN

      END SUBROUTINE solar




