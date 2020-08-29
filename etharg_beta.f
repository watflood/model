      SUBROUTINE etharg_beta(mon,ju)
      
!***********************************************************************
!    Copyright (C) 2004 by Todd Neff and Nicholas Kouwen  
        
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

! THIS SUBROUTINE CALCULATES HARGREAVES PET FOR EACH ELEMENT
! HARGREAVES EQUATION FROM A HANDBOOK OF HYDROLOGY
! - NOTE THAT THERE ARE SOME CHANGES - SEE ORIGINAL PAPER
!   REFERENCED IN MASc THESIS T.Neff 1996

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90
!     rev. 9.9.03  Dec.  15/13  - NK: Change to gridded latitude for etharg
!     rev. 9.9.04  Dec.  17/13  - NK: Change over to gridded clamate normals to diff

! akt - reduction due to excessive humidity (>54%)

!***********************************************************************

      USE area_watflood
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      real*4   :: akt1,petn,dlta,sindlta,cosdlta,tandlta,dr
      integer  :: mon,ju,n,ii,ju_last,iAllocateStatus
      logical  :: firstpass
      
      real*4, dimension(:), allocatable :: so,ws
      
      
      data ju_last/-999/
      data firstpass/.true./

c      if(firstpass)THEN
c        write(811,81100)'    n        dlta          dr     sindlta  
c     *     cosdlta     tandlta    sinlat(n)    coslat(n)    tanlat(n)   
c     *ws(n)       so(n)'
c      ENDIF
c81100 FORMAT(A126)

      if(.not.allocated(so))then
        allocate(so(na),ws(na),stat=iAllocateStatus)
        if (iAllocateStatus .ne. 0) STOP 
     *    '**Allocation failed for dr in etharg @ 32'
      endif

      if(hu(mon).lt.54)then
         akt = .125
      else
         akt = 0.035*(100.-hu(mon))**(1./3.)
      endif

c      call solar(ju)
! ju is the Julian day
! phi is the latitude (+ve in the N. hemisphere)  << not used?
! dlta is the solar declination (radians)
! dr is the relative distance between the earth and the sun
! ws is the sunset hour angle (radians)
! so is the solar radiation in equivalent water evaporation (mm/d)
! pi = re: circles
!     only need to do this once a day (NK)
      if(ju.ne.ju_last)then
        dlta=0.4093*sin(0.0172142*float(ju)-1.405)
        dr=1.0+0.033*cos(0.0172142*float(ju))
        sindlta=sin(dlta)
        cosdlta=cos(dlta)
        tandlta=sindlta/cosdlta
 
 !      Fix tanlat(n),coslat(n),sinlat(n)
 !      to get a value for each row in the map file
 !      can be done in the first pass
        
        do n=1,naa
          ws(n)=acos(-tanlat(n)*tandlta)
          so(n)=15.392*dr*
     *           (ws(n)*sinlat(n)*sindlta+coslat(n)*cosdlta*sin(ws(n)))
        end do
c        print*,ju
      endif
      
      
c      n=nnprint
c	  write(811,9091)n,dlta,dr,sindlta,cosdlta,tandlta,
c     *             sinlat(n),coslat(n),tanlat(n),ws(n),so(n)
c9091    format(i6,12f12.4)


      do n=1,naa
        if(tempv(n).ge.-5.)then
!         SB e-mail Jan. 03/14
!         Hargreaves & Samani 1985
!         ET0 = 0.0023(Tmax - Tmin)1/2(Tavg + 17.8)Ra 
!         Where ET is in mm/d, temps are in Celsius and radiation in mm/d. 

          if(dlyflg)then
            petn=0.0023*so(n)*sqrt(dly_diff(n))*(tempv(n)+17.8)/24.
          else
            petn=0.0023*so(n)*sqrt(diff(mon)*1.8)*(tempv(n)+17.8)/24.
          endif

c         if(n.eq.nnprint)print*,'h&S 1985',day_now,dly_diff(n)
c        if(n.eq.nnprint)print*,ju,mon,diff(mon),mean_dly_diff(n,mon),petn
        else
          petn=0.
        endif
c        if(n.eq.nnprint)print*,ju,ws(n),so(n),petn,tempv(n)
        do ii=1,classcount
          pet(n,ii)=petn
        end do
      end do
      
      ju_last=ju
      firstpass=.false.
     
      RETURN

      END SUBROUTINE etharg_beta

!***********************************************************************
c      SUBROUTINE solar(ju)
!***********************************************************************

! THIS SUBROUTINE CALCULATES THE EXTRATERRESTRIAL SOLAR RADIATION
! TAKEN FROM THE HANDBOOK OF HYDROLOGY - CHAPTER 4

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90

! ju is the Julian day
! phi is the latitude (+ve in the N. hemisphere)  << not used?
! ws is the sunset hour angle (radians)
! dr is the relative distance between the earth and the sun
! so is the solar radiation in equivalent water evaporation (mm/d)
! dlta is the solar declination (radians)
! pi = re: circles

!***********************************************************************

c      USE area_watflood
c      implicit none

c      real*4 :: akt1,petn,dlta,dr,sindlta,cosdlta,tandlta,ws
c      integer  :: mon,ju,n,ii

c      dlta=0.4093*sin(0.0172142*float(ju)-1.405)
c      dr=1.0+0.033*cos(0.0172142*float(ju))
c      sindlta=sin(dlta)
c      cosdlta=cos(dlta)
c      tandlta=sindlta/cosdlta
c      ws=acos(-tanlat(n)*tandlta)
c      so=15.392*dr*(ws*sinlat*sindlta+coslat*cosdlta*sin(ws))

!	write(91,9091)dlta,dr,sindlta,cosdlta,tandlta,
!     *             sinlat,coslat,tanlat,ws,so
!9091  format(12f12.4)

c      RETURN

c      END SUBROUTINE solar




