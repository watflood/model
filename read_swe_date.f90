subroutine read_swe_date

!***********************************************************************
!    Copyright (C) 2019 by Nicholas Kouwen  
        
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
     
!---------------------------------------------------------------
! xmlparse.f90 - Simple, limited XML parser in Fortran
    use xmlparse
! $Id: xmlparse.f90,v 1.14 2007/12/07 10:10:19 arjenmarkus Exp $
!
! Arjen Markus
! Adopted & modified for WATFLOOD: Nicholas Kouwen June 2019 
! to read swe and uzs modifiers for FEWS    
!---------------------------------------------------------------
! This s/r opens & reads the adjustment factor for swe in FEWS
    
!     rev. 10.2.55 June  09/19  - NK Added read_swe_update.f90 for read in swe adjustment factors
    use area_watflood
    implicit none
    
    character(len=20) :: fname
    character(len=20) :: structname
    character(len=30):: line,junk(10)
    logical :: mustread
    type(XML_PARSE) :: info
    character(len=80) :: tag
    logical :: endtag
    character(len=80), dimension(1:2,1:20) :: attribs
    integer :: no_attribs
    character(len=200), dimension(1:100) :: data
    integer :: no_data,dayrad(12),dayradL(12)
    integer :: i,Zoffset
    
!   NK    
    logical                  :: ign_whitespace
    character(2)             :: units
    character(6)             :: parameterID

    data dayrad/31,28,31,30,31,30,31,31,30,31,30,31/    
    data dayradL/31,29,31,30,31,30,31,31,30,31,30,31/    
    
! get name of XML file (store in fname)
    fname='snow1\swe_date.xml'
   
    mustread = .true.
    call xml_open( info, fname, mustread )
!
! Check for errors
!
    if ( xml_error(info) ) then
        print*,'Error opening ',fname
        print*,info
        stop 'Program aborted @ 43'
    else
!
! Start reading the file
!
        call xml_options( info, ignore_whitespace = ign_whitespace )
        do
            call xml_get( info, tag, endtag, attribs, no_attribs, data, &
            no_data )
            if ( xml_error(info) ) then
                print*,'Error reading ',fname
                print*,info
                stop 'Program aborted @ 43'
            endif
! Just write out the tag, its attributes and the character
!data we found
            write( 20,*)'x', (i, '>',trim(data(i)), '<', i=1,no_data)
            
            if(tag.eq.'stringValue')then
               write( 20,*)'x', (i, '>',trim(data(i)), '<', i=1,no_data)
!               write(20,*)trim(data(1))
!               write(20,*)data(1)
               line=trim(data(1))
!    integer :: yy_swe_update,mm_swe_update,dd_swe_update,hh_swe_update
!               read(data(1),6000)junk,line,dd_swe_update,hh_swe_update,junk,yy_swe_update
               read(line,*)(junk(i),i=1,6)
!               write(20,*)(junk(i),i=1,6)
               read(junk(3)(1:3),*)dd_swe    !day
               read(junk(4)(1:2),*)hh_swe    !hour
               read(junk(5)(4:6),*)junk(10)
               read(junk(10),*)Zoffset
               read(junk(6)(1:4),*)yy_swe    !year
               
               if(junk(2)(1:3).eq.'Jan')then
                   mm_swe=1
               elseif(junk(2)(1:3).eq.'Feb')then
                   mm_swe=2
               elseif(junk(2)(1:3).eq.'Mar')then
                   mm_swe=3
               elseif(junk(2)(1:3).eq.'Apr')then
                   mm_swe=4
               elseif(junk(2)(1:3).eq.'May')then
                   mm_swe=5
               elseif(junk(2)(1:3).eq.'Jun')then
                   mm_swe=6
               elseif(junk(2)(1:3).eq.'Jul')then
                   mm_swe=7
               elseif(junk(2)(1:3).eq.'Aug')then
                   mm_swe=8
               elseif(junk(2)(1:3).eq.'Sep')then
                   mm_swe=9
               elseif(junk(2)(1:3).eq.'Oct')then
                   mm_swe=10
               elseif(junk(2)(1:3).eq.'Nov')then
                   mm_swe=11
               elseif(junk(2)(1:3).eq.'Dec')then
                   mm_swe=12
               endif
               
               hh_swe=hh_swe-Zoffset  ! note -(-) here for going west
               
               if(mod(yy_swe,4).eq.0)then                                     ! Feb in a leap year
                   if(dd_swe.gt.dayradL(mm_swe))then
                        write(98,*)'Error: date in swe_data.xml file outside bounds'
                        stop 'Error: date in swe_data.xml file outside bounds'
                   endif
               else                                                      ! not a leap year
                  if(dd_swe.gt.dayrad(mm_swe))then
                        write(98,*)'Error: date in swe_data.xml file outside bounds'
                        stop 'Error: date in swe_data.xml file outside bounds'
                  endif
               endif
               
               if(hh_swe.ge.24)then
                    hh_swe=hh_swe-24
                    if(mod(yy_swe,4).eq.0)then                                     ! Feb in a leap year
                       if(mm_swe.eq.2.and.dd_swe.eq.dayradL(2))then           ! Feb. 29
                           dd_swe=1
                           mm_swe=3
                           print*,'a'
                       elseif(mm_swe.eq.12.and.dd_swe.eq.dayradL(12))then     ! Dec. 31
                           dd_swe=1
                           mm_swe=1
                           yy_swe=yy_swe+1
                           print*,'b'
                       elseif(dd_swe.eq.dayrad(mm_swe))then         
                           dd_swe=1
                           mm_swe=mm_swe+1
                           print*,'c'
                       else
                           dd_swe=dd_swe+1
                           print*,'d'
                       endif
                    else                                                      ! not a leap year
                      if(mm_swe.eq.2.and.dd_swe.eq.dayrad(2))then            ! Feb. 28
                           dd_swe=1
                           mm_swe=3
                           print*,'e'
                       elseif(mm_swe.eq.12.and.dd_swe.eq.dayrad(12))then      !Dec. 31
                           dd_swe=1
                           mm_swe=1
                           yy_swe=yy_swe+1
                           print*,'f'
                       elseif(dd_swe.eq.dayrad(mm_swe))then         
                           dd_swe=1
                           print*,'g'
                           mm_swe=mm_swe+1
                       else
                           dd_swe=dd_swe+1
                           print*,'h'
                       endif
                    endif
               endif
               
!               write(20,*) yy_swe,mm_swe,dd_swe,hh_swe,Zoffset   
!               write(20,*) yy_swe,mm_swe,dd_swe,hh_swe,Zoffset   
!               write(20,*) yy_swe,mm_swe,dd_swe,hh_swe,Zoffset  
               exit
            endif
                
            if ( .not. xml_ok(info) ) exit
        enddo
    endif
    call xml_close( info )

    return
    
end subroutine