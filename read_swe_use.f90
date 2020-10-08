subroutine read_swe_use    

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
    use area_debug
    implicit none
    
    
    character(len=20) :: fname
    character(len=20) :: structname
    character(len=256):: line
    logical :: mustread
    type(XML_PARSE) :: info
    character(len=80) :: tag
    logical :: endtag
    character(len=80), dimension(1:2,1:20) :: attribs
    integer :: no_attribs
    character(len=200), dimension(1:100) :: data
    integer :: no_data
    integer :: i
    
!   NK    
    logical                  :: ign_whitespace
    character(2)             :: units
    character(6)             :: parameterID
    
! get name of XML file (store in fname)
    fname='snow1\SWE_use.xml'
   
    mustread = .true.
    call xml_open( info, fname, mustread )
!
! Check for errors
!
    if ( xml_error(info) ) then
        write(98,*)'Error opening ',fname
        write(98,*)info
        write(98,*)'Program aborted in read_swe_use @ 43'
        print*,'Error opening ',fname
        print*,info
        stop 'Program aborted in read_swe_use @ 43'
    else
!
! Start reading the file
!
        call xml_options( info, ignore_whitespace = ign_whitespace )
        do
            call xml_get( info, tag, endtag, attribs, no_attribs, data, &
            no_data )
            if ( xml_error(info) ) then
                write(98,*)'Error reading ',fname
                write(98,*)info
                write(98,*)'Program aborted in read_swe_use @ 84'
                print*,'Error reading ',fname
                print*,info
                stop 'Program aborted in read_swe_use @ 84'
            endif
! Just write out the tag, its attributes and the character
!data we found
            
            if(tag(1:9).eq.'boolValue')then
!               write( 20,*)'x', (i, '>',trim(data(i)), '<', i=1,no_data)
!               write(20,*)trim(data(1))
               if(trim(data(1)).eq.'true')then
                   use_swe_update=.true.
                   write(98,*)'Info: SWE update = TRUE'
               else
                   use_swe_update=.false.
                   write(98,*)'Info: SWE update = FALSE'
               endif
               write(20,*)use_swe_update
               exit
            endif
                
            if ( .not. xml_ok(info) ) exit
        enddo
    endif
    call xml_close( info )
    if(debug_output)print*,'FEWS SWE update =',use_swe_update

    return
    
end  subroutine read_swe_use