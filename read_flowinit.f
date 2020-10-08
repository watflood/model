
      SUBROUTINE read_flowinit(flnNum)

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen and dave Watson
        
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
     
!     borrowed swe reader for this s/r
!     nk  Oct. 9/06


!***********************************************************************
!  based on:  READ_SWE_EF -  written Mar/06 by Dave Watson, CHC
!     - Derived from rdswe written by Nick Kouwen
!     - This subroutine reads the ensim compatible gridded SWE file (r2c format)
!***********************************************************************
   
      use area_watflood

! R2C data module
	USE EF_Module

	implicit none
	TYPE(SWEParam) :: header

! SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      integer       :: nclasses,ios,iallocate,iDeallocate
      real*4        :: deffactor

      LOGICAL exists
 
! parameter type definitions
	integer*4 unitNum, flnNum, iStat
	integer*4 n,i,j

! Local variables
	character*4096 line, subString, tmpString
	character*128 keyword, value
	integer lineLen, keyLen, wordCount
	logical rStat, lineType, foundEndHeader
	integer xCountLoc, yCountLoc, attCountLoc
	integer ai, xi, yi, vi, error

! Set unit and fln number
c	unitNum = 99
c	flnNum = 99
 
! Open the file
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	IF(exists)then
		open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
 		if(ios.ne.0)then
		  write(98,*)'Error: Problems opening ',fln(flnNum)(1:40)
		  STOP ' Stopped in read_flowinit @ 52_ef'
		endif
      else
          print*
		write(98,*)'Error: the file: ',	fln(flnNum)(1:40)
          write(98,*)'Error: Unit # =',unitNum
		write(98,*)'Error:is NOT found'
		STOP ' Program STOPPED in read_flowinit @ 57'
	endif
	write(98,*)'Info: Opened ',fln(flnnum)(1:40)

! Initialize default values
	CALL InitSWEParam(header)	

! Search for and read r2c file header
	line(1:1) = '#'
	foundEndHeader = .false.

	do WHILE((.NOT.foundEndHeader) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(LEN_TRIM(line) .eq. 0))) 	

		read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
		if(ios .eq. -1)then
			write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
			STOP ' Stopped in read_swe_ef'
		end if

		rStat = Detab(line)				! replace tabs with spaces
		line = ADJUSTL(line)		! Get rid of leading white space
		lineLen = LEN_TRIM(line)		! Find the length excluding trailing spaces

		if(line(1:1) .eq. ':')then
			wordCount = SplitLine(line, keyword, subString)	! find the keyword
			rStat = ToLowerCase(keyword)
			KeyLen = LEN_TRIM(keyword)

			if(keyword(1:KeyLen) .eq. ':endheader')then
				foundEndHeader = .TRUE.

			else
				iStat = ParseSWEParam(header,keyword,keyLen,
     &													subString)
				if(iStat .lt. 0) then
					write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
					write(*,'(2(A))') '   in line: ',line					
					STOP ' Stopped in read_swe_ef'
					return
				else if(iStat .eq. 0) then
!					write(*,'((A), (A))')  'Unrecognized keyword line: ',
!     &										line
				endif
			end if
		end if
	end do
!***************************************
!	Finished reading the header
!***************************************

! Assign the parsed parameters to the model variables		
	xCountLoc = header%r2cp%xCount
	yCountLoc = header%r2cp%yCount
	attCountLoc =header%r2cp%ep%attCount
c	imax=header%r2cp%yCount
c	jmax=header%r2cp%xCount 
	deffactor = header%initHeatDeficit !Set InitHeatDeficit


! Read the data section
	CALL LoadAttributeData(header%r2cp%ep, xCountLoc,
     &						yCountLoc, unitNum)	
	
      
!      write(63,*)xCountLoc,yCountLoc,unitNum
      
      	
! Validate parameters
	if(abs(header%r2cp%xOrigin-xorigin).gt.0.001)
     *    write(63,*)'xorig_rte.ne.xorigin'
	if(abs(header%r2cp%yOrigin-yorigin).gt.0.001)
     *    write(63,*)'yorig_rte.ne.yorigin'
      if(abs(header%r2cp%xDelta-xdelta).gt.0.001)
     * 	write(63,*)'xdelta_rte.ne.xdelta'
	if(abs(header%r2cp%yDelta-ydelta).gt.0.001)
     * 	write(63,*)'ydelta_rte.ne.ydelta'
      if(header%r2cp%xCount.ne.xcount)
     * 	      write(63,*)'xcount_rte.ne.xcount'
      if(header%r2cp%yCount.ne.ycount)
     * 	      write(63,*)'ycount_rte.ne.ycount'
      if(abs(header%r2cp%xOrigin-xorigin).gt.0.001.or.
     &	abs(header%r2cp%yOrigin-yorigin).gt.0.001.or.
     &	abs(header%r2cp%xDelta-xdelta).gt.0.001.or.
     &	abs(header%r2cp%yDelta-ydelta).gt.0.001.or.
     &	header%r2cp%xCount.ne.xcount.or.
     &	header%r2cp%yCount.ne.ycount) then
            print*
	      write(63,*)'Mismatch between ',fln(flnNum)
            write(63,*)'    and SHD files'
            write(63,*)'Check files for origins, deltas and counts'
            write(63,*)'Could be due to # significant digits in header' 
		STOP 'Program aborted in read_flowinit_ef @ 159'
	endif

	nclasses = header%r2cp%ep%attCount

c         imax=ycount
c     	  jmax=xcount

!     check added Jul. 13/11 nk
      if(yCountLoc.ne.ycount)then
        write(63,*)'ycount in ',fln(flnNum)
        write(63,*)'does not match the ycount in the shd file'
        stop 'Program aborted in read_flowinit.f @ 165'
      endif
      if(xCountLoc.ne.xcount)then
        write(63,*)'xcount in ',fln(flnNum)
        write(63,*)'does not match the xcount in the shd file'
        stop 'Program aborted in read_flowinit.f @ 171'
      endif

      if(allocated(inarray))then
        deallocate(inarray)
        write(63,*)'deallocated inarray in read_flowinit'
      endif
      if(.NOT.allocated(inarray))then
        allocate(inarray(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0)then
          write(63,*)'Error with allocation of inarray in read_flowinit'    
          stop 'Program aborted in read_flowinit @ 178'
        endif
      endif

! Copy Attribute data over to Nick's global variables
	do ai=1,attCountLoc
		vi = 0
!		do yi=yCountLoc,1,-1      <<<<<<<top to bottom reversed !!!
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  inarray(yi,xi) = header%r2cp%ep%attList(ai)%val(vi)
!	write(63,*)yi,xi,header%r2cp%ep%attList(ai)%val(vi)
			end do
c            write(55,55555)(inarray(yi,xi),xi=1,xCountLoc)          !?????????????????????
		end do
	      write(55,*)
55555       format(999e10.3)
	    do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
	        if(ai.eq.1)then
                  qi1(n)=inarray(i,j)
	            qi2(n)=inarray(i,j)
	        elseif(ai.eq.2)then
	            qo1(n)=inarray(i,j)
	            qo2(n)=inarray(i,j)
      	    elseif(ai.eq.3)then
	            store1(n)=inarray(i,j)
	            store2(n)=inarray(i,j)
	        elseif(ai.eq.4)then
	            over(n)=inarray(i,j)
	        elseif(ai.eq.5)then
	            lzs(n)=inarray(i,j)
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
	        elseif(ai.eq.6)then
	            qold(n)=inarray(i,j)
              endif
          end do
                  
        if(wetflg.eq.'y')then
	    do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
	        if(ai.eq.7)then
	            qiwet1(n)=inarray(i,j)
	            qiwet2(n)=inarray(i,j)
	        elseif(ai.eq.8)then
	            qowet1(n)=inarray(i,j)
	            qowet2(n)=inarray(i,j)
	        elseif(ai.eq.9)then
	            wstore1(n)=inarray(i,j)
	            wstore2(n)=inarray(i,j)
	        elseif(ai.eq.10)then
	            hcha1(n)=inarray(i,j)
	            hcha2(n)=inarray(i,j)
	        elseif(ai.eq.11)then
	            hwet1(n)=inarray(i,j)
	            hwet2(n)=inarray(i,j)
	        endif
          end do
        endif
      end do

! Deallocate the attribute data now that global attributes have been set
	do ai=1,attCountLoc
		deallocate ( header%r2cp%ep%attList(ai)%val, STAT = error )
		if (error.ne.0) STOP 'deallocation error in read_gsm_ef()' 
	end do

      close (unit=unitNum)
      
      deallocate(inarray,stat=iDeallocate)
      if(iDeallocate.ne.0)then
        write(63,*)'Error with deallocation of p array in read_flowinit'    
d       stop 'Program aborted in read_flowinit @ 247'
      endif
      
	write(98,*)'Info: Finished reading ',fln(flnNum)(1:40)

      RETURN

      END SUBROUTINE read_flowinit



