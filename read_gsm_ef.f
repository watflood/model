      SUBROUTINE read_gsm_ef()

!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen and Dave Watson (NRC)
        
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
!  READ_GSM_EF - written Mar/06 by Dave Watson, CHC
!	- Derived from read_gsm written by Nick Kouwen.
!	- This subroutine reads the ensim compatible gridded soil moisture
!     (GSM) file (r2c format)
!***********************************************************************

!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001

      use area_watflood

! R2C data module
	USE EF_Module

	implicit none
	TYPE(GSMPARAM) :: header


! SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      logical       :: exists
	integer i,j,n,ios,iAllocate,iDeallocate,classcount_local,ii
	real*4 conv

! Parameter type definitions
	integer*4 unitNum, flnNum, iStat

! Local variables
      character*1 firstpass
	character*128 keyword, value
	character*4096 line, subString, tmpString
	integer lineLen, keyLen, wordCount
	logical rStat, lineType, foundEndHeader
	integer xCountLoc, yCountLoc, attCountLoc
	integer ai, xi, yi, vi, error

      data firstpass/'y'/

! Set unit and fln number
	unitNum = 287
	flnNum = 37

! Open the file
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	if(exists)then
!       gsm file will be read even in a run with resume = 'y'
	  open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
	    if(ios.ne.0)then
	  	print*,'Problems opening ',fln(flnNum)
	 	print*
		STOP ' Stopped in read_gsm_ef'
	  endif
	else
	  if(resumflg.ne.'y')then
          print*,'ERROR: the gridded soil moisture (GSM) file ',
     *			fln(flnNum)
	    print*,'is NOT found'
	    STOP ' Program STOPPED in read_gsm_ef'
	  else
!         this means using the swe from the soilinit.r2c file
	    return
	  endif
	endif

! Initialize default values
	CALL InitGSMParam(header)	

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
			STOP ' Stopped in read_gsm_ef'
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
				iStat = ParseGSMParam(header,keyword,keyLen,
     &													subString)
				if(iStat .lt. 0) then
					write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
					write(*,'(2(A))') '   in line: ',line					
					STOP ' Stopped in read_gsm_ef'
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
	conv = header%r2cp%unitConv 
	classcount_local = header%r2cp%ep%attCount

! Read the data section
	CALL LoadAttributeData(header%r2cp%ep, xCountLoc,
     &						yCountLoc, unitNum)

!     This section copied from read_swe_ef by NK  Jul. 5/07
! Validate parameters
!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001
	if(abs(header%r2cp%xOrigin-xorigin).gt.0.001)
     *    print*,'xorig_swe.ne.xorigin'
	if(abs(header%r2cp%yOrigin-yorigin).gt.0.001)
     *    print*,'yorig_swe.ne.yorigin'
      if(abs(header%r2cp%xDelta-xdelta).gt.0.001)
     * 	print*,'xdelta_swe.ne.xdelta'
	if(abs(header%r2cp%yDelta-ydelta).gt.0.001)
     * 	print*,'ydelta_swe.ne.ydelta'
      if(header%r2cp%xCount.ne.xcount)
     * 	print*,'xcount_swe.ne.xcount'
      if(header%r2cp%yCount.ne.ycount)
     * 	print*,'ycount_swe.ne.ycount'
      if(abs(header%r2cp%xOrigin-xorigin).gt.0.001.or.
     &	abs(header%r2cp%yOrigin-yorigin).gt.0.001.or.
     &	abs(header%r2cp%xDelta-xdelta).gt.0.001.or.
     &	abs(header%r2cp%yDelta-ydelta).gt.0.001.or.
     &	header%r2cp%xCount.ne.xcount.or.
     &	header%r2cp%yCount.ne.ycount) then
	      PRINT*,'Mismatch between ',fln(flnNum) 
            print*,'     and SHD files'
            print*,'Check files for origins, deltas and counts'
            print*,'Could be due to # significant digits in header' 
		STOP 'Program aborted in read_gsm_ef @ 141'
	endif

! We have finished reading the file so close it
	close(unit=unitNum,status='keep',iostat=ios)

!	if(classcount_local.ne.classcount)then
	if(classcount_local.ne.classcount)then     ! <<<<<<<<<<<<<<nk sep 28/06
		print*,'No of classes in the gridded soil moisture file'
		print*,'do not match the number of classes in the shed'
		print*,'file.'
		print*,'no of classes in the shd file =',classcount
        print*,'Including impervious'
		print*,'no of classes in the gsm file =',classcount_local
		print*,'All classes assumed to have same init soil moisture'
		print*
		pause 'Hit enter to continue or ctrl C to abort.'	      
	endif

! Allocate the ssmc array - deallocate if one already exists
	if(ALLOCATED(ssmc)) then
		deallocate(ssmc,stat=iDeallocate)
		if(iDeallocate.ne.0) STOP
     *		'ERROR: Failed to deallocate ssmc array in read_gsm_ef'
	endif
		

      if(firstpass.eq.'y')then
	  allocate(api(na,classcount),stat=iAllocate)
	  if(iAllocate.ne.0) STOP
     *   'ERROR: Failed to allocate ssmc array in read_gsm_ef'
      endif

	allocate(ssmc(ycountloc,xcountloc),stat=iAllocate)
	  if(iAllocate.ne.0) STOP
     *   'ERROR: Failed to allocate ssmc array in read_gsm_ef'

!      print*,'ssmc allocated for ',ycountloc,xcountloc
!	pause

!     check added Jul. 17/12 nk
      if(yCountLoc.ne.ycount)then
        print*,'ycount in ',fln(flnNum)
        print*,'does not match the ycount in the shd file'
        print*
        stop 'Program aborted in read_gsm_ef.f @ 199'
      endif
      if(xCountLoc.ne.xcount)then
        print*,'xcount in ',fln(flnNum)
        print*,'does not match the xcount in the shd file'
        print*
        stop 'Program aborted in read_gsm_ef.f @ 205'
      endif
      if(attCountLoc.ne.classcount)then
        print*,'classcount in ',fln(flnNum)
        print*,'does not match the classcount in the shd file'
        print*
        stop 'Program aborted in read_gsm_ef.f @ 211'
      endif

! Copy Attribute data over to Nick's global variables
	do ai=1,attCountLoc
	  vi = 0
!	    do yi=yCountLoc,1,-1   !  backwards  nk  Mar. 20/07
	    do yi=1,yCountLoc
		do xi=1,xCountLoc
	        vi = vi+1
	        ssmc(yi,xi) =  header%r2cp%ep%attList(ai)%val(vi)*conv
!             ................... 
!             added Jan. 21/07  nk
!             need to vectorize ssmc here to save memory
!             just one 2D array needed this way
              ssmc(yi,xi)=amin1(0.55,ssmc(yi,xi))
              ssmc(yi,xi)=amax1(0.001,ssmc(yi,xi))
              n=s(yi,xi)
              if(n.gt.0)then     ! i.e. we're in the watershed
!               api here is the initial soil moisture
!               but it will be the antecedent precip index in runoff
!               It is converted in soilinit
	          api(n,ai)=ssmc(yi,xi)
	        endif
!             ................... 
	    end do
        end do
	end do

! Deallocate the attribute data now that the global variable have been set
	do ai=1,attCountLoc
	  deallocate ( header%r2cp%ep%attList(ai)%val, STAT = error )
	  if (error.ne.0) STOP 'deallocation error %val in read_gsm_ef()' 
	end do

 	deallocate (ssmc,STAT=error )
	if (error.ne.0) STOP 'deallocation error ssmc in read_gsm_ef()' 

      firstpass='n'

	if(iopt.ge.1)print*,'Finished reading ',fln(flnNum)(1:30)

      return
  
      end SUBROUTINE read_gsm_ef










