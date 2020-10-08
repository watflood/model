      SUBROUTINE read_sweinit(unitNum,flnNum)
	      
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
!    GNU Lesser General Public License for more detiils.

!    You should have received a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
     
!***********************************************************************
!  READ_SWE_EF -  written Mar/06 by Dave Watson, CHC
!     - Derived from rdswe written by Nick Kouwen
!     - This subroutine reads the ensim compatible gridded SWE file (r2c format)
!***********************************************************************
   
!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001

      use area_watflood

! R2C data module
	USE EF_Module

	implicit none
	TYPE(SWEParam) :: header

! SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      integer       :: nclasses,ios,i,j,ii,n,iallocate
      real*4        :: deffactor,class_sum

      LOGICAL exists
 
! parameter type definitions
	integer*4 unitNum, flnNum, iStat

! Local variables
!      chracter*1  firstpass
	character*4096 line, subString, tmpString
	character*128 keyword, value
	integer lineLen, keyLen, wordCount
	logical rStat, lineType, foundEndHeader
	integer xCountLoc, yCountLoc, attCountLoc
	integer xi, yi, vi, error

      real*4,  dimension(:,:),allocatable :: ratio2

      allocate(ratio2(ycount,xcount),stat=iAllocate)
      if(iAllocate.ne.0)then
         print*,'Error with allocation of ratio2 in read_sweinit.f @ 63'
         STOP 'Program aborted in sub @ 595'
      endif
      
!      data /firstpass/'y'/

! Set unit and fln number
c	unitNum = 99
c	flnNum = 36   ! gridded swe.r2c file
 
! Open the file
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	IF(exists)then
!       sweinit file will be read even in a run with resume = 'y'
	  open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
 	  if(ios.ne.0)then
	    print*,'Problems opening ',fln(flnNum)(1:40)
	    print*
	    STOP ' Stopped in read_sweinit_ef @ 53'
	  endif
	else
	  if(resumflg.ne.'y')then
	    print*,'ERROR: the SWE file: ',	fln(flnNum)(1:40)
	    print*,'is NOT found'
	    STOP ' Program STOPPED in read_sweinit @ 59'
	  else
!         this means reading the swe from the soilinit.r2c file
	    return
	  endif
	endif


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
!				write(*,'((A), (A))')  'Unrecognized keyword line: ',
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
            print*,'       and SHD files'
            print*,'Check files for origins, deltas and counts'
            print*,'Could be due to # significant digits in header' 
		STOP 'Program aborted in read_swe_ef @ 141'
	endif

	nclasses = header%r2cp%ep%attCount
	if(nclasses.lt.classcount)then
!         there are fewer classes in the swe file than needed
		print*
		print*,' The swe file has ',nclasses,' land cover classes'
		print*,' The shd file stipulates ',classcount,' which'
		print*,' includes the impervious class'
		print*,' Please fix the swe file for proper # of classes'
		print*
          STOP 'Program aborted in read_swe_ef'
      endif

!     rev. 10.3.07 Mar.  04/20  = NK Fixed weighted swe in wfo file for grids with water
!     rev. 10.3.08 Mar.  06/20  = NK For FEWS, add snow1\swe.nc for swe updating
!     Start by finding the weighted swe for each grid   
      do n=1,naa
          totsnw(n)=0.0
          class_sum=0.0
          do ii=1,classcount
             if(aclass(n,ii).gt.0.0)then
                if(ii.ne.ii_water)then
                  class_sum=class_sum+aclass(n,ii)
                  totsnw(n)=totsnw(n)+snowc(n,ii)*aclass(n,ii)*sca(n,ii)
                endif
             endif
          end do
          totsnw(n)=totsnw(n)/class_sum
      end do     

c        imax=ycount
c	  jmax=xcount

!     check added Jul. 13/11 nk
      if(yCountLoc.ne.ycount)then
        print*,'ycount in ',fln(flnNum)
        print*,'does not match the ycount in the shd file'
        print*
        stop 'Program aborted in read_swe_ef.f @ 165'
      endif
      if(xCountLoc.ne.xcount)then
        print*,'xcount in ',fln(flnNum)
        print*,'does not match the xcount in the shd file'
        print*
        stop 'Program aborted in read_swe_ef.f @ 171'
      endif
      if(attCountLoc.ne.classcount)then
        print*,'classcount in ',fln(flnNum)(1:50)
        print*,'does not match the classcount in the shd file'
        print*
        stop 'Program aborted in read_sweinit.f @ 187'
      endif
      
! Copy Attribute data over to Nick's global variables
	do ii=1,attCountLoc
		vi = 0
!		do yi=yCountLoc,1,-1   !  backwards  nk  Mar. 20/07
		do yi=1,yCountLoc
			do xi=1,xCountLoc
!                 put the same swe in all classes                  
				vi = vi+1
			    snw(yi,xi) = header%r2cp%ep%attList(ii)%val(vi)
			end do
          end do
          end do

!     rev. 10.3.09 Mar.  07/20  = NK Revise swe updating to maintain relative swe in classes
!     Find the ratio of obs swe / weighted model swe
      do n=1,naa
			i=yyy(n)
			j=xxx(n)
              if(totsnw(n).gt.10.0)then
!                  If there is not much snow in the model, ratio is too erratic   
!                  and just use the snow course data                  
                   ratio2(i,j)=snw(i,j)/totsnw(n)
              else
                   ratio2(i,j)=1.0000
              endif
c      if(i.eq.ycount/2.and.j.eq.xcount/2)
c     *      write(666,*)id,totaltime,snw(i,j),totsnw(n),
c     *                ratio2(ycount/2,xcount/2)
      end do
      
!   Adjust and put INTO VECTOR FORMAT
	do ii=1,attCountLoc
		do n=1,naa
			i=yyy(n)
			j=xxx(n)
              
            if(id.eq.1)then
!             Old way: swe in all classes replaced by swe from the r2c file              
			snowc(n,ii)=max(0.0,snw(i,j))
            else
!     rev. 10.3.09 Mar.  07/20  = NK Revise swe updating to maintain relative swe in classes
!             Multiply by the ratio so relative swe in classes stay the same
              snowc(n,ii)=snowc(n,ii)*ratio2(i,j)
            endif
			if(snowc(n,ii).gt.0.0)then
				sca(n,ii)=1.0
				oldsca(n,ii)=1.0
!				set the initial heat deficit for the snow:
				def(n,ii)=deffactor*snowc(n,ii)
			else
!     May 10, 2002 Added these lines to avoid underflows later! AB
				sca(n,ii)=0.0
				oldsca(n,ii)=0.0
				def(n,ii)=0.0
              endif
		end do
      end do

! Deallocate the attribute data now that global attributes have been set
	do ii=1,attCountLoc
		deallocate ( header%r2cp%ep%attList(ii)%val, STAT = error )
		if (error.ne.0) STOP 'deallocation error in read_sweinit()' 
      end do
      
      deallocate(ratio2,STAT = error)
	if (error.ne.0) STOP 'deallocation error in read_sweinit()' 

      close (unit=unitNum)

	if(iopt.ge.1)print*,'Finished reading ',fln(flnNum)(1:30)

      RETURN

      END SUBROUTINE read_sweinit



