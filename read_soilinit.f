      SUBROUTINE read_soilinit()
	      
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
     
!     borrowed swe reader for this s/r
!     nk  Oct. 9/06

!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001
!     rev. 9.8.35  Oct.  23/12  - NK: Fixed bug in read_soilinit_ef

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

      integer       :: classcount_local,ios
      real*4        :: deffactor

      LOGICAL exists
 
! parameter type definitions
	integer*4 unitNum, flnNum, iStat

! Local variables
	character*4096 line, subString, tmpString
	character*128 keyword, value
	integer lineLen, keyLen, wordCount
	logical rStat, lineType, foundEndHeader
	integer xCountLoc,yCountLoc,attCountLoc
	integer ai,xi,yi,vi,error,i,j,n,ii

! Set unit and fln number
	unitNum = 99
	flnNum = 99
 
! Open the file
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	IF(exists)then
		open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
 		if(ios.ne.0)then
			print*,'Problems opening ',fln(flnNum)(1:40)
			print*
			STOP ' Stopped in read_soil_init @ 53'
		endif
	else
		print*,'ERROR: the soil_init file: ',fln(flnNum)(1:40)
		print*,'is NOT found'
		STOP ' Program STOPPED in read_soilinit @ 58'
	endif


! Initialize default values
	CALL InitSWEParam(header)	

! Search for and read r2c file header
	line(1:1) = '#'
	foundEndHeader = .false.
      write(63,*)'Reading the header in read_soilinit'
	do WHILE((.NOT.foundEndHeader) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(LEN_TRIM(line) .eq. 0))) 	

		read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
		
d		write(63,*)line(1:72)
		
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
C					write(*,'((A), (A))')  'Unrecognized keyword line: ',
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
	
      
!      print*,xCountLoc,yCountLoc,unitNum
     	
! Validate parameters
!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001
	if(abs(header%r2cp%xOrigin-xorigin).gt.0.001)
     *        print*,'xorig_swe.ne.xorigin'
	if(abs(header%r2cp%yOrigin-yorigin).gt.0.001)
     *         print*,'yorig_swe.ne.yorigin'
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
            print*,'    and SHD files'
            print*,'Check files for origins, deltas and counts'
            print*,'Could be due to # significant digits in header' 
		STOP 'Program aborted in read_soilinit_ef @ 145'
	endif

c	classcount_local = header%r2cp%ep%attCount
	attCountLoc = header%r2cp%ep%attCount
d	print*,'attribute count=',attCountLoc

!     the number of state variables = 10  < change if needed
	classcount_local=(attCountLoc-3)/10   !changed fron 9 -> 10 Mar15/11 nk

d	print*,'header%r2cp%ep%attCount=',header%r2cp%ep%attCount
d	print*,'classcount_local=',classcount_local
d	print*
!	pause 'in read_soilinit_ef.for @ 179'

c        imax=ycount
c	  jmax=xcount

!     check added Jul. 13/11 nk
      if(yCountLoc.ne.ycount)then
        print*,'ycount in ',fln(flnNum)
        print*,'does not match the ycount in the shd file'
        print*
        stop 'Program aborted in read_flowinit.f @ 165'
      endif
      if(xCountLoc.ne.xcount)then
        print*,'xcount in ',fln(flnNum)
        print*,'does not match the xcount in the shd file'
        print*
        stop 'Program aborted in read_flowinit.f @ 171'
      endif

! Copy Attribute data over to Nick's global variables
	do ii=1,classcount_local
	      ai=ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
!	              print*,yi,xi,header%r2cp%ep%attList(ai)%val(vi)
			end do
d             if(iopt.ge.2)write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
55555       format(999e10.3)

	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              v(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

	do ii=1,classcount_local
	      ai=classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
!	              print*,yi,xi,header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              d1(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

	do ii=1,classcount_local
	      ai=2*classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              d1fs(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

	do ii=1,classcount_local
	      ai=3*classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              uzs(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

	do ii=1,classcount_local
	      ai=4*classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              uzsfs(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do
	do ii=1,classcount_local
	      ai=5*classcount_local+ii
		vi = 0
d         write(63,*)'snowc ',ii
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              snowc(n,ii)=p(i,j)
            end do
      end do
	do ii=1,classcount_local
	      ai=6*classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              sca(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

	do ii=1,classcount_local
	      ai=7*classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              wcl(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

	do ii=1,classcount_local
	      ai=8*classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              def(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

	do ii=1,classcount_local
	      ai=9*classcount_local+ii
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
			  vi = vi+1
			  p(yi,xi)=header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(63,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(63,*)
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
              api(n,ii)=p(i,j)
!             print*,n,qi1(n)
            end do
      end do

!             SET INITIAL UNSATURATED ZONE EFFECTIVE POROSITY:
!             THE EFFECTIVE POROSITY DEPENDS ON SOIL MOISTURE
!             WHICH IS ENTERED FOR EACH GRID
!             effpor is for the intermediate zone so same for bare & sca
      do n=1,naa
	  do ii=1,classcount
              effpor(n,ii)=spore(ii)-api(n,ii)
              effpor(n,ii)=amax1(0.0001,effpor(n,ii))
              effpor(n,ii)=amin1(spore(ii),effpor(n,ii))
              if(api(n,ii).lt.0.0)then
                api(n,ii)=0.0
              endif
        end do
	end do

! Copy Attribute data over to Nick's global variables
	do ai=1,attCountLoc
		vi = 0
		do yi=1,yCountLoc
			do xi=1,xCountLoc
				vi = vi+1
			 p(yi,xi) = header%r2cp%ep%attList(ai)%val(vi)
			end do
d                  write(55,55555)(p(yi,xi),xi=1,xCountLoc)
		end do
d	      write(55,*)
!     rev. 9.8.35  Oct.  23/12  - NK: Fixed bug in read_soilinit_ef
	      do n=1,naa
	        i=yyy(n)
	        j=xxx(n)
c	        if(ai.eq.61)then
	        if(ai.eq.10*classcount_local+1)then
                  tto(n)=p(i,j)
	        elseif(ai.eq.62)then
	        elseif(10*classcount_local+2.eq.62)then
	            ttomin(n)=p(i,j)
c      	      elseif(ai.eq.63)then
      	      elseif(10*classcount_local+3.eq.63)then
	            ttomax(n)=p(i,j)
	        endif
            end do
      end do

      do n=1,naa
	  do ii=1,classcount
              effpor(n,ii)=spore(ii)-api(n,ii)
              effpor(n,ii)=amax1(0.0001,effpor(n,ii))
              effpor(n,ii)=amin1(spore(ii),effpor(n,ii))
              if(api(n,ii).lt.0.0)then
                api(n,ii)=0.0
              endif
        end do
	end do

! Deallocate the attribute data now that global attributes have been set
	do ai=1,attCountLoc
		deallocate ( header%r2cp%ep%attList(ai)%val, STAT = error )
		if (error.ne.0) STOP 'deallocation error in read_gsm_ef()' 
	end do

      close (unit=unitNum)

	if(iopt.ge.1)then
	  print*,'Finished reading ',fln(flnNum)(1:20)
	  print*,'Earlier initialization overwritten'
	endif

      RETURN

      END SUBROUTINE read_soilinit



