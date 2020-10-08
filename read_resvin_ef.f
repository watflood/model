      SUBROUTINE read_resvin_ef()

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
     
!*****************************************************************************
!* READ_RESVIN_EF - written Sep/06 by Dave Watson, CHC
!	- Derived from rdresvin written by Nick Kouwen
!	This subroutine reads the measured reservoir inflow (RIN) file 
!	(tb0 format)
!*****************************************************************************

!     rev. 9.5.68  Oct.  07/09  - NK: debugged read_resvin_ef.f
!     rev. 9.5.78  Nov.  04/09  - NK: matched resvin locations to reach numbers
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
!     rev. 9.5.82  Jan.  26/09  - NK: replaced error check for inflow locations
!     rev. 9.7.12  Nov.  10/10  - NK: fix array bugs for reservoir inflows
!     rev. 9.8.23  Aug.  03/12  - NK: Added resinid1flg to use resinflg for id=1
!     rev. 9.8.41  Jan.  28/13  - NK: fixed bug in lst for level print statement
!     rev. 9.8.41  Jan.  31/13  - NK: fixed bug in read_resvin: nopti int conversion

      USE area_watflood

!*TB0 data module
	USE EF_Module

	implicit none
	TYPE(ResvinParam) :: header
	TYPE(ResvinColumnMetaData) :: colHeader

      Integer  :: ios,j,k,i,n,l
	integer  :: nrel_in_max,iAllocate,iDeallocate
      real*4   :: factor

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      LOGICAL exists
	CHARACTER(1)  :: firstpass

      data firstpass/'y'/

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

!*Parameter type definitions
	integer*4 unitNum, flnNum, iStat

!*Local variables
	character*4096 line, subString, tmpString
	character*128 keyword, value
	integer lineLen, keyLen, wordCount
	logical rStat, lineType, foundEndHeader, insideColMetaData
	logical errflg1

! Set unit and fln number
	unitNum = 99
	flnNum = 8

! Open the file
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	if(exists)then
		open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
		if(ios.ne.0)then
			print*,'Problems opening ',fln(flnNum)(1:40)
			print*
			STOP ' Stopped in read_resvin_ef @ 62'
		endif
	  if(iopt.ge.1)print*,'Opened unit'
     *                      ,unitNum,' filename  ',fln(flnNum)(1:40)
	else
		print*,'ERROR: reservoir inflows (rin) file ',fln(flnNum)(1:40)
		print*,'is NOT found'
		STOP ' Program STOPPED in read_resvin_ef @ 67'
	endif

! Initialize default values
	CALL InitResvinParam(header)	

! Search for and read tb0 file header
	line(1:1) = '#'
	foundEndHeader = .false.
	insideColMetaData = .false.

	do WHILE((.NOT.foundEndHeader) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(LEN_TRIM(line) .eq. 0))) 	

		read(UNIT=unitNum, FMT='((A))', iostat=ios) line	! read a line
		if(ios .eq. -1)then
			write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
			STOP ' Stopped in read_resvin_ef'
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

			else if(keyword(1:KeyLen) .eq. ':columnmetadata')then
				insideColMetaData = .TRUE.
			else if(keyword(1:KeyLen) .eq. ':endcolumnmetadata')then
				insideColMetaData = .FALSE.
			else if(insideColMetaData) then
				iStat = ParseResvinColumnMetaData(colHeader,keyword,
     &												keyLen,subString)
				if(iStat .lt. 0) then
					write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
					write(*,'(2(A))') '   in line: ',line					
					STOP ' Stopped in read_resvin_ef'
					return
				end if
			else
				iStat = ParseResvinParam(header,keyword,keyLen,subString)
				if(iStat .lt. 0) then
					write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
					write(*,'(2(A))') '   in line: ',line					
					STOP ' Stopped in read_resvin_ef'
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

! Assign the variables from the types
	ktri =  header%tb0p%deltaT    !data time step in hours
	factor = header%tb0p%unitConv ! conversion to cms
	noresvi = colHeader%tb0cmd%colCount !no of reservoirs

d     print*
d     print*,'In read_resvin'
d     print*,'No of reservoirs found =',noresvi
d     print*


!     rev. 9.5.68  Oct.  07/09  - NK: debugged read_resvin_ef.f
      if(allocated(delta))then
        deallocate(delta,qpeakh,qpeaks,dpeakh,dpeaks,volsyn,volhyd,
     *               stat=iAllocate)
        if(iDeallocate.ne.0)then
          print*,'Error with deallocation of resvin stuff'
          print*
          stop 'Program aborted in read_resvin_ef @ 141'
        endif

        allocate(
     *    delta(nnch,no+noresvi),     ! fixed bug 24/07/00
     *    qpeakh(no+noresvi),qpeaks(no+noresvi),
     *    dpeakh(no+noresvi),dpeaks(no+noresvi),
     *    volsyn(no+noresvi),volhyd(no+noresvi),
     *    stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *       'Error with allocation error arrays in read_resvin_ef'

      endif


!     nreli    =     no of hours of data
! Scan lines of data
      rewind unitNum
	nreli = CountDataLinesAfterHeader(unitNum)
	nreli=nreli*ktri         ! want no of hours, not no of lines
	rewind unitNum
	CALL GoToStartOfData(unitNum)

!     allocate stuff here
!     array size defined by noresv which is already defined in rdresv     
	if(firstpass.eq.'y')then
		firstpass='n'
		nrel_in_max=nreli
		if(iopt.ge.1)then
		  print*,'no days of resv. inflows found =',nreli/ktri
		endif

!     rev. 9.7.12  Nov.  10/10  - NK: fix array bugs for reservoir inflows
		if(noresvi.gt.0)then
			allocate(resnamei(noresv),qinfl(noresv,nreli),
     *                        xresin(noresv),yresin(noresv),
     *                        jresin(noresv),iresin(noresv),
     *                        resin_reach(noresv),
     *                        nopti(noresv),stat=iAllocate)
			if(iAllocate.ne.0) STOP
     *			'Error with allocation in rdresvin @172'
d	        print*,'noresvi=',noresvi
		else
			allocate(resnamei(1),qinfl(1,nreli),nopti(1),
     *              stat=iAllocate)
			if(iAllocate.ne.0) STOP
     *			'Error with allocation in rdresvin @178'
		endif
	  if(iopt.ge.1)then
		  print*,'res. inflow  allocation for ',noresv,nreli
		endif
	else	
		if(nreli.gt.nrel_in_max)then
!			this happens when the current event is longer than a
!			previoius one
			nrel_in_max=nreli
!			the file fegnth is longer than any of the previous events so 
!			more memory has to be allocated
	
!			DEALLOCATION OF ARRAYS FROM AREA5A:
			deallocate(qinfl,stat=iDeallocate)     
			if (iDeallocate.ne.0) print*,    
     *			'Error with deallocation of area5a arrays'
!			re-allocate for larger arrays
!     rev. 9.8.41  Jan.  28/13  - NK: fixed bug in lst for level print statement
	        if(iopt.ge.1)then
		      print*,'re-allocating res. inflow to ',noresv,nreli
		    endif
c			allocate(qinfl(noresvi,nreli),stat=iAllocate)
			allocate(qinfl(noresv,nreli),stat=iAllocate)
			if(iAllocate.ne.0) STOP
     *			'Error with allocation of area10a arrays in read_resvin'
		endif
	endif

 
!     RESERVOIR INFLOWS
!     RESERVOIR INFLOWS
!     RESERVOIR INFLOWS

!       READ RESERVOIR INFLOWS (FOR OPTIMIZATIONS ETC.)
!       SET ALL FLAGS TO 0:
	do i=1,noresvi
		nopti(i)=0
	end do

!     NOTE:  we use unit 99 (scratch) here because unit 38 is used for 
!     snow data elsewhere

!	ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
	do i=1,noresvi
		resnamei(i) = colHeader%tb0cmd%colName(i) ! reservoir name
          xresin(i) = colHeader%tb0cmd%colLocX(i) ! x coordinate
          yresin(i) = colHeader%tb0cmd%colLocY(i) ! y coordinate
		  nopti(i) = int(colHeader%colValue1(i)) ! coefficient 1
!     rev. 9.8.41  Jan.  31/13  - NK: fixed bug in read_resvin: nopti int conversion
	end do
	deallocate(colHeader%tb0cmd%colName)
	deallocate(colHeader%colValue1)

!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
	if(noresvi.gt.noresv)then    !  nk  Jul. 28/04
		print*,'Number of lakes or reservours = ',noresv
		print*,'Number of reservoir inflows   = ',noresvi
		print*,'This can not be. no of inflows must be '
		print*,'.le. no of lakes. We will try to just ignore the'
		print*,'extra inflow.'
		noresvi=noresv
		print*
c		pause 'Hit enter to continue '
	endif

!     rev. 9.5.78  Nov.  04/09  - NK: matched resvin locations to reach numbers
!     find out what reach inflows are associated with
	do i=1,noresvi
        jresin(i)=int((xresin(i)-xorigin)/xdelta)+1
        iresin(i)=int((yresin(i)-yorigin)/ydelta)+1
	end do

      do i=1,noresv
        jres(i)=int((xres(i)-xorigin)/xdelta)+1
        ires(i)=int((yres(i)-yorigin)/ydelta)+1
        do j=1,noresvi
	      if(jres(i).eq.jresin(j).and.ires(i).eq.iresin(j))then
	        resin_reach(j)=i
            endif
	  end do
      end do

!	DATA CHECKING:
	if(iopt.ge.1.and.noresvi.gt.0)then
		write(53,5301)
		write(53,4901)(nopti(i),i=1,noresvi)
		write(53,4902)noresvi,nreli,ktri
		write(53,5303)(resnamei(i),i=1,noresvi)
	    write(53,5304)(xresin(i),i=1,noresvi)
	    write(53,5304)(yresin(i),i=1,noresvi)
 	    write(53,5305)(jresin(i),i=1,noresvi)
 	    write(53,5305)(iresin(i),i=1,noresvi)
	    write(53,5305)(resin_reach(i),i=1,noresvi)
d		write(*,5301)
d		write(*,4901)(nopti(i),i=1,noresvi)
d		write(*,4902)noresvi,nreli,ktri
d		write(*,5303)(resnamei(i),i=1,noresvi)
d	    write(*,5304)(xresin(i),i=1,noresvi)
d	    write(*,5304)(yresin(i),i=1,noresvi)
d	    write(*,5305)(jresin(i),i=1,noresvi)
d	    write(*,5305)(iresin(i),i=1,noresvi)
d	    write(*,5305)(resin_reach(i),i=1,noresvi)
	endif

!     rev. 9.5.82  Jan.  26/09  - NK: replaced error check for inflow locations
      errflg1=.false.
	do i=1,noresvi
	  if(resin_reach(i).le.0.or.resin_reach(i).gt.noresv)then
	    print*,'No matching lake/reservoir for inflow '
	    print*,'location ',i,resnamei(i)
	    print*,'Please check the location of the inflow'
	    print*,'in the yyyymmdd_rin.tb0 file'
	    print*,'and in inflow\inflow_locations.xyz'
	    print*
	    errflg1=.true.
	  endif
	end do
	if(errflg1)stop 'Program aborted in read_resvin_ef.f @ 288'

	if(noresvi.gt.0)then
		do j=ktri,nreli,ktri
			read(unitNum,*)(qinfl(k,j),k=1,noresvi)
			write(53,4904)(qinfl(k,j),k=1,noresvi)
cd			write(*,4904)(qinfl(k,j),k=1,noresvi)
		end do
	endif
	close(unit=unitNum,status='keep')

 4901 format(<noresvi>i1)
 4902 format(3i5)
 4904 format(256f10.3)
 5301 format(' ','Reservoir inflow data echoed:')
 5303 format(256(' ',a12))
 5304 format(256(f13.4))
 5305 format(256(i13))
	end subroutine read_resvin_ef