      SUBROUTINE read_swe()

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
     
!***************************************************************************
! PROGRAM BY: NK  June 2012
!
      use area_watflood
	USE EF_Module

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
	implicit none
	save

	TYPE(CourseParam) :: header
	TYPE(CourseColumnMetaData) :: colHeader

      CHARACTER(1) :: new
      INTEGER      :: jan,ju,iallcnt2,iallocate,n,ii,i,j,nold
      integer      :: n_max,ideallocate,no_hdr_lines
      integer      :: ios,linecount,unitnum,flnNum,iStat
      integer      :: nlines,nlines_old
	REAL*4       :: time,ttime,taold,temmmp  
	character*16 :: swe_file              
	logical      :: exists,firstpass

! Local variables
	character*128 keyword, value
	character*4096 line, subString, tmpString
      CHARACTER(1)  :: col0(13),col1(13),col2(27)
	integer lineLen, keyLen, wordCount
	logical rStat, lineType, foundEndHeader, insideColMetaData

      save
	
      DATA new/'f'/
      DATA firstpass/.true./
	DATA iallcnt2/0/

      DATA col0/'a','c','e','g','i','k','m','o','q','s','u','w','y'/
      DATA col1/'b','d','f','h','j','l','n','p','r','t','v','x','z'/
      DATA col2/' ','a','b','c','d','e','f','g','h','i','j','k','l','m',
     *         'n','o','p','q','r','s','t','v','u','w','x','y','z'/

d      if(iopt.eq.2)print*,'in melt @ 0'

      stop

      if(iopt.ge.1)then
        print*,'~~~~~~~~~~~~~~~~~~~~~~~~'
        print*,'NEW  <<<<<'
        print*,'SWE report'
      endif
! Set unit and fln number
	unitNum = 284
	flnNum = 54      ! ! unit=284  fln(54)- swe time series pillows & crs  yyyymmdd_swe.tb0

      write(line,10000)year_now
      if(iopt.ge.1)print*,'year_now=',year_now
10000 format('snowg\',i4,'0101_swe.tb0')
      read(line,10001)fln(flnNum)
10001 format(a22)            
      if(iopt.ge.1)then
        print*,'fl#=',flnNum
        print*,'file name =',fln(flnNum)(1:50)
      endif
     
c      pause 1
! Open the file
!     rev. 10.1.81 May   05/17  - NK: Added snowg\yyyymmdd_swe.tb0 obs. swe
	INQUIRE(FILE=fln(flnNum),EXIST=exists)
	if(exists)then
	  open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
	    if(ios.ne.0)then
	  	print*,'Problems opening:',fln(flnNum)(1:40)
	   	print*
		   STOP ' Stopped in read_swe @ 39'
	  endif
	  print*,'Opened ',fln(flnNum)(1:40)
	else
        if(dds_flag.eq.0)then
          print*,fln(flnNum)(1:60),'not found'
	    print*,'program continues without swe analysis'
	  endif
	  courseflg=.false.
	  return
	endif

c      pause 2

! Initialize default values
	CALL InitCourseParam(header)	
	
c	pause 3
	
! Search for and read tb0 file header
      linecount=0
      line(1:1) = '#'
      foundEndHeader = .false.
      insideColMetaData = .false.
      no_hdr_lines=0

      write(51,*)
      write(51,*)'SWE - observed snow course/pillow file'
      write(51,*)

      do WHILE((.NOT.foundEndHeader) .AND.
     &        ((line(1:1) .eq. '#') .OR.
     &        (line(1:1) .eq. ':') .OR.
     &        (LEN_TRIM(line) .eq. 0)))   
                  linecount=linecount+1
!     if(iopt.eq.2)print*,'reading line ',linecount,' in read_swe'
          read(UNIT=unitNum, FMT='((A))', iostat=ios) line    ! read a line
          write(51,*)line(1:60)
          if(ios .eq. -1)then
              write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
              STOP ' Stopped in read_swe @ 88'
          end if
          no_hdr_lines=no_hdr_lines+1
          
          rStat = Detab(line)             ! replace tabs with spaces
          line = ADJUSTL(line)        ! Get rid of leading white space
          lineLen = LEN_TRIM(line)        ! Find the length excluding trailing spaces

          if(line(1:1) .eq. ':')then
              wordCount = SplitLine(line, keyword, subString) ! find the keyword
              rStat = ToLowerCase(keyword)
              KeyLen = LEN_TRIM(keyword)

              if(keyword(1:KeyLen) .eq. ':endheader')then
                  foundEndHeader = .TRUE.
              else if(keyword(1:KeyLen) .eq. ':columnmetadata')then
                  insideColMetaData = .TRUE.
              else if(keyword(1:KeyLen) .eq. ':endcolumnmetadata')then
                  insideColMetaData = .FALSE.
              else if(insideColMetaData) then
                  iStat = ParseCourseColumnMetaData(colHeader,
     &                        keyword,keyLen,subString)
                  if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_swe @113'
                      return
                  endif
              else
                  iStat = ParseCourseParam
     *                (header,keyword,keyLen,subString)
                  if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                      write(*,'(2(A))') '   in line: ',line                   
                      STOP ' Stopped in read_swe @122'
                      return
                  else if(iStat .eq. 0) then
!                     write(*,'((A), (A))')  'Unrecognized keyword line: ',
!     &                                       line
                  endif
              end if
          endif
          
      end do  
!***************************************
!     Finished reading the header
!***************************************

!     nswe    = no of snow course data
	nswe    = colHeader%tb0cmd%colCount
      print*,'no of snow courses found =',nswe
      print*,'no header lines =',no_hdr_lines
      
! Scan lines of data
      rewind unitNum
	nlines = CountDataLinesAfterHeader(unitNum)
	print*,'no of data lines found =',nlines
      print*,'~~~~~~~~~~~~~~~~~~~~~~~~'
	
      rewind unitnum
      do n=1,no_hdr_lines
        read(unitnum,*)line
c        print*,line(1:60)
      end do
      

!       rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!       allocate stuff      
c		if(id.eq.1)then
		if(firstpass)then
			n_max=n
			nold=n
			nlines_old=nlines
	      allocate(course_obs(nswe,max(366,nlines)),
     *              course_calc(nswe,max(366,nlines)),
     *             	      indomainflg(nswe),stat=iAllocate)
!           corrected bug above - replaced na by nswe   Nov. 19/14  nk     
            If(iopt.ge.1)print*,'course_obs dimensioned for:',
     *                        max(nlines,366)
			if(iAllocate.ne.0) STOP
     *		'Error with allocation of  arrays in read_swe @ 147'
		else                    !  firstpass
!         check to see memory allocated is adequate      

			if(n.ne.nold)then
				print*,'No of swe stations has been changed in'
				print*,'nold=',nold,' n=',n
				print*,'in file ',fln(99)(1:60)
				print*,'This is not permitted'
				print*
				stop 'Program aborted in read_swe @ 156'
			endif
			if(nlines.gt.nlines_old)then
				print*,'No of data lines have changed in'
				print*,'nlines_old=',nlines_old,' nlines=',nlines
				print*,'in file ',fln(99)(1:60)
				print*
c			if(n.gt.n_max)then
				n_max=n

!           the file length is longer than any of the previous events so 
!           more memory has to be allocated

!           DEALLOCATION OF ARRAYS FROM AREA10A:
				deallocate(course_obs,course_calc,stat=iDeallocate)
				if(iDeallocate.ne.0)then
					print*,'Error with deallocation ofswe)obs'
					print*
					stop 'Program aborted in read_swe @ 170'
				endif

                print*,'reallocation for more swe locations',nswe

				allocate(course_obs(nswe,max(366,nlines)),
     *                    course_calc(nswe,max(366,nlines)),
     *         				stat=iAllocate)
				if(iAllocate.ne.0) STOP
     *			'Allocation Error: arrays in read_swe @177'
                nlines_old=nlines
d               print*,'course_obs re-dimensioned for:',max(nlines,366)
			endif
		endif                   !  firstpass
		
      if(.not.allocated(xswe))then
          allocate(xswe(nswe),yswe(nswe),gname_swe(nswe),
     *                   ewg_swe(nswe),sng_swe(nswe),stat=iAllocate)
          if(iAllocate.ne.0)then
            print*,'Error with allocation of station descripters'
            print*,' in read_swe'
            STOP 'Program aborted in read_swe @ 213'
          endif
        endif
!     rev. 9.1.68  Dec.  19/04  - NK: rewrote read_tbo c/w memory allocation 
  
!     ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
      do n=1,nswe
        gname_swe(n) = colHeader%tb0cmd%colName(n) ! station name
        xswe(n) = colHeader%tb0cmd%colLocX(n) ! x coordinate
        yswe(n) = colHeader%tb0cmd%colLocY(n) ! y coordinate
      end do
      deallocate(colHeader%tb0cmd%colName)
      deallocate(colHeader%tb0cmd%colLocX)
      deallocate(colHeader%tb0cmd%colLocY)

c      pause '261'

!     write the header on the swe.csv file
c      filename(951)='results\swe.csv'       !snow course comparison

c      if(id.eq.1)then
      if(firstpass)then
        write(951,91101)(gname_swe(n),gname_swe(n),n=1,nswe)
91101   format(<2*nswe>(a9,','))      
      endif

c      pause '269'

c      if(id.eq.1)then
      if(firstpass)then
!       turn into local coordinates
        do n=1,nswe
          ewg_swe(n)=int((xswe(n)-xorigin)/xdelta+1.0)
          sng_swe(n)=int((yswe(n)-yorigin)/ydelta+1.0)
        end do

!       write snow course & pillow the location file
        open(unit=99,file='swe_location.xyz',
     *            status='unknown',iostat=ios)
        i=0
        j=-3
        do n=1,nswe
          i=n-n/13*13
          if(i.eq.0)then
            i=13
          endif
          j=(n-1)/13+1
          if(sng_swe(n).ge.1.and.ewg_swe(n).ge.1.and.
     *       sng_swe(n).le.ycount.and.ewg_swe(n).le.xcount)then
            if(s(sng_swe(n),ewg_swe(n)).gt.0)then
!             print only if in the watershed
              write(99,6012)xswe(n),yswe(n),n,gname_swe(n),
     *                col2(j),col0(i),col2(j),col1(i)
 6012         format(2f12.3,i5,1x,a12,3x,2a1,3x,2a1)
            endif
          endif
        end do
        close(unit=99,status='keep')
      endif

c      pause 297


      do j=1,nlines   
        read(unitNum,*,iostat=ios)(course_obs(n,j),n=1,nswe)
c        print*,j,(course_obs(n,j),n=1,4)
        if(ios.ne.0)then
          print*, 'In read_swe.f'
          write(*,*)j,(course_obs(n,j),n=1,4)
          print*,'Got as far as data line',j
        endif
      end do
      
      close(unit=99,status='keep')

!     check to see if stations are in the model domain
      do n=1,nswe
        if(ewg_swe(n).ge.1.and.ewg_swe(n).le.xcount.and.
     *     sng_swe(n).ge.1.and.sng_swe(n).le.ycount)then
          if(s(sng_swe(n),ewg_swe(n)).gt.0)then
            indomainflg(n)=.true.
d           print*,n,xswe(n),yswe(n),s(sng_swe(n),ewg_swe(n)),
d    *         indomainflg(n)
          else
            indomainflg(n)=.false.
          endif
        else
          indomainflg(n)=.false.
        endif
      end do
c      pause 'domain check'

      close(unit=unitNum,status='keep')

      firstpass=.true.  ! reset for next event to read header only
      firstpass=.false.  ! reset for next event to read header only

      return

      END SUBROUTINE read_swe
