      SUBROUTINE read_temp_ef(hdrflg,jan,jz)

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
     
!***********************************************************************
!  READ_TEMP_EF - written May/06 by Dave Watson
!	- Derived from rdtemp written by Nick Kouwen
!	- This subroutine reads the ensim compatible gridded temperature
!	  (TEM) file (r2c format)
!***********************************************************************
!     rev. 9.7.04  Aug.  30/10  - NK: added to error message in read_rain & read_temp
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
!     rev. 9.7.18  Jan.  17/11  - NK: Changed tolerance on the grid check in rear_rain & read_temp

      use area_watflood
      use EF_Module

      implicit none
      type(TempParam) :: header
      type(FrameRecord) :: frameRec

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      DIMENSION :: ntoc(9000)

      real*4  :: x1,deltatmp,teltatmp
      integer :: idlast,mhlast,jan,ntocorrect,ios,ios1,
     *           ntoc,i,iend,j,jz,n,
     *           ndeltatt,tem_hr,iostat,iAllocate,iDeallocate,
     *           nh_count,
     *           nolines,xcount_max,ycount_max
!      character(10) ::  fileformat,starttime,startdate
      character(10) ::  fileformat
      real*4	timeHrs

      character*20  :: junk
      character*1   :: junk2,newfmtflg,firstpass,lineflg,hdrflg
      character*6   :: junk0
      DATA idlast/0/
      data firstpass/'y'/
      LOGICAL      :: exists

     
! Local variables
      integer*4 unitNum, flnNum, iStat
      character*4096 line, subString
      character*128 keyword, value
      integer lineLen, keyLen, wordCount
      logical rStat, foundEndHeader,foundFirstFrame,foundFrame
      integer frameCount

      integer modelHour !used to be k
      integer modelHourLast !used to be klast
      integer dataHour !used to be nr_count

! Dave's debug variables
!	integer daveInt
!	real daveReal

! Initialize default values within frame module
      CALL InitFrameRecord(frameRec)

! Set unit and fln number
      unitNum = 45
      flnNum = 15

      foundEndHeader = .false.

! Debug line
!	print*, 'Inside read_temp_ef at time = ',timeHrs, 'hrs'

! If hdrflg==1 ,then it's the first time through for this event
      if(hdrflg.eq.'1')then
c        open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
c        if(ios.ne.0)then
c            print*, 'Problems opening the tem file',fln(flnNum)
c            STOP ' Stopped in read_temp_ef'
c        endif


		inquire(FILE=fln(flnNum),EXIST=exists)
		if(exists) then
		  open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
	    if(ios.ne.0)then
		    print*, 'Problems opening file',fln(flnNum)(1:40)
		    STOP ' Stopped in read_temp_ef @ 85'
		  endif
	    if(iopt.ge.1)print*,'Opened unit'
     *                      ,unitNum,' filename  ',fln(flnNum)(1:40)
		else
		  print*, 'Attempting to open file name ',fln(flnNum)(1:40)
	    print*, 'Unit number ',unitNum
		  print*, 'but it is not found (in this location)'
		  STOP 'Program aborted in read_temp_ef @95'
		endif



! Initialize default values within rain module
        call InitTempParam(header)	

! Search for and read r2c file header
        line(1:1) = '#'
        do while((.NOT.foundEndHeader) .AND.
     &	    ((line(1:1) .eq. '#') .OR.
     &		(line(1:1) .eq. ':') .OR.
     &		(len_trim(line) .eq. 0))) 	

            read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line
            if(ios .eq. -1)then
                print*, 'ERROR: Premature EndOfFile encountered'
                STOP ' Stopped in read_temp_ef'
            end if

            rStat = Detab(line)				! replace tabs with spaces
            line = adjustl(line)		! Get rid of leading white space
            lineLen = len_trim(line)		! Find the length excluding trailing spaces

            if(line(1:1) .eq. ':')then
                wordCount = SplitLine(line, keyword, subString)	! find the keyword
                rStat = ToLowerCase(keyword)
                KeyLen = len_trim(keyword)

                if(keyword(1:KeyLen) .eq. ':endheader')then
                    foundEndHeader = .TRUE.
                else
                    ! parse the header
                    iStat = ParseTempParam(header,keyword,keyLen,
     &													subString)
                    if(iStat .lt. 0) then
                       write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                        write(*,'(2(A))') '   in line: ',line					
                        STOP ' Stopped in read_temp_ef'
                        return
                    else if(iStat .eq. 0) then
!						write(*,'((A), (A))')
!     &								'Unrecognized keyword line: ',line
                    endif
                end if
            end if
        end do

! Assign the parsed parameters to the model variables		
        xcount3 = header%r2cp%xCount
        ycount3 = header%r2cp%yCount
        xorigin3 = header%r2cp%xOrigin
        yorigin3 = header%r2cp%yOrigin
        xdelta3 = header%r2cp%xDelta
        ydelta3 = header%r2cp%yDelta
        
c!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
c        if(xcount3.ne.xcount_precip.or.ycount3.ne.ycount_precip)then
c          print*
c          print*,'Fatal error:'
c          print*,'xcount3 =',xcount3,'   xcount_precip =',xcount_precip
c          print*,'ycount3 =',ycount3,'   ycount_precip ='ycount_precip
c          print*,'Precip and temperature grids do not match'
c          stop 'Program aborted in read_temp_ef @ 141'
c        endif

! Scan data frames
        frameCount = 0
        do while((.NOT.EOF(unitNum)))
            read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line
            if(ios .eq. -1)then
             write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
                STOP ' Stopped in read_temp_ef'
            end if

            rStat = Detab(line)				! replace tabs with spaces
            line = adjustl(line)		! Get rid of leading white space
            lineLen = len_trim(line)		! Find the length excluding trailing spaces

            if(line(1:1) .eq. ':')then
                wordCount = SplitLine(line, keyword, subString)	! find the keyword
                rStat = ToLowerCase(keyword)
                KeyLen = len_trim(keyword)
                
                if(keyword(1:KeyLen).eq.':frame')then
                    iStat = ParseFrameLine(frameRec,keyword,keyLen,
     &													subString)
                    
!	Identify the first and last frame's timestamp 
                    if(iStat .lt. 0) then
                       write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                        write(*,'(2(A))') '   in line: ',line					
                        STOP ' Stopped in read_temp_ef'
                        return
                    else if(iStat .eq. 0) then
!						write(*,'((A), (A))')
!     &								'Unrecognized keyword line: ', line
                    else if(frameRec%frame.EQ.1) then
                        header%startJulianDay =
     &						JDATE(frameRec%tStamp%year,
     &							frameRec%tStamp%month,
     &							frameRec%tStamp%day)
                        header%startHour = frameRec%tStamp%hour
                    else
                        header%endJulianDay =
     &						JDATE(frameRec%tStamp%year,
     &							frameRec%tStamp%month,
     &							frameRec%tStamp%day)
                        header%endHour = frameRec%tStamp%hour
                    end if
                    header%r2cp%frameCount = header%r2cp%frameCount+1
                end if
            end if
        enddo	
      
!	deltat2 = model timestep in hours 
        deltat2 = 1	
!	njtemp = number of hours spanned by this tem file
!	convert to julian days to properly calculate hours spanned
        nhtemp = (header%endJulianDay - header%startJulianDay)*24
     &				 + (header%endHour - header%startHour) + 1

!		nr = (header%endJulianDay - header%startJulianDay)*24
!     &				 + (header%endHour - header%startHour) + 1

!     Dave used nr instead of nhtemp here becasue file was copied from read_rain_ef
!     we need a different variable for the # hours of data for the tem file
!     nk  Nov. 23/06

        
!	Position to start of data (immediately after first frame record but before first frame data)
        REWIND (unitNum)
        foundFirstFrame = .false.
        do WHILE(.NOT.foundFirstFrame)
            read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
            temp_frame_header=line
            new_temp_frame_header=line
            temp_flg=.true.
            if(line(1:1) .eq. ':')then
                wordCount = SplitLine(line, keyword, subString)	! find the keyword
                rStat = ToLowerCase(keyword)
                KeyLen = len_trim(keyword)
                if(keyword(1:KeyLen).eq.':frame')then
                    iStat = ParseFrameLine(frameRec,keyword,
     &												keyLen,subString)
                    foundFirstFrame = .true.
                end if
            end if		      
        enddo	
        
        modelHour=1 ! set model hour to first hour
        dataHour = 0 ! set data hour to just before first hour

! Are we LatLong or Projected
        if(IsLatLong(header%r2cp%csp).eq.(.true.))then
            rgrde=xdelta2*60.
            rgrdn=ydelta2*60.  
            rads=int(yorigin2*60.0)
            radn=int((yorigin2+ycount2*ydelta2)*60.0) 
            radw=int(xorigin2*60.0)  
            rade=int((xorigin2+xcount2*xdelta2)*60.0)  
        else
            rgrde=xdelta2/1000.
            rgrdn=ydelta2/1000. 
            rads=int(yorigin2/1000.)
            radn=iymin+grde*(ycount2-1)
            radw=int(xorigin2/1000.)
            rade=jxmin+grdn*(xcount2-1)
            if(abs(xdelta2-ydelta2).lt.0.001)then  ! changed jan20/08 nk
                rgrd=xdelta2/1000.
            else
              print*, 'ERROR - xdelta .ne. ydelta'
              print*,' xdelta =',xdelta
              print*,' ydelta =',ydelta
              STOP 'Program aborted in read_temp_ef @ 249'
            endif
        endif

! firstpass means this is the first tem file of the first event...which is the only time it needs to be done
        if(firstpass.eq.'y')then
            firstpass='n'      
            xcount_max=max(xcount+1,xcount3)
            ycount_max=max(ycount+1,ycount3)
            allocate(ttemp(ycount_max,xcount_max),tempv(na),
     *						tempvmin(na),rh(na),stat=iAllocate)
            if(iAllocate.ne.0) STOP
     *	   'Error with allocation in read_temp_ef'
c			xcount_max=xcount3
c			ycount_max=ycount3
        endif

 
! Nick initializes stuff here...may want to move this to process_temp()
        do n=1,na
            tempvmin(n)=99.99
            rh(n)=.50
        enddo

!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability

        if(jan.eq.1.and..not.netCDFflg)then     
!       CALCULATE THE GRID OFFSETS:
            if(iymin.ge.rads.and.jxmin.ge.radw)then
!         COMPARE THE SOUTH-WEST CORNER COORDINATES
                if(IsLatLong(header%r2cp%csp).eq.(.true.))then
                    iyoffset=int(float(iymin-rads)/rgrdn)+iyshift
                    jxoffset=int(float(jxmin-radw)/rgrde)+jxshift
                else
                    iyoffset=int(float(iymin-rads)/grdn)+iyshift
                    jxoffset=int(float(jxmin-radw)/grde)+jxshift
                endif
            elseif(abs(iymin-rads).gt.1.or.abs(jxmin-radw).gt.1)then
!     rev. 9.7.18  Jan.  17/11  - NK: Changed tolerance on the grid check in rear_rain & read_temp
!         WATERSHED GRID does not coincide with met gris
                  print*
                  print*,'WARNING:'
                  print*,' In read_temp_ef:'
                  print*,'tmp grid does not coincide with shd file'
                  print*,'This could be due to round off when using'
                  print*,'fractions of degrees as the grid size'
                  print*,'iymin,jxmin,rads,radw/',iymin,jxmin,rads,radw
                  print*,'Please check the headers on the shd & tmp'
                  print*,'for compatibility'
                  print*                     
                  STOP 'Program aborted in read_tmp_ef @ 298'
            else
                iyoffset=0
                jxoffset=0
                ! DO NOTHING
            endif
c           what is this doing here??????????????????????????????
c			call precip_adjust()
        endif   ! END OF JAN=1
        
! Check for change of grid size
        if(xcount3.gt.xcount_max.or.ycount3.gt.ycount_max)then
            deallocate(ttemp,stat=iDeallocate)
            if(iDeallocate.ne.0) then
                print*,	'Warning: error with deallocation of outarray'
            endif

            xcount_max=xcount2
            ycount_max=ycount2
!       outarray is in areawfo
!       this has to be allocated before calling write_r2c
            allocate(ttemp(ycount_max,xcount_max),stat=iAllocate)
            if(iAllocate.ne.0) then
                STOP 'Error with allocation of outarray in tr_sub'      
            end if
        endif


        return
      endif          ! (hdrflg.eq.'1')

!***********************************************************
!***********************************************************
      
!     Finished reading the header      

!***********************************************************
!***********************************************************

! Moved from above        
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
        if(xcount3.ne.xcount_precip.or.ycount3.ne.ycount_precip)then
          print*
          print*,'Fatal error:'
          print*,'xcount3 =',xcount3,'   xcount_precip =',xcount_precip
          print*,'ycount3 =',ycount3,'   ycount_precip =',ycount_precip
          print*,'Precip and temperature grids do not match'
          stop 'Program aborted in read_temp_ef @ 141'
        endif


!     We are looking for temp data for this step
      dataHour = dataHour + 1

      if(modelHour.eq.modelHourLast.and.dds_flag.eq.0)then
c        print*,'Hr=',modelhour,' Error in ',
c     *                 fln(flnNum)(1:40),' met data repeated??'
         write(*,*)'Error in the tem file @ ',frameRec%frame
c     &            frameRec%tStamp%year,'/',frameRec%tStamp%month,'/',
c     &            frameRec%tStamp%day,'/',frameRec%tStamp%hour
         write(51,*)'Error in the tem file @ frame #',frameRec%frame
c     &            frameRec%tStamp%year,'/',frameRec%tStamp%month,'/',
c     &            frameRec%tStamp%day,'/',frameRec%tStamp%hour
51000     format(a24,i4,a1,i2,a1,i2,a1,i2)    
      endif    
 
!	If there is a record for this hour than modelHour = dataHour
!	print*,	'modelHour:',modelHour, 'dataHour:',dataHour
!	if(modelHour.eq.dataHour)then

!     This is changed so if we've passed the time this data should have
!     been read - we'd better read it now      
      if(modelHour.eq.jz)then
          
!	Go ahead and read the data for this frame 
        do i=1, ycount3
!           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            read(45,*,iostat=ios)(ttemp(i,j),j=1,xcount3)
!           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if(ios.ne.0)then
                write(*,9993)modelHour,i
                print*,'iostat =',ios
                print*,' last data read:'
                print*,' jradw= 1 jrade= ',jrade
                print*,' iradn= ',iradn,' irads= 1'
                print*,'Event #',id
                PRINT*,'DataHour,ModelHour(jz)',modelHour,jz
                STOP ' program aborted - read_temp_ef @ 396'
            endif
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
            temp_frame_header=new_temp_frame_header
            temp_flg=.true.
            
        end do
        
!     rev. 10.2.48 Feb.  25/19  - NK: Moved temp correction to read_temp from process_temp
!     RAISE OR LOWER THE TEMPERATURE FILED:
      if(scaletem.gt.0.00001.or.scaletem.lt.-0.00001)then
         do i=1,ycount3
            do j=1,xcount3
               ttemp(i,j)=ttemp(i,j)+scaletem
            end do
         end do
      endif
      

!	  Read the next frame. If we are at the end of the file...close
        foundFrame = .false.
        do WHILE(.NOT.foundFrame)
            read(unit=unitNum, FMT='((A))', iostat=ios) line	! read a line
            if(EOF(unitNum)) then 
                close(unit=unitNum,status='keep',iostat=ios)
                RETURN
            endif

            if(line(1:1) .eq. ':')then
                wordCount = SplitLine(line, keyword, subString)	! find the keyword
                rStat = ToLowerCase(keyword)
                KeyLen = len_trim(keyword)
                if(keyword(1:KeyLen).eq.':frame')then
                    iStat = ParseFrameLine(frameRec,keyword,
     &												keyLen,subString)
                    foundFrame = .true.
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
                    new_temp_frame_header=line
                    new_temp_flg=.true.
                end if
            end if		      
        enddo	

        modelHourLast=modelHour

!	Determine the next modelhour		
        modelHour =	(JDATE(frameRec%tStamp%year,frameRec%tStamp%month,
     &				frameRec%tStamp%day) - header%startJulianDay) * 24
     &				+ (frameRec%tStamp%hour - header%startHour) + 1

!       check added May 2/12  nk
c        write(555,55500)modelHour,
c     &                frameRec%tStamp%year,frameRec%tStamp%month,
c     &				frameRec%tStamp%day,frameRec%tStamp%hour,
c     &                header%startJulianDay,header%startHour,
c     &                JDATE(frameRec%tStamp%year,frameRec%tStamp%month,
c     &				frameRec%tStamp%day) 
c55500   format(10i10)
!     rev. 10.2.24 May   21/18  - NK: Added error message in Read_rain & read_tmp
        if(modelHour-modelHourLast.lt.0)then
          print*
          print*,'Error:'
          print*,'file name ',fln(flnNum)(1:50)
          print*,'modelHour =',modelHour
          print*,'modelHourLast =',modelHourLast
          print*,'Note: Frame numbers must be monotomically increasing'
          print*,'in an r2c time series file.'
          print*,'year ',frameRec%tStamp%year
          print*,'month',frameRec%tStamp%month
          print*,'day  ',frameRec%tStamp%day
          print*,'hour ',frameRec%tStamp%hour
          print*
          print*,'possible cause: month repeated in the r2c file'
          print*
          stop 'Program aborted in read_temp_ef @ 411'
        endif

!		do nothing...keep last temperature values...

      endif

      RETURN

!     FORMATS

 9993 format(' error reading temp at hour/line ',2i10/)
      

      END SUBROUTINE read_temp_ef

