c      SUBROUTINE read_rain_ef(hdrflg,jz,jan,timeHrs)
      SUBROUTINE read_rain_ef(hdrflg,conv,jz,jan)

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
     
c      arguments changed  Nov. 9/06 nk
c      SUBROUTINE read_rain_ef(hdrflg,conv,scale,jan,timeHrs)
!     rev. 9.5.30  May.  26/08  - NK: conv back in read_rain & process_rain arg. list
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced Precip & temp files for sub-basins
!     rev. 9.7.18  Jan.  17/11  - NK: Changed tolerance on the grid check in rear_rain & read_temp
!     rev. 9.8.59  May   14/13  - NK: REmoved psmear & punused from the program
!     rev. 9.9.07  Jan.  10/14  - NK: Overhaul of the frame numbers to EnSim specs
 
C*****************************************************************************
C  READ_RAIN_EF - written Mar/06 by Dave Watson, CHC
C     - Derived from rdrain written by Nick Kouwen
C     This subroutine reads the ensim compatible gridded rain (MET) file 
C     (r2c format)
C*****************************************************************************

c     use area1
c     use area2
c      use area3
c     use area6   ! added for sum_precip  nk Dec. 29/06 nk
c     use area10
c     use area12
c     use area16

      use area_watflood

C R2C data module
      use EF_Module
      implicit none
      type(RainParam) :: header
      type(FrameRecord) :: frameRec


C SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      save

      character(1) :: hdrflg,firstpass
      logical      :: exists
      integer  :: jan,i,j,n,ios,jz,
     *              iAllocate,iDeallocate,
     *              xcount_max,ycount_max
      real*4   :: conv,scale
      real*4      timeHrs
      data firstpass/'y'/

C Local variables
      integer*4 unitNum, flnNum, iStat
      character*4096 line, subString
      character*128 keyword, value 
      character*10 junk
      integer lineLen, keyLen, wordCount
      logical rStat, foundEndHeader,foundFirstFrame,foundFrame
      integer frameCount
      integer framesSkipped !# of frames skipped beforefirst data /nk
      integer modelHour !used to be k
      integer modelHourLast !used to be klast
      integer dataHour !used to be nr_count
      integer yearLLL
      integer monthLLL 
      integer dayLLL 
      integer hourLLL

C Dave's debug variables
C     integer daveInt
C     real daveReal

C Initialize default values within frame module
      CALL InitFrameRecord(frameRec)

C Set unit and fln number
      unitNum = 40
      flnNum = 10
      
C Debug line
C     print*, 'Inside read_rain_ef at time = ',timeHrs, 'hrs'

C If hdrflg==1 ,then it's the first time through for this event
      if(hdrflg.eq.'1')then
       
c            open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
c            if(ios.ne.0)then
c              print*, 'Problems opening the met file',fln(flnNum)
c              STOP ' Stopped in read_rain_ef @ 82'
c            endif
            
	  	inquire(FILE=fln(flnNum),EXIST=exists)
		    if(exists) then
		      open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
	          if(ios.ne.0)then
		        print*, 'Problems opening file',fln(flnNum)(1:40)
		        STOP ' Stopped in read_rain_ef @ 89'
		      endif
	          if(iopt.ge.1)print*,'Opened unit'
     *                      ,unitNum,' filename  ',fln(flnNum)(1:40)
                if(iopt.ge.2)then
                  print*,'If program dies here, check to see'
                  print*,'there`s data in the met file'
                endif
		    else
		      print*, 'Attempting to open file name ',fln(flnNum)(1:40)
	          print*, 'Unit number ',unitNum
		      print*, 'but it is not found (in this location)'
		      STOP 'Program aborted in read_rain @95'
		    endif

C Initialize default values within rain module
            call InitRainParam(header)    

C Search for and read r2c file header
            line(1:1) = '#'
            foundEndHeader = .false.
            do while((.NOT.foundEndHeader) .AND.
     &          ((line(1:1) .eq. '#') .OR.
     &            (line(1:1) .eq. ':') .OR.
     &            (len_trim(line) .eq. 0)))     

                  read(unit=unitNum, FMT='((A))', iostat=ios) line      ! read a line
                  if(ios .eq. -1)then
                        print*, 'ERROR: Premature EndOfFile encountered'
                        STOP ' Stopped in read_rain_ef'
                  end if

d      if(iopt.eq.3)print*,line(1:40)

                  rStat = Detab(line)                       ! replace tabs with spaces
                  line = adjustl(line)          ! Get rid of leading white space
                  lineLen = len_trim(line)            ! Find the length excluding trailing spaces

                  if(line(1:1) .eq. ':')then
                        wordCount = SplitLine(line, keyword, subString) ! find the keyword
                        rStat = ToLowerCase(keyword)
                        KeyLen = len_trim(keyword)

                        if(keyword(1:KeyLen) .eq. ':endheader')then
                              foundEndHeader = .TRUE.
                        else
                              ! parse the header
                          iStat = ParseRainParam(header,keyword,keyLen,
     &                                                      subString)
                              if(iStat .lt. 0) then
                        write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                                  write(*,'(2(A))') '   in line: ',line                             
                                    STOP ' Stopped in read_rain_ef'
                                    return
                              else if(iStat .eq. 0) then
C                                   write(*,'((A), (A))')
C     &                                               'Unrecognized keyword line: ',line
                              endif
                        end if
                  end if
            end do

C Assign the parsed parameters to the model variables       
            xcount2 = header%r2cp%xCount
            ycount2 = header%r2cp%yCount
            xorigin2 = header%r2cp%xOrigin
            yorigin2 = header%r2cp%yOrigin
            xdelta2 = header%r2cp%xDelta
            ydelta2 = header%r2cp%yDelta
            
            conv = header%r2cp%unitConv   
            
!           for grid consistency checking
            xcount_precip=xcount2
            ycount_precip=ycount2            
            
C Scan data frames
            frameCount = 0
            do while((.NOT.EOF(unitNum)))
                  read(unit=unitNum, FMT='((A))', iostat=ios) line      ! read a line
d                 if(iopt.eq.3)print*,line(1:40)
                  if(ios .eq. -1)then
               write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
                        STOP ' Stopped in read_rain_ef'
                  end if

                  rStat = Detab(line)                       ! replace tabs with spaces
                  line = adjustl(line)          ! Get rid of leading white space
                  lineLen = len_trim(line)            ! Find the length excluding trailing spaces

                  if(line(1:1) .eq. ':')then
                        wordCount = SplitLine(line, keyword, subString) ! find the keyword
                        rStat = ToLowerCase(keyword)
                        KeyLen = len_trim(keyword)
                        
                        if(keyword(1:KeyLen).eq.':frame')then
                         iStat = ParseFrameLine(frameRec,keyword,keyLen,
     &                                subString)
                              frameCount = frameCount+1
                              
C     Identify the first and last frame's timestamp 
                              if(iStat .lt. 0) then
                         write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                                 write(*,'(2(A))') '   in line: ',line                             
                                    STOP ' Stopped in read_flow_ef'
                                    return
                              else if(iStat .eq. 0) then
C                                   write(*,'((A), (A))')
C     &                                               'Unrecognized keyword line: ', line
                              else if(frameCount.EQ.1) then
                                    header%startJulianDay =
     &                                    JDATE(frameRec%tStamp%year,
     &                                          frameRec%tStamp%month,
     &                                          frameRec%tStamp%day)
                                 header%startHour = frameRec%tStamp%hour
                              else
                                    header%endJulianDay =
     &                                    JDATE(frameRec%tStamp%year,
     &                                          frameRec%tStamp%month,
     &                                          frameRec%tStamp%day)
                                   header%endHour = frameRec%tStamp%hour
                              end if
                              header%r2cp%frameCount = frameCount
                        end if
                  end if
            enddo 
      
C     deltat2 = model timestep in hours 
            deltat2 = 1 
C     nr = number of hours spanned by this met file
C     convert to julian days to properly calculate hours spanned
            nr = (header%endJulianDay - header%startJulianDay)*24
     &                         + (header%endHour - header%startHour) + 1
      

C     Position to start of data (immediately after first frame record but before first frame data)
            REWIND (unitNum)
            foundFirstFrame = .false.
            do WHILE(.NOT.foundFirstFrame)
                  read(unit=unitNum, FMT='((A))', iostat=ios) line      ! read a line
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
                   precip_frame_header=line
                  new_precip_frame_header=line
                  precip_flg=.true.
                  if(line(1:1) .eq. ':')then
                       wordCount = SplitLine(line, keyword, subString) ! find the keyword
                        rStat = ToLowerCase(keyword)
                        KeyLen = len_trim(keyword)
                        if(keyword(1:KeyLen).eq.':frame')then
                              iStat = ParseFrameLine(frameRec,keyword,
     &                                              keyLen,subString)
                              foundFirstFrame = .true.
                        end if
                  end if                  
            enddo 


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     changed Dec. 4/07 nk to allow for missing data at start of met file
      backspace unitNum
      read(unitNum,*)junk,modelHour
      framesSkipped=modelHour-1
c           modelHour=1 ! set model hour to first hour
c      print*,junk,modelhour,framesSkipped
c      pause 1
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            dataHour = 0 ! set data hour to just before first hour

C Are we LatLong or Projected
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
                  if(abs(xdelta2-ydelta2).lt.0.001)then  !changed jan21/08 nk
                        rgrd=xdelta2/1000.
                  else
                        print*, 'ERROR - xdelta .ne. ydelta'
	                  print*,' xdelta =',xdelta
	                  print*,' ydelta =',ydelta
                        STOP 'Program aborted in read_rain_ef @ 249'
                  endif
            endif

C firstpass means this is the first met file of the first event...which is the only time it needs to be done
            if(firstpass.eq.'y')then
                  firstpass='n'        
                  xcount_max=max(xcount+1,xcount2)
                  ycount_max=max(ycount+1,ycount2)
                allocate(p(ycount_max,xcount_max),
     *                 sum_precip(ycount_max,xcount_max),
     *                 stat=iAllocate)
c     *                 psmear(na),punused(na),stat=iAllocate)
                       if(iAllocate.ne.0) STOP
     *       'Error with allocation of area6/16a arrays in read_rain_ef'
c                  xcount_max=xcount2
c                  ycount_max=ycount2
            endif

      if(jan.eq.1)then     
C       CALCULATE THE GRID OFFSETS:
            if(iymin.ge.rads.and.jxmin.ge.radw)then
C         COMPARE THE SOUTH-WEST CORNER COORDINATES
                  if(IsLatLong(header%r2cp%csp).eq.(.true.))then
                        iyoffset=int(float(iymin-rads)/rgrdn)+iyshift
                        jxoffset=int(float(jxmin-radw)/rgrde)+jxshift
                  else
                        iyoffset=int(float(iymin-rads)/grdn)+iyshift
                        jxoffset=int(float(jxmin-radw)/grde)+jxshift
                  endif
            elseif(abs(iymin-rads).gt.1.or.abs(jxmin-radw).gt.1)then
!     rev. 9.7.18  Jan.  17/11  - NK: Changed tolerance on the grid check in rear_rain & read_temp
!         WATERSHED GRID IS NOT WITHIN MET GRID
                  print*
                  print*,'WARNING:'
                  print*,' In read_rain_ef:'
                  print*,'Met grid does not coincide with shd file'
                  print*,'This could be due to round off when using'
                  print*,'fractions of degrees as the grid size'
                  print*,'iymin,jxmin,rads,radw/',iymin,jxmin,rads,radw
                  print*,'Please check the headers on the shd & met'
                  print*,'for compatibility'
                  print*                     
                  STOP 'Program aborted in read_rain_ef @ 301'
            else
                  iyoffset=0
                  jxoffset=0
                  ! DO NOTHING
            endif

!     rev. 9.5.16  Feb.  28/08  - NK: moved precip_adjust to sub
c            call precip_adjust()

!       initialize sum_precip
        do i=1,ycount2  
          do j=1,xcount2
            sum_precip(i,j)=0.0
          end do
        end do
      endif   ! END OF JAN=1

C Check for change of grid size
      if(xcount2.gt.xcount_max.or.ycount2.gt.ycount_max)then
            deallocate(p,stat=iDeallocate)
            if(iDeallocate.ne.0) then
              print*,     'Warning: error with deallocation of outarray'
            endif

            xcount_max=xcount2
            ycount_max=ycount2
C       outarray is in areawfo
C       this has to be allocated before calling write_r2c
            allocate(p(ycount_max,xcount_max),stat=iAllocate)
            if(iAllocate.ne.0) then
                  STOP 'Error with allocation of outarray in tr_sub'      
            end if
      endif

!     rev. 10.2.33 Sep.  14/18  - NK: Changed unit=42  fln(12) from clutter to model\*.r2c file '
!           needed when there is only one hour of data in an r2c file      
!           probably should be added to all read r2c s/r's      
            modelHourLast=0

      
            return
      endif          ! (hdrflg.eq.'1')

C***********************************************************
C************************************
C     Finished reading the header 
C************************************   
C***********************************************************

C     We are looking for rain data for this step
      dataHour = dataHour + 1
c        print*,'modelHourLast =',modelHourLast
c        print*,'modelHour =',modelhour,'  ',fln(flnNum)(1:40)
c      if(modelHour.eq.modelHourLast)then
      if(modelHour.eq.modelHourLast)then
        print*,'modelHourLast =',modelHourLast
        print*,'modelHour =',modelhour,' Error in ',fln(flnNum)(1:40)
     *                 ,' met data repeated??'
        write(*,51000)'Error in the met file @ frame #',frameRec%frame
     &            ,frameRec%tStamp%year,'/',frameRec%tStamp%month,'/',
     &            frameRec%tStamp%day,'/',frameRec%tStamp%hour

        write(*,*)'Same time stamp as previous?'
        write(51,*)'Error in the met file @ frame #',frameRec%frame,
     &            frameRec%tStamp%year,'/',frameRec%tStamp%month,'/',
     &            frameRec%tStamp%day,'/',frameRec%tStamp%hour
51000   format(a34,i10,i8,a1,i4,a1,i4,a1,i4) 
        print*
        print*,'Possible cause:'
        print*,'If this file is based on model data, this error can'
        print*,'result from overlapping frames - e.g. last frame in one'
        print*,'monthly CaPA file = same as first frame in next month'
        print*,line(1:62) 
        
        
        
        
        stop 'Program aborted in read_rain @ 379'
        
        
        
        

      endif    
 
C     If there is a record for this hour than modelHour = dataHour

c      print*,'modelHour:',modelHour, 'dataHour:',dataHour
c      print*,'modelHour:',modelHour, 'spl hour:',jz
c      pause


c     if(modelHour.eq.dataHour)then

c      print*,'modelHour,jz',modelHour,jz

      if(modelHour.eq.jz)then
!           the modelHour is found after reading the previous data. 
!           I.e. it is the hour
!           of the next set of data. So here we just wait until 
!           we get the next data set.
C           Go ahead and read the data for this frame 
            do i=1,ycount2
!                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  read(unitNum,*,iostat=ios)(p(i,j),j=1,xcount2)
99955             format(999f6.1)                  
!                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  if(ios.ne.0)then
                        write(*,9993)modelHour,i
                        print*,' last data read:'
                        print*,' jradw= 1 jrade= ',jrade
                        print*,' iradn= ',iradn,' irads= 1'
                        STOP ' program aborted - read_rain_ef'
                  endif
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
                  precip_frame_header=new_precip_frame_header
                  precip_flg=.true.
            end do
!           sum_precip as read in
!           see also addition of snow in read_gsn.f
            do i=1,ycount2  
              do j=1,xcount2
                 sum_precip(i,j)=sum_precip(i,j)+p(i,j)*conv
              end do
            end do

C           Read the next frame. If we are at the end of the file...close
            foundFrame = .false.
            do WHILE(.NOT.foundFrame)
                  read(unit=unitNum, FMT='((A))', iostat=ios) line      ! read a line
c                  print*,'rain ',line(1:60)
                  if(EOF(unitNum)) then 
                        close(unit=unitNum,status='keep',iostat=ios)
                        RETURN
                  endif
                  if(line(1:1) .eq. ':')then
                        wordCount = SplitLine(line, keyword, subString) ! find the keyword
                        rStat = ToLowerCase(keyword)
                        KeyLen = len_trim(keyword)
                        if(keyword(1:KeyLen).eq.':frame')then
                              iStat = ParseFrameLine(frameRec,keyword,
     &                                    keyLen,subString)
                              foundFrame = .true.
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
                              new_precip_frame_header=line
                              new_precip_flg=.true.
                        end if
                  end if            
                  
            enddo 

            modelHourLast=modelHour

            yearLLL=frameRec%tStamp%year
            monthLLL=frameRec%tStamp%month
            dayLLL=frameRec%tStamp%day
            hourLLL=frameRec%tStamp%hour
            
C           Determine the next modelhour        
            modelHour = 
     &           (JDATE(frameRec%tStamp%year,frameRec%tStamp%month,
     &            frameRec%tStamp%day) - header%startJulianDay) * 24
     &            + (frameRec%tStamp%hour - header%startHour) + 1
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 changed Dec. 4/07 nk to allow for missing data at start of met file
     &            +framesSkipped
!                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!           check added May 2/12  nk
!     rev. 10.2.24 May   21/18  - NK: Added error message in Read_rain & read_tmp
            if(modelHour-modelHourLast.lt.0)then
              print*,'modelHourLast =  ',modelHourLast
              print*,'modelHour =      ',modelHour
              print*,'headerStartHour =',header%startHour
              print*,'Note: Frame date & time must be monotomically '
              print*,'increasing in an r2c time series file.'
              print*,'Stopped at:'
              print*,line(1:72)
              print*,'year  ',frameRec%tStamp%year
              print*,'month ',frameRec%tStamp%month
              print*,'day   ',frameRec%tStamp%day
              print*,'hour  ',frameRec%tStamp%hour
              print*,'Previour frame @:'
              print*,'year  ',yearLLL
              print*,'month ',monthLLL
              print*,'day   ',dayLLL
              print*,'hour  ',hourLLL
              print*
              print*,'possible cause: month repeated in ',
     *                                     fln(flnNum)(1:40)
              print*
              stop 'Program aborted in read_rain_ef @ L411'
            endif
c           print*,'in read_rain modelHour=',modelHour
c	      print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      else                          !IF(K.EQ.1)...
c           print*,'                                skipped read'
!           No data this hour
!           Nothing was read in this call in this time step
!           The model is continuing with zero precip until 
!           an hour with precip is found
            do i=1,ycount2  
              do j=1,xcount2
                p(i,j)=0.0
              end do
            end do
      endif
      RETURN

!     FORMATS

 9993 format(' error reading precip at hour/line ',2i10/)
      
      END SUBROUTINE read_rain_ef
