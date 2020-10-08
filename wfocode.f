       SUBROUTINE write_both_headers(coordsys,jan)

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
!     REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90
!     rev. 9.1.44  Jun.  11/03  - Added Cumulative precip to the wfo file
!     rev. 9.9.02  Dec.  12/13  - NK: Changed format for origin in wfo code

!***********************************************************************


      use area_watflood
	implicit none


!  INPUT PARAMETERS: FILENAME FOR WFO FILE, UNIT NUMBER, NUMER OF 
!                    COLUMNS, NUMBER OF ROWS, GRID SIZE, GRID ORIGINS
!                    (X,Y), NUMBER OF LANDCOVER TYPES


!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

        CHARACTER(64):: appname
        CHARACTER(2) :: num_to_char(16)
        CHARACTER(*) :: coordsys
!        INTEGER      :: xcount,ycount,classcount,i,j,attcount,seq,step,
!     *          nsteps,attnum,iAllocate,iDeallocate,jan,nrvr,nj,npick
        INTEGER      :: i,j,attcount,seq,step,
     *          nsteps,attnum,iAllocate,iDeallocate,jan,nj,npick
!        REAL         :: xdelta,ydelta,xorigin,yorigin

!        CHARACTER(64), DIMENSION(:), ALLOCATABLE :: attname
!        CHARACTER(32), DIMENSION(:), ALLOCATABLE :: attunits

      DATA num_to_char/'01','02','03','04','05','06','07','08',
     *                 '09','10','11','12','13','14','15','16'/
      
!       WFO IO FUNCTIONS
        INTEGER wfo_open_file,wfo_write_header
        INTEGER wfo_write_attribute_header
        INTEGER wfo_close_header

      appname = 'Watflood'
!     coordsys='UTM'       !UTM or LATLONG

!     SAMPLE ATTRIBUTES
!     SET UP ATTRIBUTE INFORMATION (NAMES AND UNITS)

!     Note:  These are also read from the wfo_spec.txt file
!     where the names are just for information. They are hard wired
!     here so they can not be changed accidentally by the user.

      j=0
      if(wfo_pick(1).eq.1)then
        j=j+1
        attname(j)='Temperature'
        attunits(j)='degrees Celsius'
      endif

      if(wfo_pick(2).eq.1)then
        j=j+1
        attname(j)='Precipitation'
        attunits(j)='mm'
      endif

      if(wfo_pick(3).eq.1)then
        j=j+1
        attname(j)='Cumulative Precipitation'
        attunits(j)='mm'
      endif

      if(wfo_pick(4).eq.1)then
        j=j+1
        attname(j)='Lower Zone Storage'
        attunits(j)='mm'
      endif

      if(wfo_pick(5).eq.1)then
        j=j+1
        attname(j)='Lower Zone Discharge'
        attunits(j)='m^3/s'
      endif

      if(wfo_pick(6).eq.1)then
        j=j+1
        attname(j)='Grid Runoff'
        attunits(j)='Cms'
      endif

!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
      if(wfo_pick(7).eq.1)then
        j=j+1
        attname(j)='Observed Grid Outflow'
        attunits(j)='Cms'
      endif

      if(wfo_pick(8).eq.1)then
        j=j+1
        attname(j)='Computed Grid Outflow'
        attunits(j)='Cms'
      endif

      if(wfo_pick(9).eq.1)then
        j=j+1
        attname(j)='Weighted SWE'
        attunits(j)='mm'
      endif

!     rev. 9.1.35  Dec. 26/02  - Added wetland & channel height to the wfo file
      if(wfo_pick(10).eq.1)then
        j=j+1
        attname(j)='Wetland Internal Depth'
        attunits(j)='meters'
      endif

      if(wfo_pick(11).eq.1)then
        j=j+1
        attname(j)='Channel Depth'
        attunits(j)='metres'
      endif

!     rev. 9.1.21  Jun.  28/02  - Added wetland storage & outflow to the wfo file
      if(wfo_pick(12).eq.1)then
        j=j+1
        attname(j)='Wetland Storage'
        attunits(j)='cubic meters'
      endif

      if(wfo_pick(13).eq.1)then
        j=j+1
        attname(j)='Wetland Outflow'
        attunits(j)='Cms'
      endif

      if(wfo_pick(14).eq.1)then
        j=j+1
        attname(j)='Weighted Cumm ET'
        attunits(j)='mm'
      endif
      
!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
      if(wfo_pick(15).eq.1)then
        j=j+1
        attname(j)='Bankfull'
        attunits(j)='Percent'
      endif
      
!     rev. 9.9.70  Jun.  12/15  - NK: Add del_rain, and dSTRconc2  to the wfo file
!     rev. 9.9.70  Jun.  12/15  - NK: Add del_rain, and dSTRconc2  to the wfo file
      if(frcflg.eq.'y')then
c      if(wfo_pick(12).eq.1)then
        j=j+1
        attname(j)='del_rain'
        attunits(j)='‰'
c      endif
c      if(wfo_pick(13).eq.1)then
        j=j+1
        attname(j)='dSTRconc2'
        attunits(j)='‰'
c      endif
      endif
      

      npick=15         ! change this if items are added above
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~` 

      do i=1,classcount
        if(wfo_pick(npick+i).eq.1)then
          j=j+1
          attname(j)='Depression Storage Class '//num_to_char(i)
          attunits(j)='mm'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Depression Storage (Snow) Class '//
     *                       num_to_char(i)
          attunits(j)='mm'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+2*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Snow Water Equivalent Class '//
     *                         num_to_char(i)
          attunits(j)='mm'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+3*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Snow Covered Area Class '//
     *                         num_to_char(i)
          attunits(j)='fraction'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+4*(classcount)+i).eq.1)then
          j=j+1
!     rev. 10.4.21 Apr.  21/20  = NK Add UZS deficit to wfo file = UZS(class=classcount)
          if(i.ne.classcount)then
            attname(j)='Upper Zone Storage Class '//
     *                          num_to_char(i)
          else
            attname(j)='Upper Zone Storage Deficit '//
     *                          num_to_char(i)
          endif    
          attunits(j)='mm'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+5*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Upper Zone Storage (Snow) Class '//
     *                         num_to_char(i)
          attunits(j)='mm'
        endif
      end do

!     rev. 9.1.30  Nov.  08/02  - added q1, qint, drng & qlz to the wfo file
      do i=1,classcount
        if(wfo_pick(npick+6*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Surface Flow (bare) Class '//
     *                         num_to_char(i)
          attunits(j)='m^3/s'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+7*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Surface Flow (snow) Class '//
     *                         num_to_char(i)
          attunits(j)='m^3/s'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+8*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Interflow (bare) Class '//
     *                         num_to_char(i)
          attunits(j)='m^3/s'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+9*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Interflow (snow) Class '//
     *                         num_to_char(i)
          attunits(j)='m^3/s'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+10*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Recharge Class '//
     *                         num_to_char(i)
          attunits(j)='mm'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+11*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Recharge (snow) Class '//
     *                         num_to_char(i)
          attunits(j)='mm'
        endif
      end do

      do i=1,classcount
        if(wfo_pick(npick+12*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Pot. Evapotranspiration Class '//
     *                         num_to_char(i)
          attunits(j)='mm'
        endif
      end do
      do i=1,classcount
        if(wfo_pick(npick+13*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Evapotranspiration Class '//
     *                         num_to_char(i)
          attunits(j)='mm'
        endif
      end do
      do i=1,classcount
        if(wfo_pick(npick+14*(classcount)+i).eq.1)then
          j=j+1
          attname(j)='Sublimation Class '//
     *                         num_to_char(i)
          attunits(j)='mm'
        endif
      end do

!     attcount redefined here
      attcount=j

!        OPEN THE FILE FOR WRITING
      if(wfo_open_file().ne.1)then
        write(*,'(A)') ' '
        write(*,'(A)') '  *** FATAL ERROR ***     '
        write(*,'(A)') ' Unable to Open File for Writing'

        write(*,'(9X,(A))') filename(65)
        STOP 'Program terminated in write_both_headers @89'
      endif

!     WRITE THE MAIN FILE HEADER


      if(wfo_write_header(appname,
     *               coordsys) .ne. 1)then

         write(*,'(A)') ' '
         write(*,'(A)') '  *** FATAL ERROR ***     '
         write(*,'(A)') ' Unable to Write Main Header'
         STOP 'Program terminated in write_both_headers @98'
      endif

!     WRITE THE ATTRIBUTE HEADER
      if(wfo_write_attribute_header(attcount,attname,attunits) 
     *                       .ne. 1)then
         write(*,'(A)') ' '
         write(*,'(A)') '  *** FATAL ERROR ***     '
         write(*,'(A)') ' Unable to Write Attribute Header'
         STOP 'Program terminated in write_both_headers @107'
      endif	

!     CLOSE THE HEADER
      if(wfo_close_header().ne. 1)then
         write(*,'(A)') ' '
         write(*,'(A)') '  *** FATAL ERROR ***     '
         write(*,'(A)') ' Unable to Close Header'
         STOP 'Program terminated in write_both_headers @115'
      endif

! DEALLOCATIONS OF ARRAYS:
      deallocate(attname,attunits,stat=iDeallocate)
      if (iDeallocate.ne.0) STOP   
     *    'Error with deallocation of ensim arrays in wfocodea'

      RETURN

      END SUBROUTINE write_both_headers

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

!***********************************************************************
      INTEGER FUNCTION wfo_open_file()
!***********************************************************************
!    REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90

! THIS FUNCTION OPENS A FILE AS BINARY

! RETURN:        = 1 SUCCESS
!                = 0 ERROR
!
!***********************************************************************

      use area_watflood
	implicit none

!     PARAMETER TYPE DEFINITIONS
!        CHARACTER(*) :: fln
!        INTEGER      :: lun
      integer*4     ::  ios
      logical       :: exists

      wfo_open_file = 0


!       OPEN THE FILE
        open(unit=65,file=filename(65),status='unknown',
     *                form='binary',iostat=ios)
	  if(ios.ne.0)then
	    print*,' Problems opening unit 65'
	    print*,' Problem ignored. Look for file name for65'
	    print*
	  else
          print*,' Opened the wfo file for ENSIM'
          wfo_open_file = 1
	  endif

      RETURN

      END FUNCTION wfo_open_file 


!***********************************************************************
!      INTEGER FUNCTION wfo_write_header(appname, xorigin, yorigin,
!     *	               xcount, ycount, xdelta, ydelta, coordsys)
      INTEGER FUNCTION wfo_write_header(appname,coordsys) 
!***********************************************************************
!    REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90

! THIS FUNCTION WRITES A PROPER WFO FILE HEADER

! RETURN:       = 1 SUCCESS
!               = 0 ERROR
!ed basin\
!***********************************************************************

!        USE DFLIB

      use area_watflood
	implicit none

!     PARAMETER TYPE DEFINITIONS
        CHARACTER(*) :: appname,coordsys
!        INTEGER      :: xcount,ycount
!        REAL         :: xorigin,yorigin,xdelta,ydelta

!     LOCAL VARIABLES
        CHARACTER(128) :: line
        CHARACTER(64)  :: username
        CHARACTER(10)  :: ctime
        CHARACTER(8)   :: fmt_string,cday
        CHARACTER(5)   :: czone
        CHARACTER(4)   :: data_type
        INTEGER        :: mi,ss,ms,ivalues(8)

      print*,'Writing the wfo header for ENSIM' 
      wfo_write_header = 1

      fmt_string = 'BINARY'
      data_type = 'wfo'

!     WRITE DATA
      call wfo_write_header_record(
     *'###########################################################')

      write(line,101) data_type, fmt_string
  101 format(':FileType ',(A),' ',(A),' EnSim 1.0')

      call wfo_write_header_record(line)
      call wfo_write_header_record('#')
      call wfo_write_header_record(
     *	'# DataType       Binary Watflood Output')
      call wfo_write_header_record('#')

      write(line,102) appname
  102 format(':Application       ', (A))

      call wfo_write_header_record(line)

!      call getlog(username)
!      write(line,103) username
! 103 format(':WrittenBy         ', (A))
!      call wfo_write_header_record(line)

!     GETDAT AND GETTIM REPLACED BY DATE_AND_TIME INTRINSIC
!      call getdat(yy,mm,dd)
!      call gettim(hh,mi,ss,ms)
      call date_and_time(cday,ctime,czone,ivalues)
      yy=ivalues(1)
      mm=ivalues(2)
      dd=ivalues(3)
      hh=ivalues(5)
      mi=ivalues(6)
      ss=ivalues(7)
      ms=ivalues(8)

      write(line,104) yy,mm,dd,hh,mi,ss
  104 format(':CreationDate        ',I4.4,"/",I2.2,"/",I2.2,2x,
     *  I2.2,":",I2.2,":",I2.2)
      call wfo_write_header_record(line)
      
      write(line,114) year_start,mo_start,day_start,hour_start,0,0
  114 format(':StartTime ',I4.4,"/",I2.2,"/",I2.2,2x,
     *  I2.2,":",I2.2,":",I2.2)
      call wfo_write_header_record(line)
     
      call wfo_write_header_record(          
     *'#-----------------------------------------------------------')
      call wfo_write_header_record('#')

!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
      write(line,1005) -999.0                !noDataValue
 1005 format(':NoDataValue ', (f))
      call wfo_write_header_record(line)

      write(line,105) coordsys1
  105 format(':Projection ', (A))
      call wfo_write_header_record(line)

!     REV. 10.1.37 Jul   28/16  - NK: Added "Ellipsoid to the WFO header
      write(line,1006) datum1
 1006 format(':Ellipsoid ', (A))
      call wfo_write_header_record(line)

      write(line,106) xorigin
!     rev. 9.9.02  Dec.  12/13  - NK: Changed format for origin in wfo code
c  106 format(':xorigin ',(f))
  106 format(':xorigin ',(e))
      call wfo_write_header_record(line)

      write(line,107) yorigin
!     rev. 9.9.02  Dec.  12/13  - NK: Changed format for origin in wfo code
c  107 format(':yorigin ',(f))
  107 format(':yorigin ',(e))
      call wfo_write_header_record(line)

      write(line,108) xcount
  108 format(':xcount ',(i))
      call wfo_write_header_record(line)

      write(line,109) ycount
  109 format(':ycount ',(i))
      call wfo_write_header_record(line)

      write(line,110) xdelta
  110 format(':xdelta ',(f))
      call wfo_write_header_record(line)

      write(line,111) ydelta
  111 format(':ydelta ',(f))
      call wfo_write_header_record(line)

      RETURN

      END FUNCTION wfo_write_header


!***********************************************************************
      SUBROUTINE wfo_write_header_record(line)
!***********************************************************************
!    REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90

! THIS SUBROUTINE WRITES A HEADER RECORD TO THE OPEN LUN
! APPEND A NEWLINE CHAR SINCE THE FILE IS BINARY
!
!***********************************************************************

		
!     PARAMETER TYPE DEFINITIONS
        CHARACTER(*) :: line   	
!       can't use (*) cause we are adding something to the end
!        INTEGER :: lun
      integer*4    :: ios

!     LOCAL VARIABLES
        CHARACTER(256) :: lline
        INTEGER :: llen

      llen = LEN_TRIM(line)
      lline = line(1:llen)// ' '               ! Append a space

      llen = LEN_TRIM(lline)
      if(llen .GT. 255)then
         llen = 255
      end if

      llen = llen+1
      lline(llen:llen) = CHAR(10)               ! Append a line feed

      write(65,iostat=ios) lline(1:llen)
	if(ios.ne.0)then
	  print*,' Error writing to unit 65 at approx line 310 in wfocode'
	  print*
	  stop 'in wfocode'
	endif


      RETURN

      END SUBROUTINE wfo_write_header_record


!***********************************************************************
      INTEGER FUNCTION wfo_write_attribute_header(attcount,
     *                 attname,attunits)
!***********************************************************************
!    REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90

! THIS FUNCTION WRITE THE ATTRIBUTE NAMES AND UNITS INTO THE HEADER  

! RETURN:       = 1 SUCCESS
!               = 0 ERROR
!
!***********************************************************************

!        USE DFLIB

!     PARAMETER TYPE DEFINITIONS
        INTEGER       :: attcount
        CHARACTER(64) :: attname(attcount)
        CHARACTER(32) :: attunits(attcount)


!     LOCAL VARIABLES
        CHARACTER(128) :: line
        INTEGER        :: attnum

      wfo_write_attribute_header = 1
      call wfo_write_header_record('#')

      write(line,121) attcount
  121 format(':AttributeCount ',(i))
      call wfo_write_header_record(line)

      do attnum=1, attcount
         write(line,122) attnum, attname(attnum)
         call wfo_write_header_record(line)
         write(line,123) attnum, attunits(attnum)
         call wfo_write_header_record(line)
      end do
  122 format(':AttributeName ',(i),' ',(A))
  123 format(':AttributeUnits ',(i),' ',(A))

      call wfo_write_header_record('#')

      RETURN

      END FUNCTION wfo_write_attribute_header


!***********************************************************************
      INTEGER FUNCTION wfo_write_timestamp(seq,step,yy,mm,dd,
     *                 hh,mi,ss,ms)
!***********************************************************************
!    REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90

! THIS FUNCTION WRITES A BINARY TIMESTAMP FOR THIS STEP
! YR, MONTH, DAY, HR, MiN, SEC, MiLLISEC FOR THIS STEP
 
! RETURN:       = 1 SUCCESS
!               = 0 ERROR
!
!***********************************************************************


!     PARAMETER TYPE DEFINITIONS
        INTEGER :: seq,step
        INTEGER :: yy,mm,dd,hh,mi,ss,ms,ios

      wfo_write_timestamp = 1

!     WRITE THE RECORD

      write(65,iostat=ios) seq,step,yy,mm,dd,hh,mi,ss,ms
	if(ios.ne.0)then
	  print*,' Error writing to unit 65 at approx line 391 in wfocode'
	  print*,' Disk full maybe?'
	  print*
	  stop 'in wfocode'
	endif

      RETURN

      END FUNCTION wfo_write_timestamp


!***********************************************************************
      INTEGER FUNCTION wfo_close_header()
!***********************************************************************
!    REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90

! THIS FUNCTION CLOSES THE HEADER

! RETURN:       = 1
!
!***********************************************************************


!     PARAMETER TYPE DEFINITIONS
!        INTEGER :: lun

      wfo_close_header=1

      call wfo_write_header_record('#')

      call wfo_write_header_record(':EndHeader')

      RETURN

      END FUNCTION wfo_close_header


!***********************************************************************
!      INTEGER FUNCTION wfo_write_attribute_data(nx,ny,arr,ensimflg)
      INTEGER FUNCTION wfo_write_attribute_data(nx,ny)
!***********************************************************************
!    REV. 9.00   Mar.  2000  -  TS: CONVERTED TO FORTRAN 90

! THIS FUNCTION WRITES A BINARY WFO DATA RECORD
! WRITES GRID DATA FROM LEFT TO RIGHT STARTING FROM BOTTOM ROW
 
! RETURN:       = 1 SUCCESS
!               = 0 ERROR
!
!***********************************************************************

      use area_watflood
	implicit none

!     PARAMETER TYPE DEFINITIONS
      
      INTEGER :: nx,ny,i,j,ios

!      integer :: nx,ny,i,j
!      REAL*4    :: arr(nx,ny)

c      character*1  :: ensimflg

!       temp fix until allocated properly


      wfo_write_attribute_data = 1

!     WRITE THE RECORD
!	write(65,iostat=ios) ((arr(I,J),I=1,nx),J=1,ny)
      write(65,iostat=ios) ((outwfo(i,j),i=1,nx),j=1,ny)
	if(ios.ne.0)then
	  print*,' Error writing to unit 65 at approx line 458 in wfocode'
	  print*
        ensimflg='n'
	  
	  write(51,5101)
        write(98,5101)
5101  format(' Error writing to unit 65 at approx line=458 in wfcode'/
     * ' File closed and progarm continues without writing the rest'/
     * ' of the watflood.wfo file')
	endif

      RETURN

      END FUNCTION wfo_write_attribute_data
