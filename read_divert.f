      SUBROUTINE read_divert(unitNum,flnNum)

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
     
     
!*****************************************************************************
!  read_divert_ef - written Jan/09 by Nick Kouwen, UW
!     - Derived from rdresv written by Dave Watson CHC
!     This subroutine reads the diversion (div) file 
!     (tb0 format)
!*****************************************************************************

!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
!     rev. 9.5.81  Jan.  05/11  - NK: allow reservoirs outside watershed in resv file
!     rev. 9.7.24  Apr.  20/11  - NK: Added diverflg to indicate if a diversion is in grid
!******************
      use area_watflood

! TB0 data module
      USE EF_Module

      implicit none
      TYPE(DivParam) :: header
      TYPE(DivColumnMetaData) :: colHeader

      Integer  :: ios,j,k,i,n,l,jz,jj,ng,nt
      integer  :: nodivert_firstpass,ndiv_max,
     *            iAllocate,iDeallocate,
     *            ntake,ngive
      real*4   :: factor

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      LOGICAL exists,diversion_local
      CHARACTER(1)  :: firstpass
      

      data firstpass/'y'/

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

! Parameter type definitions
      integer*4 unitNum, flnNum, iStat

! Local variables
      character*4096 line, subString, tmpString
      character*128 keyword, value 
      character*12  outlet_type
      integer lineLen, keyLen, wordCount
      logical rStat, lineType, foundEndHeader, insideColMetaData

! Open the file
      INQUIRE(FILE=fln(flnNum),EXIST=exists)
      if(exists)then
            open(unit=unitNum,file=fln(flnNum),status='old',iostat=ios)
	      if(iopt.ge.1)print*,'opened unit=',unitNum,fln(flnNum)(1:40)
            if(ios.ne.0)then
                  print*,'Problems opening ',fln(flnNum)
                  print*
                  STOP ' Stopped in read_divert_ef'
            endif
        nodivert=1   ! assume there is at least one diversion if there is a file
      else
        if(numa.eq.0)then
          print*,'WARNING:'
          print*,'Diversion (div) file ',fln(flnNum)(1:40)
          print*,'is NOT found'
          print*,'Program continues with no diversions'
          print*
        endif
        nodivert=0
        ndiv=0
      endif

d      if(iopt.eq.2)print*,'in read_divert_ef passed 70'
cd      if(iopt.eq.2)pause 'in read_divert @ 70'


      if(nodivert.ge.1)then


! Initialize default values
        CALL InitDivParam(header)   

d       if(iopt.ge.2)print*,'in read_divert_ef passed 76'
cd       if(iopt.ge.2)pause 'passed 76'

! Search for and read tb0 file header
        line(1:1) = '#'
        foundEndHeader = .false.
        insideColMetaData = .false.

        do WHILE((.NOT.foundEndHeader) .AND.
     &        ((line(1:1) .eq. '#') .OR.
     &          (line(1:1) .eq. ':') .OR.
     &          (LEN_TRIM(line) .eq. 0)))     

          read(UNIT=unitNum, FMT='((A))', iostat=ios) line      ! read a line
c      print*,line(1:72)
          if(ios .eq. -1)then
             write(6,'((A))') 'ERROR: Premature EndOfFile encountered'
             STOP ' Stopped in read_divert_ef'
          end if

          rStat = Detab(line)                       ! replace tabs with spaces
          line = ADJUSTL(line)          ! Get rid of leading white space
          lineLen = LEN_TRIM(line)            ! Find the length excluding trailing spaces

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

!             this is not needed I think as there are no coefficients
             else if(insideColMetaData) then
                iStat = ParseDivColumnMetaData
     &                      (colHeader,keyword,keyLen,subString)
                if(iStat .lt. 0) then
                   write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                   write(*,'(2(A))') '   in line: ',line                             
                   STOP ' Stopped in read_divert_ef'
                   return
                endif
             else
                iStat = ParseDivParam(header,keyword,keyLen,subString)
                if(iStat .lt. 0) then
                      write(*,'(2(A))') 'ERROR parsing ', fln(flnNum)
                      write(*,'(2(A))') '   in line: ',line                             
                      STOP ' Stopped in read_divert_ef'
                      return
                else if(iStat .eq. 0) then
!                     write(*,'((A), (A))')  'Unrecognized keyword line: ',
!     &                                                   line
                endif
             end if
          end if
        end do    
!***************************************
!       Finished reading the header
!***************************************

d      if(iopt.ge.2)print*,'in read_divert_ef passed 132'
cd      if(iopt.ge.2)pause 'passed 132'


! Assign the variables from the types
        ktr =  header%tb0p%deltaT    !data time step in hours
!     rev. 9.8.54  Apr.  02/13  - NK: deltat conversion seconds to hours
!       In the past, deltat has bee in hours but GK wants them in seconds
!       This converts deltat to hours as needed in spl if a large deltat is found        
        if(ktr.gt.1000)then
          ktr=ktr/3600
        endif
        factor = header%tb0p%unitConv ! conversion to cms
        nodivert = colHeader%tb0cmd%colCount !no of reservoirs

!     ndiv    =     no of hours of data
! Scan lines of data
        rewind unitNum

!        added *ktr     Nov. 27/06  nk
        ndiv = CountDataLinesAfterHeader(unitNum)*ktr
        rewind unitNum
        CALL GoToStartOfData(unitNum)

      endif

!       allocate stuff      
      if(firstpass.eq.'y')then
        nodivert_firstpass=nodivert
        diversion_local=.true.
!       nl comes from the .str file and is the # hour of the event

c        ndiv_max=max0(ndiv,nl)
        if(ndiv.gt.nl)then
          print*,'div file longer than the str file'
          print*,'data past ',nl,' hours ignored'
        endif

!       but we need to provide enough memory to simulate a whole event
!       sometimes users specify the duration in the rel to be just 1 hr.
!       when a rule is given. However, we need memory of all the variables
        if(nodivert.gt.0)then
          allocate(divertname(nodivert),
     *    xtake(nodivert),ytake(nodivert),
     *    xgive(nodivert),ygive(nodivert),
     *    jtake(nodivert),itake(nodivert),
     *    jgive(nodivert),igive(nodivert),
     *    gridtake(nodivert),gridgive(nodivert),
     *    qdivert(nodivert,ndiv*ktr),
     *    divert_inbasin_flg(nodivert),
     *    upstrda(nodivert),stat=iAllocate)

          if(iAllocate.ne.0)then
            print*,'Error with allocation in read_divert_ef @172'
            print*,'firstpass =',firstpass
            STOP 'Program aborted in read_divert @ 198'
          endif
!         initialize upstrda
          do l=1,nodivert
            upstrda(l)=0.0
          end do
                    
        endif
!     rev. 9.9.43  Nov.  26/14  - NK: Allocation for divertflg = 'g' 
        if(divertflg.eq.'g')return   ! ' only here to allocate




!       ASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
!       this is done only during the first pass if coefficient values 
!       are set to -1 for subsequent events. This makes tweaking easy
!       as only the values in the first event need to be adjusted.
        do n=1,nodivert
          divertname(n) = colHeader%tb0cmd%colName(n) ! diversion name
          xtake(n) = colHeader%tb0cmd%colLocX(n) ! x coordinate
          ytake(n) = colHeader%tb0cmd%colLocY(n) ! y coordinate
          xgive(n) = colHeader%tb0cmd%colLocX1(n) ! x coordinate
          ygive(n) = colHeader%tb0cmd%colLocY1(n) ! y coordinate
!     rev. 9.9.55  Jan.  22/15  - NK: Added diversion upstream drainage area in div file
          if(allocated(colHeader%colValue1))then
!           means upstream drainage area is present          
            upstrda(n) = colHeader%colValue1(n)     ! upstream drainage area
          else
            upstrda(n)=0.0
            print*,
     *          'Warning: upstream drainage area at diversion not found'
          endif
    
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
!     check to see if diversion locations are in the watershed
!     first see if it is in the grid
          if(xtake(n).lt.xorigin.or.
     *       xtake(n).gt.xorigin+xcount*xdelta.or.
     *       ytake(n).lt.yorigin.or.
     *       ytake(n).gt.yorigin+ycount*ydelta)then
            print*
            print*,'WARNING:'
            print*,'Diversion origin and/or destination'
            print*,'not in watershed'
            print*,'Diversion ',n,' is ignored'
            print*
            write(51,*)'WARNING:'
            write(51,*)'diversion origin and/or destination'
            write(51,*)'not in watershed'
            write(51,*)'diversion ',n,' is ignored'
          endif

!     rev. 9.9.41  Nov.  20/14  - NK: Added check if diversion = in-basin
          write(51,*) 
          write(51,*)'Check if diversion locations are in the grid'
          write(51,*)'There are ',na,' grids' 
          write(51,*)
	    divert_inbasin_flg(n)=.false.
          i=int((ytake(n)-yorigin)/ydelta)+1
	    j=int((xtake(n)-xorigin)/xdelta)+1
	    write(51,*)'i,ycount,j,xcount',i,ycount,j,xcount
		  if(i.ge.1.and.i.le.ycount.and.j.ge.1.and.j.le.xcount)then
      	    if(s(i,j).ge.1.and.s(i,j).le.na)then
              gridtake(n)=s(i,j)
              write(51,*)'gridtake,n,i,j',gridtake(n),n,i,j,s(i,j)
            else
              gridtake(n)=-1
            endif
          else
            gridtake(n)=-1  
          endif
c          gridtake(n)=s(i,j)
          i=int((ygive(n)-yorigin)/ydelta)+1
          j=int((xgive(n)-xorigin)/xdelta)+1
	    if(i.ge.1.and.i.le.ycount.and.j.ge.1.and.j.le.xcount)then
      	    if(s(i,j).ge.1.and.s(i,j).le.naa)then
              gridgive(n)=s(i,j)
            else
              gridgive(n)=-1
            endif
          else
            gridgive(n)=-1  
          endif
          if(gridtake(n).ge.1.and.gridgive(n).ge.1)then
            divert_inbasin_flg(n)=.true.
          endif
         write(51,*)n,divert_inbasin_flg(n)
         write(51,*)xgive(n),ygive(n),n,'gridgive,x,y,n'
         write(51,*)
         write(51,*)
      end do

      
!     rev. 9.9.55  Jan.  22/15  - NK: Added diversion upstream drainage area in div file
!     REV. 10.1.40 Oct   11/16  - NK: Fixed bug in read_divert for missing u/s DA
        do l=1,nodivert
          n=gridgive(l)
d         print*,'source grid =',l,' in grid # ',n          
!         add the u/s area to the grid receiving the diversion
          if(n.ge.1)then
d            print*,n,da(n),upstrda(l),da(n)+upstrda(l)
            da(n)=da(n)+upstrda(l)
          endif
d          print*,'naa=',naa,'  in read_divert @ 309'
          i=0
!     rev. 10.2.19 Mar.  13/18  - NK: Fixed array fault read_divert.f
c          do while(n.lt.naa.and.n.gt.1)
          if(n.lt.naa.and.n.gt.1)then
            if(next(n).ge.1)then
              da(next(n))=da(next(n))+upstrda(l)
d              print*,n,next(n),da(next(n))
              n=next(n)
            endif
          endif
c          end do
c          do i=1,naa
c            if(next(n).ge.1)then
c              da(next(n))=da(next(n))+upstrda(l)
cd              print*,l,n,next(n),da(next(n))
c              n=next(n)
c            endif
c          end do
        end do
      
        firstpass='n'

      endif      !   firstpass
      
d      if(iopt.ge.2)print*,'End first pass in read_divert'
cd      if(iopt.ge.2)pause 'at 320'
      
!     rev. 10.1.85 May   17/17  - NK: Level_station_location.xyx for iopt > 0 only
!     rev. 10.1.86 May   17/17  - NK: Diversion_location.xyx for iopt > 0 only
      if(iopt99)then
      open(unit=99,file='diversion_location.xyz',status='unknown')
        do n=1,nodivert
          write(99,*)xtake(n),ytake(n),n
          write(99,*)xgive(n),ygive(n),n
        end do
        close(unit=99,status='keep')
      endif

c	if(nodivert_firstpass.eq.0)return
	if(.not.diversion_local)then
!       rule: if no diversion in the first event, then not later either
	  diversion=.false.
	  print*,'WARNING  <<<<<<<<<<<<<<<<<<<<'
	  print*,'A diversion file was found but the data can not be used'
	  print*,'for some reason'
	  print*,'Possible problem: locations not in watershed grids'
	  print*,'Program continues without diversion flows'
	  print*
	  return
	endif

!     subsequent passes
!     check to see memory allocated is adequate      
d      if(iopt.eq.2)print*,'In read_divert_ef @ 258'
      if(nodivert.ne.nodivert_firstpass)then
	  print*
        print*,'No of diversions has been changed in'
        print*,'in file ',fln(7)
        print*,'This is not permitted'
        print*
        stop 'Program aborted in read_divert @ 264'
      endif
d      if(iopt.eq.2)print*,'In read_divert_ef @ 266'

      if(nl.gt.ndiv_max)then
        ndiv_max=max0(ndiv,nl)
!       the event is longer than any of the previous events so 
!       more memory has to be allocated
!       DEALLOCATION OF ARRAYS FROM AREA10A:

d       if(iopt.eq.2)print*,'in read_divert_ef @ 213'

!       DEALLOCATION OF ARRAYS FROM AREA5A:
        deallocate(qdivert,stat=iDeallocate)     
        if (iDeallocate.ne.0) print*,    
     *    'Error with deallocation of area5a arrays'

!       re-allocate for larger arrays
        allocate(qdivert(nodivert,ndiv_max*ktr),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Error with allocation of area10a arrays in read_divert'
!     rev. 9.5.23  Mar.  12/08  - NK: fixed allocation error in read_divert_ef
          endif

!       REASSIGN THE GAUGE LOCATIONS NAMES AND FUNCTION COEFFICIENTS
!       This part is used only if coefficient values area >0
!       Used if values change over time & need to be reassigned.
        do i=1,nodivert
          divertname(i) = colHeader%tb0cmd%colName(i) ! reservoir name
          xtake(i) = colHeader%tb0cmd%colLocX(i) ! x coordinate
          ytake(i) = colHeader%tb0cmd%colLocY(i) ! y coordinate
          xgive(i) = colHeader%tb0cmd%colLocX1(i) ! x coordinate
          ygive(i) = colHeader%tb0cmd%colLocY1(i) ! y coordinate
        end do

!       rev. 9.1.69  Dec.  19/04  - NK: rewrote rdresv c/w memory allocation 

d      if(iopt.eq.2)print*,'In read_divert_ef @ 295'

      if(nodivert.eq.0)return   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----

      deallocate(colHeader%tb0cmd%colName)
      deallocate(colHeader%tb0cmd%colLocX)
      deallocate(colHeader%tb0cmd%colLocY)
      deallocate(colHeader%tb0cmd%colLocX1)
      deallocate(colHeader%tb0cmd%colLocY1)

d      if(iopt.eq.2)print*,'in read_divert_ef passed 274'


      if(nodivert.gt.0)then
!       FIND THE LOCAL COORDINATES FOR THE RESERVOIRS
!       THE VALUES FOR idiv AND jdiv ARE CHANGED

        do i=1,nodivert
!         convert to local coordinate unit system for new .res file
          jtake(i)=int((xtake(i)-xorigin)/xdelta)+1
          itake(i)=int((ytake(i)-yorigin)/ydelta)+1
          jgive(i)=int((xgive(i)-xorigin)/xdelta)+1
          igive(i)=int((ygive(i)-yorigin)/ydelta)+1
        end do
        if(iopt.ge.1)then
          write(53,*)
          write(53,*)'In read_divert_ef:'
          write(53,*)'Note: set iopt = 0 to not write this.'
          write(53,1011)
          write(53,1013)(i,itake(i),jtake(i),
     *           igive(i),jgive(i),divertname(i),i=1,nodivert)
        endif

!       THE ORDER OF READING THE COORDINATES OF THE RESERVOIRS
!       MUST BE THE SAME AS READING THE CORRESPONDING FLOWS
!       IN S/R REROUT.
!       READ RELEASES
!       THE RESERVOIR OUTFLOWS ARE ALL READ IN THE FIRST TIME
!       REROUT IS CALLED. THEY ARE THEN STORED AND USED EACH TIME
!       REROUT IS CALLED.
!       IF NATURAL CONTROL, FLOWS ARE SET TO -1.0 ABOVE

!       initialize releases
        do k=1,nodivert
          do j=1,ndiv
            qdivert(k,j)=-1.0
          end do
        end do
        
c!     rev. 9.7.24  Apr.  20/11  - NK: Added diverflg to indicate if a diversion is in grid
c!       set flags
c!       diverflg is initialized in read_shed_ef
c        do i=1,nodivert
c          ng=s(igive(i),jgive(i))
c          nt=s(itake(i),jtake(i))   
c          if(ng.eq.0.or.nt.eq.0)then     ! added Jan. 10/11 nk
c            if(ng.eq.0)then
c              print*,'WARNING:'
c              print*,'Diversion #',i,'receiving grid'
c              print*,'is not in the watershed.'
c              print*,'row=',igive(i),' column=',jgive(i)
cc              pause 'hit enter to continue & check others'
c            endif
c            if(nt.eq.0)then
c              print*,'WARNING:'
c              print*,'Diversion #',i,'source grid'
c              print*,'is not in the watershed.'
c              print*,'row=',itake(i),' column=',jtake(i)
cc              pause 'hit enter to continue & check others'
c            endif
c          else    
c            diverflg(ng)=.true.
c          endif
c        end do
            
c        if(b1(1).eq.0.0)then
        if(ndiv.gt.0)then
!         case with reservoir releases
!         do j=ktr,ndiv,ktr

!     rev. 9.1.13  Mar.  23/02  - fixed resv. timing, moved to beginning of dt
          do j=ktr,ndiv,ktr

c	print*,ktr,ndiv,ktr

            read(unitNum,*,iostat=ios)(qdivert(k,j),k=1,nodivert)
d            if(iopt.eq.2)print*,j,(qdivert(k,j),k=1,nodivert)
              if(ios.ne.0)then
                write(98,*)' Error on unit=',unitNum,' fln=',
     *                                     fln(flnNum)
                write(98,*)' Trying to read diversions hour =',j
                print*,' Error on unit=',unitNum,' fln=',fln(flnNum)
                print*,' Trying to read diversions'
                print*,' ios= ',ios
                if(ios.eq.-1)then
                  write(98,*)'End of file in fln= ',fln(flnNum)
                  write(98,*)
     *                 'Possibly last line does not have a return'
                  print*,'End of file in fln= ',fln(flnNum)
                  print*,'Possibly last line does not have a return'
                  print*
                else
                  print*
                  STOP ' program aborted in read_divert_ef.for'
                endif
              endif
!     rev. 9.5.24  Mar.  18/08  - NK: fixed missing data in read_resl_ef.f
!             fill in the gaps to hourly data
              if(ktr.gt.1)then
                do k=1,nodivert
	              do jj=ktr-1,1,-1
                      qdivert(k,j-jj)=qdivert(k,j)
c	print*,j,jj,j-jj,qdivert(k,j-jj)
	              end do
                end do
              endif

!             fill in missing data (-ve data)
              do k=1,nodivert
!     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_divert_ef
                if(qdivert(k,ktr).lt.0.0)then
                  print*
                  print*,'WARNING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                  print*,'first diversion flow should not be < 0.0'
                  print*,'diversions set to 0.0 until +ve value found'
                  qdivert(k,ktr)=0.0
                endif
!     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_divert_ef
                if(qdivert(k,j).lt.0.0)then
                  do jj=j,j+ktr-1
                    qdivert(k,jj)=qdivert(k,jj-1)
                  end do
                endif
              end do
            end do
          endif     !  if(ndiv.gt.0)


      endif !         if(nodivert.gt.0)

d     if(iopt.eq.2)print*,'in read_divert_ef passed 187'

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 
      istep=al/1000
c     if(iopt.eq.2)print*,'in read_divert_ef passed 272'

      if(nodivert.gt.0)then                             !999999999999999
!       FIND THE LOCAL COORDINATES FOR THE RESERVOIRS
!       THE VALUES FOR idiv AND jdiv ARE CHANGED

!       THE ORDER OF READING THE COORDINATES OF THE RESERVOIRS
!       MUST BE THE SAME AS READING THE CORRESPONDING FLOWS
!       IN S/R rdresv.

!       READ RELEASES
!       THE RESERVOIR OUTFLOWS ARE ALL READ IN THE FIRST TIME
!       rdresv IS CALLED. THEY ARE THEN STORED AND USED EACH TIME
!       rdresv IS CALLED.

!       IF NATURAL CONTROL, FLOWS ARE SET TO -1.0 ABOVE
d       if(iopt.eq.2)print*,'in read_divert_ef passed 304'

      endif           !   if(nodivert.gt.0)

d     if(iopt.eq.2)print*,'in read_divert_ef passed 551'

!     WE CAN'T HAVE -VE flows WHEN WE START
      do k=1,nodivert
        if(qdivert(k,1).lt.0.0)qdivert(k,ktr)=0.001
      end do

!     SET FUTURE RELEASES = LAST NON-NEGATIVE RELEASE
!     REALLY, WE'RE WORKING IN HOURLY INTERVALS ALTHOUGH THE      
!     RELEASES MAY BE READ IN ONLY WHEN THE RES OUTFLOW IS CHANGED.

!     rev. 9.5.14  Feb.  26/08  - NK: padded rel file for missing data
!     fill in missing data if rel file is shorter than the str file
!     nl is the length in nrs of the str file
!     ndiv is the length of the rel file
!     added Feb. 26/08  -nk-

      if(ndiv.gt.0)then
!       fill data at end of file only if there are values
        if(ndiv.lt.nl)then
          do j=ndiv+1,nl
            do k=1,nodivert
              qdivert(k,j)=qdivert(k,j-1)
            end do
          end do
          
          print*
          print*,'WARNING:'
          print*,'reservoir release file is shorter than the str file'
          print*,'missing data has been set = to last recorder release'
          print*
          do j=2,nl
            do k=1,nodivert
              if(qdivert(k,j).lt.0.0)qdivert(k,j)=qdivert(k,j-1)
            end do
          end do
        endif
      endif                        

      close(unit=unitNum,status='keep')

      firstpass='n'

  999 RETURN                 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-----

! FORMATS

  500 format(256f10.3)
 1011 format(' ',3x,'  i  idiv(i) jdiv(i)    b1(i)     b2(i)',
     *      '    b3(i)     b4(i)')
 1013 format(' ',3x,i3,4i8,a12/)
 6801 format('   read_divert_ef: reservoir no =',i3,' mhtot =',i5)
 6802 format('   ',i5,256f8.2)

      END SUBROUTINE read_divert

