      SUBROUTINE lst(scale,igrdshft)

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

!  THIS SUBROUTINE PRINTS THE FINAL RESULTS AND CALCULATES THE ERROR

!     REV. 7.80   Oct.  29/96 -  spl7 added yymmdd.rin for res inflows
!                             -  unit = 39   fln = 09
!     REV. 8.61 - Dec.  12/97 -  added contflg for statistics cont'n
!     REV. 8.77 - June   1/98 -  added sub-basin error calculation
!     REV. 8.83 - sep.  23/98 -  added step to the lst argument list
!     REV. 8.84 - Sep.  28/98 -  added runoff and evap fields to spl.txt
!     REV. 8.95 - Mar.  15/99 -  computed mean flows for time increment

!     REV. 9.00 - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90
!     rev. 9.1.18  Jun.  03/02  - Added sub-watershed modelling capability
!     rev. 9.1.22  Jul.  22/02  - Added simout\error.r2s file for ENSIM_Hydrologic
!     rev. 9.1.29  Nov.  07/02  - Changed the threshold flow values for error calculations
!     rev. 9.1.41  May   15/03  - Event average flows output to unit=75
!     rev. 9.1.48  Dec.  08/03  - NK: sumrechrge() added to get total recharge
!     rev. 9.1.67  Oct.  21/04  - NK; added unit 80 for lake_stor & lake_outflow
!     rev. 9.1.79  Mar.  30/05  - NK: ktri to area2 for reservoir inflow dt
!     rev. 9.2.19  Oct.  28/05  - NK: Compute daily & monthly flows
!     rev. 9.2.41  Jun.  15/06  - NK: changed the resin.txt file to resin.csv
!     rev. 9.3.07  Dec.  29/06  - NK: added sum_precip for whole domain
!     rev. 9.5.37  Oct.  14/08  - NK: added deltat_report to lake_sd.csv file write
!     rev. 9.5.56  Mar.  26/09  - NK: Fix bug with month in yearly events
!     rev. 9.5.60  Sep.  01/09  - NK: added deltat_report for lake_sd.csv file
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.65  Sep.  23/09  - NK: change class frac to whole basin values
!     rev. 9.5.78  Nov.  04/09  - NK: matched resvin locations to reach numbers
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
!     rev. 9.7.09  Oct.  11/10  - NK: update flowflag in lst.f for subsequent events
!     rev. 9.7.25  Apr.  28/11  - NK: Fixed daily flows
!     rev. 9.8.27  Sep.  27/12  - NK: changed action on resumflg='s' - keep tbcflg='y'
!     rev. 9.8.43  Jan.  31/13  - NK: fixed bug in lst.f : undefined output for iopt=99
!     rev. 9.8.48  Feb.  12/13  - NK: Replaced spl.plt with spl.tb0 file
!     rev. 9.8.51  Mar.  11/13  - NK: Link skiphours in s/r stats to value1 in the str file
!     rev. 9.8.88  Oct.  26/13  - NK: Fixed header writing sequence for spl.tb0 
!     rev. 10.1.50 Nov.  08/16  - NK: Overhauled lst for new isotope output
!     rev. 10.2.41 Dec.  10/18  - NK: Added winter monthly peaks

!  a - mm  runoff recorded
!  b - mm of runoff computed

!***********************************************************************
 
      use area_watflood
      use areacg
	implicit none


!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(10) :: time
      character(60) :: col_name(999)
      CHARACTER(8)  :: cday,tempflg
      character(80) :: line
      character(1)   :: dsflg,firstpass,answer
      logical       :: exists,tmpflg,errflg1
      
!     the MRBflg is custom for revising the MRB_WILLISTON_INFLOWS>tb0 
!     file for the MRB model. The natural flows will be added to col 3 
!     CHARM will look for the file and if present, add the data is so directed    
      logical       :: MRBflg
      character(72) :: MRBheader(100)
      integer       :: MRBheaderlines,MRBdatalines
      integer       :: day_last,jd,jul_day1_id1
      real(4),  dimension(:,:), allocatable ::  MRBdata
 
      real(4) :: aaa,bbb,ccc,ab,scale,tttmmmpp,qmaxs,
     *            qsynmean,smc5(99),             !qdwprmd(60),
     *            diff2,qmaxh,aig,ajg,
     *            qhydevent(1000),qsynevent(1000),
     *            qdwrp_sum
      integer :: i,j,n,ii,ig,jg,ios,n1,kk,l,noo,k,nend,evt_hrs,
     *             nmissed,igrdshft,kfirst,mh,nqevent(1000),mon,
     *             mohours(24),nmonths,month_1,hour_no,hour_neg,
     *             time_sum,iDeallocate,iAllocate,level_count,
     *             strLength,un

      DATA mohours/744,672,744,720,744,720,744,744,720,744,720,744,
     *             744,672,744,720,744,720,744,744,720,744,720,744/
     
      real(4),  dimension(:),   allocatable :: lake_elv_max
      real(4),  dimension(:),   allocatable :: lake_elv_min
      real(4),  dimension(:),   allocatable :: lake_elv_max_av
      real(4),  dimension(:),   allocatable :: lake_elv_min_av
      real(4),  dimension(:),   allocatable :: sum_max
      real(4),  dimension(:),   allocatable :: sum_min
      real(4),  dimension(:),   allocatable :: mnl_diff
      real(4),  dimension(:,:), allocatable :: lake_elv_max_id
      real(4),  dimension(:,:), allocatable :: lake_elv_min_id
      real(4),  dimension(:),   allocatable :: qdwpr_sum

      data time_sum/0/
      data firstpass/'y'/

!       SET UP A FEW FILENAMES FOR BACKGROUND STREAMFLOW FILES
!       FIVE IS ALL WE NEED FOR BC HYDRO

!     THIS SECTION GETS THE COMPUTER'S DATE AND TIME
!     THIS INFO USED TO COME FROM HEADER.FI, HOWEVER F90 INTRINSICS
!     AND CALLS ARE DIFFERENT AND THEREFORE IT NEEDED TO BE MODFIED.
      call date_and_time(cday,time)

!       MAX FLOW CALCUL'D IN EACH SQUARE IS PRINTED WHEN IOPT.GT.0
!       GIVES QMAX=FN(DRAINAGE AREA) FOR IOPT.NEQ.0 ONLY

!!!!!!!!!!!!!!!!!!!!!!!      if(ni.le.1)then

!     CALUCLATE THE TOTAL PRECIP ON EACH SQUARE SUB-BASIN
!     THERE IS A PROBLEM IN THAT WE ARE CALCULATING ON EACH SUB-AREA, 
!     NOT THE TOTAL CONSTRIBUTING AREA

d	if(iopt.eq.2)print*,' checkpoint 10 in lst'

      if(firstpass.eq.'y')then
        allocate(mnl_diff(no),stat=iAllocate)
        if(iAllocate.ne.0)STOP
     *         'Error with allocation of mnl_diff @ 108'
        allocate(net_basin_supply(noresv,100),stat=iAllocate)  ! good for up to 100 days
        if(iAllocate.ne.0)STOP
     *         'Error with allocation of net_basin_supply @ 144'
         
        if(iopt99)then
          allocate(qp_month_syn(no,13),
     *             qp_month_hyd(no,13),stat=iAllocate)
          if(iAllocate.ne.0)STOP
     *         'Error with allocation of qp_month** @ 127'
      
      endif
        
        
      endif

!     MRB model custom code
      if(firstpass.eq.'y')then
        INQUIRE(FILE='MRB_WILLISTON_INFLOWS.tb0',EXIST=exists)
        if(exists)then
          MRBflg=.true.
        else
          MRBflg=.false. 
        endif  
        if(MRBflg)then
          print*
          print*,'Found the file: MRB_WILLISTON_INFLOWS.tb0'
          print*,'MRBflg =',MRBflg
          print*
          print*,'Would you like to replace col 3 with CHARM generated'
          print*,'natural flows at the Bennett dam?'
          print*,'To make this work, you will need to run'
          print*,'1960_natural.evt  as the event.evt file'
          print*,'and answer `y` next, else answer `n`'
          print*
          print*,'If you do not wish to see this option with each run'
          print*,'answer `s` to stop this run'
          print*,'Then move the MRB_WILLISTON_INFLOWS.tb0 file from'
          print*,'the working directory to somewhere safe.' 
          print*
          read*,answer
          if(answer.eq.'s')stop
          if(answer.eq.'n')MRBflg=.false.
c          if(answer.eq.'y')then
c            print*,'What is the rank for the Hudson Hope WSC'
c            print*,'Flow station?'
c            read*,HHrank
c          endif
        endif
        if(MRBflg)then
          open(unit=99,file='MRB_WILLISTON_INFLOWS.tb0',
     *                          status='old',iostat=ios)
          read(99,99501)line
99501     format(a80)  
          n=1        
          MRBheader(n)=line(1:72)
          do while(.not.line(1:10).eq.':EndHeader')
            n=n+1
            read(99,99501)line
            MRBheader(n)=line(1:72)
          end do
          MRBheaderlines=n
          do i=1,n
            print*,MRBheader(i)(1:72)
          end do
          n=0
          do while(.not.eof(99))
            n=n+1
            read(99,99501)line
          end do
          print*,'found ',n,' lines with data'
          MRBdatalines=n
          allocate(MRBdata(3,MRBdatalines),stat=iAllocate)
          if(iAllocate.ne.0)STOP
     *         'Error with allocation of gname lst & 522'
          rewind(unit=99)
          do n=1,MRBheaderlines
            read(99,*)line
          end do
          do n=1,MRBdatalines
            read(99,*,iostat=ios)(MRBdata(i,n),i=1,3)
            if(ios.ne.0)then
              print*,'Error in data line ',n
              stop 'Program aborted in lst @ 188'
            endif
            write(555,*)(MRBdata(i,n),i=1,3)
          end do
          close(unit=99)
c        endif     ! ** moved down
        
!         write the header for the file updated with CHARM computed natural flows        
          open(unit=99,file='MRB_WILLISTON_INFLOWS.new',
     *                       status='unknown',iostat=ios)
          if(ios.ne.0)then
            print*,'Problems opening MRB_WILLISTON_INFLOWS.new'
            stop 'CHARM aborted in lst @ 193'
          endif
          MRBheader(8)=':WrittenBy                WATFLOOD/CHARM'
          do n=1,MRBheaderlines
            if(n.eq.9)then
              write(99,3010)':CreationDate       ',
     *        cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
            else
               Write(99,99502)MRBheader(n)
99502         format(a72)
            endif
          end do
          close(unit=99,status='keep')
        endif    !  ** moved from above  Oct. 18/17 NK
        

!     REV. 10.1.36 Jul   12/16  - NK: Added results\LakeName.tb0
!       Write the headers for the yyyymmdd_lake.tb0 files
!       make up file names 
        if(noresv.gt.0.and.tb0flg.eq.'y')then
        do l=1,noresv 
          strLength = LEN_TRIM(resname(l))
          print*,'strLength=',strLength
                 
          if(l.le.9)then
             write(line,50001)
     *        'results\lake_00',l,'_',resname(l),'.tb0'
             write(*,50001)
     *        'results\lake_00',l,'_',resname(l),'.tb0'
50001        format(a15,i1,a1,a<strLength>,a4)
          elseif(l.le.99)then  
             write(line,50002)
     *        'results\lake_0',l,'_',resname(l),'.tb0'
             write(*,50002)
     *        'results\lake_0',l,'_',resname(l),'.tb0'
50002        format(a14,i2,a1,a<strLength>,a4)
          else
             write(line,50003)
     *        'results\lake_',l,'_',resname(l),'.tb0'
             write(*,50003)
     *        'results\lake_',l,'_',resname(l),'.tb0'
50003        format(a13,i3,a1,a<strLength>,a4)  
          endif
          read(line,50004)line
50004     format(a72)     
          print*,line(1:72)        
          un=2000+l
          open(unit=un,file=line,status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',line(1:40)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in write_tb0.f @ 52'
          endif
          PRINT*
	    if(iopt.ge.1)print*,'Opened unit=',un,' filename=',line(1:40)
          write(un,3005,iostat=ios)
     *              '##########################################'
          if(ios.ne.0)then
            print*,'Unable to write to ',line(1:40)
            print*,'Likely problem: duplicate names in the rel file'
            print*
            stop 'Program aborted in lst @ 252'
          endif
          write(un,3005)':FileType tb0  ASCII  EnSim 1.0           '
          write(un,3005)'#                                         '
	    write(un,3005)'# DataType              EnSim Table Data  '
          write(un,3005)'#                                         '
          write(un,3005)':Application            WATFLOOD          '
	    write(un,3005)':Version                2.1.23            '
	    write(un,3005)':WrittenBy              WATFLOOD/CHARM    '
          call date_and_time(cday,time)
	    write(un,3010)':CreationDate       ',
     *   cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
3010      format(a20,a4,'-',a2,'-',a2,2x,a2,':',a2)
          write(un,3005)'#                                         '
	    write(un,3005)'#-----------------------------------------'
          write(un,3005)':SourceFile             Various           '                                    '
          write(un,3005)'#                                         '
          write(un,3005)'#                                         '
!         fixed format for the date      
          if(mo_start.le.9.and.hour_start.le.9)then
            write(un,3024)':StartDate          ',
     *   year_start,mo_start,day_start,hour_start
3024        format(a20,i4,'/0',i1,'/0',i1,'  ',i2,':00:00')  
          elseif(mo_start.le.9.and.hour_start.gt.9)then
            write(un,3025)':StartDate          ',
     *     year_start,mo_start,day_start,hour_start
3025        format(a20,i4,'/0',i1,'/',i2,'  ',i2,':00:00')  
          elseif(mo_start.gt.9.and.hour_start.le.9)then
            write(un,3026)':StartDate          ',
     *     year_start,mo_start,day_start,hour_start
3026        format(a20,i4,'/',i2,'/0',i1,'  ',i2,':00:00')  
          else    
            write(un,3027)':StartDate          ',
     *         year_start,mo_start,day_start,hour_start
          endif
3027      format(a20,i4,'/',i2,'/',i2,'  ',i2,':00:00')      
          write(un,3005)':deltaT             24:00:00.000          '
!         This next line asked foe by AB @ NRC
!         but there should never be missing data in these computed 
!         lake elv. inflow & outflow          
!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
          write(un,3005)':NoDataValue -999.0                       '
3028      format(a20,i10)   
          write(un,3005)'#                                         '
	    write(un,3005)':ColumnMetaData                           '
	    write(un,3006)
     *	    '   :ColumnName       level      inflow     outflow'
          write(un,3006)
     *      '   :ColumnUnits       masl         cms         cms'
	    write(un,3006)
     *	    '   :ColumnType       float       float       float'
          write(un,3005)':endHeader                                ' 
 3001     format(a20,i16)
 3005     format(a42)
 3006     format(a50)
 3020     format(a20,a40)
 3021     format(a20,999(a42))

        end do
        endif    ! noresv > 0
!     END    REV. 10.1.36 Jul   12/16  - NK: Added results\LakeName.tb0

!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
        if(nbsflg.eq.'0')then
!         Open the file and write header in the first event only
          open(unit=957,file=filename(957),status='unknown') 
          write(957,*)'Sub-basin supply (between numbered lakes)'
	    write(957,3010)'Run time            ',
     *        cday(1:4),cday(5:6),cday(7:8),time(1:2),time(3:4)
	    write(957,3010)'Valid approx.       ',
     *        cday(1:4),cday(5:6),cday(7:8),00,00
          write(957,95701)(l,l=1,noresv)
95701     format(<noresv>i8)
          write(957,95702)(resname(l),l=1,noresv)
95702     format(<noresv>(1x,a7))
        endif
        
      endif   !first pass
      
c      stop 'in lst'
      
      do i=1,ij
         do j=1,ji
            p(i,j)=0.0
!           This used to be -9999. but it gave format errors in 
!           the newerror.txt file. This may have an effect
!            when the file is read.   nk  Apr. 10/03
            basinerr(i,j)=-999.
         end do
      end do
d	if(iopt.eq.2)print*,' checkpoint 11 in lst'

!     P(I,J) IS THE CUMMULATIVE RAINFALL NOW
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         p(i,j)=sump(n)
      end do
d	if(iopt.eq.2)print*,' checkpoint 12 in lst'

!     THIS NEEDS TO BE DONE EACH EVENT BECAUSE SUMP() AND P() ARE
!     ALREADY THE SUM OF THE PRECIP SINCE THE START      
      do l=1,no
        nnsum(l)=0
        ppsum(l)=0.0
        qhydevent(l)=0.0
        qsynevent(l)=0.0
        nqevent(l)=0
      end do
d	if(iopt.eq.2)print*,' checkpoint 13 in lst'
 
      do n=1,naa
         i=yyy(n)
         j=xxx(n)
         l=nhyd(i,j)
         if(l.gt.0)then
            ppsum(l)=ppsum(l)+p(i,j)
            nnsum(l)=nnsum(l)+1
         endif
      end do
d	if(iopt.eq.2)print*,' checkpoint 14 in lst'

      do l=1,no
         if(nnsum(l).gt.0)then
            ppsum(l)=ppsum(l)/float(nnsum(l))
         else
            ppsum(l)=-10.
         endif
      end do
d	if(iopt.eq.2)print*,' checkpoint 15 in lst'

      if(iopt.ge.1)then 	
         write(51,6080)
         write(51,6081)
         write(51,6666)(n,yyy(n),xxx(n),da(n),qmax(n),sump(n),n=1,naa)
      endif
d	if(iopt.eq.2)print*,' checkpoint 16 in lst'

!     P(I,J) IS USED HERE ONLY TO REPORT PRECIPITATION ON
!     A SQUARE GRID MATRIX, THIS P(I,J) IS NOT USED IN ANY
!     COMPUTATIONS, AND HAS DIFFERENT VALUE FROM THAT IN S/R,
!     RAIN,RUNOFF,ETC...
      write(51,6668)scale
	write(51,6667)scalesnw

!     rev. 9.9.28  Sep.  30/14  - NK: Remove unnecessary writes for watroute
      if(modelflg.eq.'n')then
!       WRITE THE PRECIP ON EACH GRID
        if(xcount.le.250) then
          do i=ycount,1,-1
            write(51,6090)(p(i,j),j=1,xcount)
          end do
	  else
	    do i=ycount,1,-1
            write(51,6092)(p(i,j),j=1,xcount)
          end do
	  endif

        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=sumrff(n)
        end do
        write(51,6672)
        if(xcount.le.250) then
          do i=ycount,1,-1
            write(51,6090)(p(i,j),j=1,xcount)
          end do
	  else
	    do i=ycount,1,-1
            write(51,6092)(p(i,j),j=1,xcount)
          end do
	  endif
!       P(I,J) IS THE CUMMULATIVE RUNOFF NOW

        do n=1,naa
          i=yyy(n)
          j=xxx(n)
!         TRY PRINTING OUT ACTUAL CUMMULATIVE EVAPORATION FRANK S FEB/99
          p(i,j)=eloss(n)
        end do
        write(51,6673)
        if(xcount.le.250) then
          do i=ycount,1,-1
            write(51,6090)(p(i,j),j=1,xcount)
          end do
	  else
	    do i=ycount,1,-1
            write(51,6092)(p(i,j),j=1,xcount)
          end do
	  endif
!       P(I,J) IS THE TOTAL EVAPORATION 
!            +/- CHANGE IN BASIN ELEMENT STORAGE

!       rev. 9.1.48  Dec.  08/03  - NK: sumrechrge() added to get total recharge
!       write the recharge totals:
!       Write THE TOTAL recharge to spl.txt
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          p(i,j)=sumrechrg(n)/frac(n)
        end do
        write(51,6677)

        if(xcount.le.250) then
          do i=ycount,1,-1
c           write(51,6090)(p(i,j),j=1,xcount)
           write(51,*)(p(i,j),j=1,xcount)
          end do
	  else
	    do i=ycount,1,-1
c            write(51,6092)(p(i,j),j=1,xcount)
            write(51,*)(p(i,j),j=1,xcount)
         end do
	  endif
!	
!       Write the runoff coefficient to spl.txt 
        do n=1,naa
           i=yyy(n)
           j=xxx(n)
           if(sump(n).gt.0.01)then
              p(i,j)=sumrff(n)/sump(n)*100.
           else
              p(i,j)=-1.0
           endif
        end do
        write(51,6674)
        if(xcount.le.250) then
          do i=ycount,1,-1
            write(51,6090)(p(i,j),j=1,xcount)
          end do
	  else
	    do i=ycount,1,-1
            write(51,6092)(p(i,j),j=1,xcount)
          end do
	  endif

!       WRITE A FILE TO PLOT PRECIP, LOSS AND RUNOFF FIELDS
        if(iopt99)then
          open(unit=99,file='debug\lossplot.txt',
     *                   status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file lossplot.txt'
            print*,'Possible causes:'
            print*,'file in use by another application'
            stop 'Program aborted in lst.f @ 312'
          endif
          write(99,6675)
          do n=1,naa
            i=yyy(n)
            j=xxx(n)
            aig=yorigin+(i-1)*ydelta+ydelta/2.0
            ajg=xorigin+(j-1)*xdelta+xdelta/2.0
            if(sump(n).gt.0.001)then
              tttmmmpp=sumrff(n)/sump(n)*100.0
            else
              tttmmmpp=-1.0
            endif
            write(99,6094)nbasin(i,j),i,j,aig,ajg,sump(n),sumrff(n),
     *                 sump(n)-sumrff(n),tttmmmpp
          end do
          close(unit=99,status='keep')
        endif
      endif    ! modelflg='n'
!     end rev. 9.9.28  Sep.  30/14  - NK: Remove unnecessary writes for watroute


!     JUST TO MAKE SURE WE DON'T GET HUGE AMOUNTS REMAINING IN MEMORY
!     FOR THE NEXT MONTH:
      do i=1,ycount
         do j=1,xcount
            p(i,j)=0.0
         end do
      end do

!!!!!!!!!!!!!!!!!!!!!!         endif

!     THIS SECTION PRINTS A SUMMARY OF THE MEASURED AND COMPUTED
!     FLOWS FOR PLOTTING
d	if(iopt.eq.2)print*,' checkpoint 30 in lst'

      open(unit=27,file=filename(27),status='unknown',
     *     iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file',filename(27)(1:40)
        print*,'Possible causes:'
        print*,'file in use by another application'
        print*,'or target directory does not exist'
        stop 'Program aborted in lst.f @ 312'
      endif

      write(51,6015)
     *  time(1:2),time(3:4),time(5:6),cday(1:4),cday(5:6),cday(7:8)
      if(writeFlg(1))write(27,6015)
     *  time(1:2),time(3:4),time(5:6),cday(1:4),cday(5:6),cday(7:8)
c      write(*,6015)
c     *  time(1:2),time(3:4),time(5:6),cday(1:4),cday(5:6),cday(7:8)

!      write(51,6015)hrs,mins,secs,day,month,year
!      write(99,6015)hrs,mins,secs,day,month,year

      write(51,1013)
      if(writeFlg(1))write(27,1013)

d            if(iopt.eq.2)then
d              print*,' In lst before writing to unit 58'
d              print*
d            endif

d	if(iopt.eq.2)print*,' checkpoint 40 in lst'

!     * * * * * * * * * * * * * * * * 
!     Calculate the statistics
!     * * * * * * * * * * * * * * * * 

c      do n=1,no,9
c         n1=n+8
c         n1=min0(n1,no)
c         noo=n1-n+1

c         write(58,1000)fln(10)
c         write(58,1015)noo,nl,mhtot,kt
c         write(58,1014)

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      do l =1,no

           if(inbsnflg(l).eq.1.and.area(l).gt.0.0)then  ! area added Feb. 11/11 nk

!           SEE S/R SUB FOR INTIALIZATION OF NQ,NQC,AA,BB,CC,ASHNUM,
!           ASHDEN,QPEAKH & QPEAKS
            if(id.eq.1)then
               qmaxh=0.1e-10
               qmaxs=0.0e-10
               qbarobs(l)=-1.0
            endif

!           RESET THE DAILY PEAKS FOR EACH EVENT FOR THE DAILY PEAK 
!           TABLE
c            dpeakh(l)=0.1e-10
c            dpeaks(l)=0.1e-10
            qpeakh(l)=0.1e-10
            qpeaks(l)=0.1e-10
            volhyd(l)=0.0
            volsyn(l)=0.0
	  
!   REV. 8.95 - Mar.  15/99 -  COMPUTED MEAN FLOWS FOR TIME INCREMENT
!           CALCULATE THE MEAN FLOWS FOR THE TIME INCREMENT:
!           QHYD(L,K) WILL BE MODIFIED TO BECOME THE MEAN DAILY FLOW
!           INSTEAD OF THE INSTANTANEOUS FLOW (BUT ONLY FOR DAILY'S)

            if(kt.eq.24)then
               do k=kt,mhtot,kt
                  qsynmean=0.0
                  do kk=1,kt
                     qsynmean=qsynmean+qsyn(l,k-kk+1)  
                  end do
                  qsyn(l,k)=qsynmean/float(kt)         ! mean daily <<
c                  dpeakh(l)=amax1(dpeakh(l),qhyd(l,k))
c                  dpeaks(l)=amax1(dpeaks(l),qsyn(l,k))
                  qpeakh(l)=amax1(qpeakh(l),qhyd(l,k))
                  qpeaks(l)=amax1(qpeaks(l),qsyn(l,k))
              end do
            else
               do k=kt,mhtot,kt
c                  dpeakh(l)=amax1(dpeakh(l),qhyd(l,k))
c                  dpeaks(l)=amax1(dpeaks(l),qsyn(l,k))
                  qpeakh(l)=amax1(qpeakh(l),qhyd(l,k))
                  qpeaks(l)=amax1(qpeaks(l),qsyn(l,k))
               end do
            endif

!     rev. 10.2.15 Feb.  05/18  - NK: Added 'results\monthly_peaks'
!           pick our the monthly peak flow pairs (Obs/Comp)
!           qp_month_sym(no,12),qp_month_hyd(no,12),   
!           in first day of the event set peaks to 0.001
!           (Needed for log plots)
            if(iopt99)then
              do j=1,13
                qp_month_hyd(l,j)=0.001
                qp_month_syn(l,j)=0.001
              end do
              do k=1,mhtot
                j=k/732+1
                if(qhyd(l,k).gt.qp_month_hyd(l,j))then
                  qp_month_hyd(l,j)=qhyd(l,k)
!                 pick the max flow in a 12 hour window
                  qp_month_syn(l,j)=qsyn(l,k)
                  if(k.gt.25.and.k.lt.mhtot-25)then
                      do i=k-24,k+24
                        qp_month_syn(l,j)=
     *                   amax1(qp_month_syn(l,j),qsyn(l,i))
                      end do
c                      qp_month_syn(l,j)=amax1(qp_month_syn(l,j),
c     *                qsyn(l,k-5),qsyn(l,k-4),qsyn(l,k-3),
c     *                qsyn(l,k-2),qsyn(l,k-1),                  
c     *                qsyn(l,k+5),qsyn(l,k+4),qsyn(l,k+3),
c     *                qsyn(l,k+2),qsyn(l,k+1))     
                  endif
                endif
c                write(450+l,45300)k,j,qhyd(l,k),qsyn(l,k),
c     *                  qp_month_hyd(l,j),qp_month_syn(l,j)
c45300           format(2i5,4f10.1)              
              end do
            endif

!     rev. 9.8.51  Mar.  11/13  - NK: Link skiphours in s/r stats to value1 in the str file
            if(.not.skipflg)then
              do k=kt,nl,kt
c               qhyd(l,k)=amax1(qhyd(l,k),0.100E-32)
!               no good for the stats program

!               CALCULATE HYDROGRAPH VOLUMES:
!     rev. 10.2.64 Sep.  09/26  - NK added min_flow_cutoff for error calculations
                if(qhyd(l,k).gt.min_flow_cutoff)then 
                  nq(l)=nq(l)+1
                  nqevent(l)=nqevent(l)+1
                  aa(l)=aa(l)+qhyd(l,k)
                  bb(l)=bb(l)+qsyn(l,k)
                  volhyd(l)=volhyd(l)+qhyd(l,k)
                  volsyn(l)=volsyn(l)+qsyn(l,k)
                  qhydevent(l)=qhydevent(l)+qhyd(l,k)
	            qsynevent(l)=qsynevent(l)+qsyn(l,k)
                endif
                cc(l)=cc(l)+qsyn(l,k)
                nqc(l)=nqc(l)+1
              end do

!	        pre-emption
!             calculated station errors:ashnum(l)
              if(nq(l).gt.0)then
                qbarobs(l)=aa(l)/nq(l)
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.0000)then
                     ashnum(l)=ashnum(l)+(qhyd(l,k)-qsyn(l,k))**2
                     ashden(l)=ashden(l)+(qhyd(l,k)-qbarobs(l))**2
                  endif
                end do
                if(ashden(l).ne.0.0)then
                  rsquare(l)=1.0-ashnum(l)/ashden(l) ! actually Nash effeciency
                else
                  rsquare(l)=-99.
                endif
              else
                rsquare(l)=-99.
              endif
            endif    !skipflg

!   REV. 8.83 - Nov. 30/96 - FIX DIVISION BY 0 - CHECK
c            if(area(l).le.0.0)then
c               write(6,6020)id,l,area(l)
c               print*,'spl.csv not written'
c               STOP 'Probram aborted in lst @ 410'  
c            endif

!           calculate the mm runoff          
            aaa=amin1(99999.0,aa(l)*float(kt)*3.6/area(l))
            bbb=amin1(99999.0,bb(l)*float(kt)*3.6/area(l))
            ccc=amin1(99999.0,cc(l)*float(kt)*3.6/area(l))
            if(aaa.gt.99998.9)aaa=-999.
            if(bbb.gt.99998.9)bbb=-999.
            if(ccc.gt.99998.9)ccc=-999.
            volhyd(l)=volhyd(l)*float(kt)*3.6/area(l)
            volsyn(l)=volsyn(l)*float(kt)*3.6/area(l)
!           calculate Dv
c            if(aaa.gt.0.0)then
            if(aaa.gt.1.0)then
               ab=(bbb-aaa)/aaa*100.
               statnerr(l)=(bb(l)-aa(l))/aa(l)*100.0
            else
               ab=-1.0
               statnerr(l)=-999.9
            endif

!           reinitialize flowflag - in flowinit is is basin only on year 1
            flowflag(l)=.false.
            if(statnerr(l).gt.-999.0)flowflag(l)=.true.

!           CALCULATION & PRINTING OF RECORDED DATA & CALCULATED  
!           HYDROGRAPHS.
!           FLOWS ARE REPORTED AT KT INTERVALS.

!		  write precip.txt
!     rev. 9.8.43  Jan.  31/13  - NK: fixed bug in lst.f : undefined output for iopt=99
!     rev. 9.8.45  Jan.  31/13  - NK: disabled some writes for iopt = 99

!     REV. 10.1.31 May   15/16  - NK: Revised output to precip.txt : include all str stations
!           section moved outside inbsnflg loop (below)

            qmaxh=amax1(qmaxh,qpeakh(l))
            qmaxs=amax1(qmaxs,qpeaks(l))

           endif   !  inbsnflg
           
           
!     REV. 10.1.31 May   15/16  - NK: Revised output to precip.txt : include all str stations
c            if(iopt.lt.99.and.iopt.gt.0)then
            if(writeFlg(1))then   ! always print precip.txt so results can be compared during dds
!     rev. 9.5.65  Sep.  23/09  - NK: change class frac to whole basin values
!             write only where there is data
	        i=(ystr(l)-yorigin)/ydelta+1
	        j=(xstr(l)-xorigin)/xdelta+1
              if(aaa.gt.ppsum(l)*0.01.and.inbsnflg(l).eq.1)then
                if(iopt.ge.2.and.iopt.le.10)
     *                    print*,'i,j,aaa,ppsum',i,j,aaa,ppsum
	          if(sta_area(l).gt.1.0)then
                  write(27,1018,iostat=ios)
c                  write(27,*)
     *            l,gage(l),area(l),ppsum(l),aaa,bbb,ccc,statnerr(l),
     *            rsquare(l),qpeakh(l),qpeaks(l),sta_area(l),area(l),
     *            (area(l)-sta_area(l))/sta_area(l)*100.0,
     *          ((areaclass(s(i,j),ii)/areasum(s(i,j))),ii=1,classcount)
	          else
                  if(iopt.ge.2.and.iopt.le.10)
     *                         print*,'i,j',i,j,areasum(s(i,j))
                  if(iopt.ge.2.and.iopt.le.10)print*,'i,j',i,j,
     *                 (areaclass(s(i,j),ii),ii=1,classcount)
c                  write(27,*)
                  write(27,1018,iostat=ios)
     *            l,gage(l),area(l),ppsum(l),aaa,bbb,ccc,statnerr(l),
     *            rsquare(l),qpeakh(l),qpeaks(l),sta_area(l),area(l),
     *                 -9999.0,
     *                 ((areaclass(s(i,j),ii)/areasum(s(i,j))),
     *                 ii=1,classcount)
                  endif
              else
!     rev. 10.1.09 Dec.  07/15  - NK: Add blank line for missing data in the precip.txt file in lst.f 
!               write a blank line to not screw up the numbering
!     rev. 10.2.18 Mar.  12/18  - NK: Fixed array fault in read_resv_ef and sub
                if(i.gt.0.and.i.le.ycount.and.
     *                             j.gt.0.and.j.le.xcount)then  
c                  if(s(i,j).gt.0)then
c                    write(27,1019,iostat=ios)l,gage(l),area(l),
                    write(27,1018,iostat=ios)l,gage(l),-999.0,
     *              -999.0,-999.0,-999.0,-999.0,-999.0,-99.0,
     *              -999.0,-999.0,-999.0,-999.0,-999.0,
     *               (-1.0,ii=1,classcount)
c                  endif
                endif
!     rev. 10.1.73 Mar.  27/17  - NK: Advisory message set in precip.txt for iopt=0 
              endif
            else  
              if(writeFlg(1))write(27,*)
              if(writeFlg(1))write(27,*)'Set iopt > 0 to get this table'
              if(writeFlg(1))write(27,*)
            endif

      end do
!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

d	   if(iopt.eq.2)print*,' checkpoint 50 in lst'

!!!!     id=1

c!        WRITE THE SPL.PLT FILES - FOR THE SPLPLT PROGRAM
cc         write(58,6011)id,(smc5(ii),ii=1,classcount)
c         write(58,6400)
c
c         if(qmaxh.lt.10.0.and.qmaxs.lt.10.0)then
c            do k=kt,nl,kt
c               write(58,6401)k,(qhyd(l,k),qsyn(l,k),l=n,n1)
c            end do
c         else   
c            do k=kt,nl,kt
c               write(58,6402)k,(qhyd(l,k),qsyn(l,k),l=n,n1)
c            end do
c         endif
c       end do

c      close(unit=99,status='keep')


c      SUBROUTINE write_tb0(un,fn,nfg,ng,no_signf)

c     write the spl.tb0 file
!     rev. 9.8.48  Feb.  12/13  - NK: Replaced spl.plt with spl.tb0 file
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
!                               - had to add xsta1 & ysta1 to allocation here as well

c      if(id.eq.1)then
      if(allocated(gname))then
c        deallocate(gname,xsta,ysta,xsta1,ysta1,stat=iDeallocate)
        deallocate(gname,stat=iDeallocate)
        if(ideallocate.ne.0)print*,
     *          'problem deallocating gname in lst @ 506'

      endif
      if(allocated(xsta))then
        deallocate(xsta,ysta,stat=iDeallocate)
        if(ideallocate.ne.0)print*,
     *          'problem deallocating xsta in lst @ 507'

      endif
      if(allocated(xsta1))then
        deallocate(xsta1,ysta1,stat=iDeallocate)
        if(ideallocate.ne.0)print*,
     *          'problem deallocating xsta1 in lst @ 508'

      endif
      
      allocate(gname(no*2),xsta(no*2),ysta(no*2),
     *              xsta1(no*2),ysta1(no*2),stat=iAllocate)
      if(iAllocate.ne.0)STOP
     *         'Error with allocation of gname lst & 522'


!     rev. 9.8.88  Oct.  26/13  - NK: Fixed header writing sequence for spl.tb0 
      author='WATFLOOD/CHARM                          '
      name='WATFLOOD model output: recorded/computed'
      coordsys_temp=coordsys1
      !     GreenKenue uses LatLong - code below uses LATLONG
      if(coordsys_temp.eq.'LatLong   ')coordsys_temp='LATLONG   '
      if(coordsys_temp.eq.'Cartesian ')coordsys_temp='CARTESIAN '
      zone_temp=zone1
      datum_temp=datum1
      xorigin_temp=xorigin
      yorigin_temp=yorigin
      xcount_temp=no*2
      ycount_temp=nl/kt
      startdate='unknown   '
      starttime='unknown   '
      unit_conversion=1.0
      deltat_temp=kt*3600
      source_file_name='last spl run'   
      do l=1,no
        gname(l*2-1)=gage(l)
        xsta(l*2-1)=xstr(l)
        ysta(l*2-1)=ystr(l)
        gname(l*2)=gage(l)
        xsta(l*2)=xstr(l)
        ysta(l*2)=ystr(l)
      end do
      do l=1,2*no
        column_units(l)='cms' 
        column_type(l)='float'
      end do
!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
      noDataValue=-999.0   
      
!     rev. 9.8.95  Nov.  20/13  - NK: Changed unit 58 to 22 for spl.tb0
      if(id.eq.1)then
!       write the header only once
!       and then just add data for each event      
!       writing to spl.tb0
        if(writeFlg(1))call write_tb0(955,955,1,no*2,8)
      endif

      if(allocated(outarray))then
        deallocate(outarray,stat=iDeallocate)
        if(ideallocate.ne.0)print*,
     *          'problem deallocating outarray in lst @ 573'
      endif
        
      allocate(outarray(nl/kt,no*2),stat=iAllocate)
      if(iAllocate.ne.0)STOP
     *         'Error with allocation of outarray lst & 511'
          
      do k=kt,nl,kt
        do l=1,no
          outarray(k/kt,l*2-1)=qhyd(l,k)
          outarray(k/kt,l*2)=qsyn(l,k)
        end do
      end do

!     rev. 9.8.95  Nov.  20/13  - NK: Changed unit 58 to 22 for spl.tb0
!     write the data to spl.tb0
      noDataValue=-999.0
      if(writeFlg(1))call write_tb0(955,955,0,no*2,12)    
        
!     close file after last event.      
      if(id.eq.ni)then
c        close(unit=22,status='keep')
        if(writeFlg(1))then
          write(51,*)'Closed unit 955 Filename=  ',fln(955)(1:40)
          write(*,*)'Closed unit 955 Filename=  ',fln(955)(1:40)
        endif
      endif

      deallocate(outarray,stat=iDeallocate)
      if(ideallocate.ne.0)print*,
     *          'problem deallocating outarray in lst @ 579'
      allocate(outarray(ycount,xcount),stat=iAllocate)
      if(iAllocate.ne.0)STOP
     *         'Error with allocation of outarray lst & 582'


       
!     write a new mean_observed_flows.txt file if .not.exists
      if(dds_flag.eq.0.and.id.eq.ni)then
        inquire(FILE='mean_observed_flows.txt',EXIST=exists)
        if(.not.exists)then
          do l=1,no+noresvi
            if(nq(l).le.0)qbarobs(l)=-1.0
          end do
          open(unit=99,file='mean_observed_flows.txt',
     *              status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file mean_observed_flows.txt'
            print*,'Possible causes:'
            print*,'file in use by another application'
            stop 'Program aborted in lst.f @ 657'
          endif
d         print*,'reading the mean flows in sub'
          write(99,*)'Station Averages - observed'
          do l=1,no+noresvi
            write(99,*)l,qbarobs(l)
d            write(*,*)i,mean_observed(l)
          end do
          close(unit=99,status='keep')
          print*
          print*,'mean_observed_flows.txt written'
          print*,'for ',no,' flow stations and'
          print*,'for ',noresvi,' reservoirs'
          print*
        endif
      endif

!     * * * * * * * * * * * * * * * * 

d	if(iopt.eq.2)print*,' checkpoint 60 in lst'

!     * * * * * * * * * * * * * * * * * * * * * * * * *
!     CALCULATE MONTHLY FLOWS:
!     * * * * * * * * * * * * * * * * * * * * * * * * *
      do l=1,no
	  if(nqevent(l).gt.1)then
          qhydevent(l)=qhydevent(l)/nqevent(l)
          qsynevent(l)=qsynevent(l)/nqevent(l)
        else
          qhydevent(l)=-999.0
          qsynevent(l)=-999.0
        endif
      end do


d	if(iopt.eq.2)print*,' checkpoint 70 in lst'

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     THIS SECTION PRINTS THE ERROR FIELD ON GAUGED WATERSHEDS ONLY
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!     rev. 9.7.09  Oct.  11/10  - NK: update flowflag in lst.f for subsequent events
!     update the flow flag incase flows were found in later events

      do l=1,no
	  suberr(l)=statnerr(l)
	end do

c!     rev. 9.7.01  Jun.  09/10  - NK: fixed error.xyz & error.r2s

!     pause 103

!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
c      if(id.eq.ni.and.abs(dds_flag).ne.1)then  ! this part only when all done after last event

      if(writeflg(1))then
      if(id.eq.ni.and.dds_flag.ne.-1)then  ! this part only when all done after last event
c        if(dsflg.eq.'n')then
c          write(98,*)'++++++++++ WARNING++++++++++'
c          write(98,*)'Subarea assumed= 1.0 where -ve areas calculated'
c          write(98,*)'Probable cause: '
c          write(98,*)' two or more stream gaugesin one grid'
c          write(98,*)'                 or'
c          write(98,*)' stations not ordered in the downstream direction'
c          write(98,*)'>>>>>>Fatal mess for PAF file<<<<<<<'
c          write(98,*)
c          write(98,*)
c          write(*,*)'++++++++++ WARNING++++++++++'
c          write(*,*)'Subarea assumed= 1.0 where -ve areas calculated'
c          write(*,*)'Probable cause: '
c          write(*,*)' two or more stream gaugesin one grid'
c          write(*,*)'                 or'
c          write(*,*)' stations not ordered in the downstream direction'
c          write(*,*)'>>>>>>Fatal mess for PAF file<<<<<<<'
c          write(*,*)
c          write(*,*)'See spl.err file for more details'
c          write(*,*)
c        endif

!       WRITE THE RUNOFF ERRORS IN A GRIDDED FORMAT FOR THE 
!       srfERROR.TXT FILE 
!           unit 66 = simout\error.xyz
!           changed to Nash efficiency Dec. 12/15 NK  nash_eff.r2c
c        open(unit=66,file=filename(66),status='unknown',iostat=ios)
c        if(ios.ne.0)then    ! added Nov. 10/14  nk
c          print*
c          print*,'Unable to open file',filename(66)(1:40)
c          print*,'Possible causes:'
c          print*,'file in use by another application'
c          print*,'or target directory does not exist'
c          stop 'Program aborted in lst.f @ 744'
c        endif
!         write(66,6676)hrs,mins,secs,day,month,year
c        do i=1,ycount
c          do j=1,xcount
c            if(basinerr(i,j).lt.999.0)then
c              if(llflg.ne.'y')then
cc               ig,jg = utm coordinates
c                aig=float(iymin+(i-1)*istep)+float(istep)/2.0
c                ajg=float(jxmin+(j-1)*istep)+float(istep)/2.0
c!                write(66,6093)nbasin(i,j),i,j,ig,jg,basinerr(i,j),
c!     *                        precadj(i,j)
c              else
cc               ig,jg = lat-long in degrees.00
c                aig=(float(iymin)+float(i-1)*grdn+grdn/2.0)/60.0
c                ajg=(float(jxmin)+float(j-1)*grde+grde/2.0)/60.0
c!                write(66,6094)nbasin(i,j),i,j,aig,ajg,basinerr(i,j),
c!     *                        precadj(i,j)
c              endif
c              write(66,6095)ajg,aig,basinerr(i,j)
c            endif
c          end do
c        end do
c        close(unit=66,status='keep')

!       rev. 9.7.01  Jun.  09/10  - NK: fixed error.xyz & error.r2s
!       rev. 9.7.09  Oct.  11/10  - NK: update flowflag in lst.f for subsequent events
!       update the flow flag incase flows were found in later events
!        check to see if all stations have and error value
!       added a bunch of checks Dec. 15/10 nk
        do i=1,no
          if(inbsnflg(i).eq.1)then
            do l=1,no
              if(inbsnflg(l).eq.1.and.ds_sta(l).gt.0)then
                if(.not.flowflag(l).and.flowflag(ds_sta(l)))then
                  statnerr(l)=statnerr(ds_sta(l))
                  flowflag(l)=.true.
                endif  
	        endif
  	      end do
  	    endif
	  end do
         
!       NOT CORRECTED FOR U/S ERRORS < < < < < < < < <   
!       to get error for each sub basin use nudging 
        do i=1,ycount
          do j=1,xcount
            basinerr(i,j)=-999.000
            nasheff(i,j)=-999.000
          end do
        end do
        do n=1,naa
          i=yyy(n)
          j=xxx(n)
          if(nbasin(i,j).ne.0)then
            l=nbasin(i,j)
!           WE ARE IN A GAUGED SUB-AREA
            basinerr(i,j)=statnerr(l)
!     rev. 10.1.12 Dec.  12/15  - NK: Added Nash Efficiency nasheff.r2c file unit-66
            nasheff(i,j)=rsquare(l)    ! Nash efficiency actually
          endif
        end do

!     rev. 9.3.07  Dec.  29/06  - NK: added error field for whole domain
        author='WATFLOOD/CHARM                          '
        name='Volume errors in %                      '
        coordsys_temp=coordsys1
        zone_temp=zone1
	  datum_temp=datum1
	  xorigin_temp=xorigin
	  yorigin_temp=yorigin
	  xcount_temp=xcount
	  ycount_temp=ycount
	  xdelta_temp=xdelta
	  ydelta_temp=ydelta
	  startdate='unknown   '
	  starttime='unknown   '
        unit_conversion=1.0
          attribute_count=1
	    attribute_name='error                                    '
	    attribute_units='percent                                 ' 
        source_file_name='last spl run'     
        do j=1,xcount
	    do i=1,ycount
            outarray(i,j)=amin1(basinerr(i,j),999.000)
            outarray(i,j)=amax1(basinerr(i,j),-999.000)
	    end do
	  end do
        
!     rev. 9.8.45  Jan.  31/13  - NK: disabled some writes for iopt = 99
        if(dds_flag.ne.1.and.iopt.lt.99)then
!         write the header for error.r2c
!       written in sub for ddsflg=1
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(67,67,0,0,0,0,11)   
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         write the data
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(67,67,0,0,0,1,11)   
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  endif


!     rev. 10.1.12 Dec.  12/15  - NK: Added Nash Efficiency nasheff.r2c file unit-66
        author='WATFLOOD/CHARM                          '
        name='Nash_Efficiency                         '
        coordsys_temp=coordsys1
	    attribute_name='Nash Efficiency                          '
	    attribute_units='dimensionless                            ' 
        source_file_name='last CHARM run                          '
        do j=1,xcount
	    do i=1,ycount
            outarray(i,j)=nasheff(i,j)
	    end do
	  end do
        
!     rev. 9.8.45  Jan.  31/13  - NK: disabled some writes for iopt = 99
        if(dds_flag.ne.1.and.iopt.lt.99)then
!         write the header for error.r2c
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(66,66,0,0,0,0,11)   
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         write the data
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call write_r2c(66,66,0,0,0,1,11)   
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  endif

!     rev. 9.9.39  Nov.  14/14  - NK: Modifications for watroute
      if(modelflg.eq.'n')then     ! this section not needed for watroute
!     rev. 9.3.07  Dec.  29/06  - NK: added sum_precip for whole domain
        fln(99)='sum_precip.r2c'
        author='watflood                                '
        name='Model run summed precip                 '
        unit_conversion=1.0
	  attribute_name='sum precipitation                       '
	  attribute_units='mm                                      ' 
!        attribute_type='Runoff                                  '  
        source_file_name='last spl run'     
        do j=1,xcount
	    do i=1,ycount
c            outarray(i,j)=sum_precip(i,j)
	    end do
	  end do
        
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
!       write the header
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(99,99,1,0,0,0,10)   
        call write_r2c(99,99,0,1,0,0,10)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       write the datawrite_r2c
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(99,99,1,1,1,1,10)   
        call write_r2c(99,99,0,1,0,1,10)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        fln(99)='sum_runoff.r2c'
        author='watflood                                '
        name='Model run summed runoff                 '
        unit_conversion=1.0
	  attribute_name='sum runoff                              '
	  attribute_units='mm                                      ' 
!        attribute_type='Runoff                                  '  
        source_file_name='last spl run'     
        do j=1,xcount
	    do i=1,ycount
              outarray(i,j)=0.0
              n=s(i,j)
              if(n.gt.0)outarray(i,j)=QQsum(n)/da(n)
	    end do
	  end do
        
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
!       write the header
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(99,99,1,0,0,0,10)   
        call write_r2c(99,99,0,1,0,0,10)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       write the datawrite_r2c
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        call write_r2c(99,99,1,1,1,1,10)   
        call write_r2c(99,99,0,1,0,1,10)   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
d         if(iopt.eq.2)then
d           print*,' In lst before writing to newerror.txt'
d           print*
d         endif

d	  if(iopt.eq.2)print*,' checkpoint 80 in lst'

         if(dds.eq.0)then
           do i=ycount,1,-1
              do j=1,xcount
	          basinerr(i,j)=amin1(999.0,basinerr(i,j))
	          basinerr(i,j)=amax1(-999.0,basinerr(i,j))
c	          basinerr(i,j)=amax1(0.0,basinerr(i,j))
                if(abs(basinerr(i,j)).ge.999.0)basinerr(i,j)=0.0
              end do
           end do
           if(pafflg.eq.'y')then
!            WRITE THE GRIDDED FLOW ERRORS FOR EACH WATERSHED
             INQUIRE(FILE='paf.r2s',EXIST=exists)
             if(.not.exists)then
               INQUIRE(FILE='newerror.txt',EXIST=exists)
!              write a newerror.txt file only when an error file already exists.   
               open(unit=99,file='newerror.txt',status='unknown',
     *                                               iostat=ios)
               if(ios.ne.0)then    ! added Nov. 10/14  nk
                 print*
                 print*,'Unable to open file newerror.txt'
                 print*,'Possible causes:'
                 print*,'file in use by another application'
                 stop 'Program aborted in lst.f @ 908'
               endif
               if(exists)then
                 do while(.not.eof(99))
                   read(99,*,iostat=ios)line
                 end do
                 backspace 99
                 backspace 99
                 read(99,*)nblock
                 backspace 99
               else 
                 nblock=0
               endif
	      
               print*
	         print*,'A newerror.txt file with',nblock
	         print*,'blocks of data has been found '
	         print*,'Do you want to continue with this old file'
	         print*,'and add to it ...`y`'
	         print*,'or delete this file ... `n`?'
	         print*,'Answer:'
	         read*,answer
	         print*,answer
	         if(answer.eq.'n')then
                 close(unit=99,status='delete')
                 print*,'old mewerror.txt file deleted'                        
    	           nblock=0   ! start a new newerror.txt file
               else
                 print*,'To continue the sequence of refining the PAF'
                 print*,'copy newerror.txt to error.txt and rerun the'
                 print*,'the program (SPLXnn)'
                 print*,'This will create a newpaf.r2s file'
                 print*,'which can be copied to paf.r2s and used as a '
                 print*,'permanent PAF'
                 print*,'If a paf.r2s file is found, it will take'
                 print*,'precedence over the error.txt file'
               endif	      

               nblock=abs(nblock)+1
               write(99,6670)nblock,time(1:2),time(3:4),time(5:6),
     *          cday(1:4),cday(5:6),cday(7:8)
!     rev. 9.5.69  Oct.  10/09  - NK: added xcount & ycount to error & paf files
               write(99,*)':xcount',xcount
	         write(99,*)':ycount',ycount
               do i=ycount,1,-1
                 write(99,6090)(basinerr(i,j),j=1,xcount)
               end do
               write(99,*)-1*nblock    ! no more data
               close(unit=99,status='keep')
               print*,'newerror.txt file has been written'
               print*
             endif
      
           endif
         endif

         endif    ! modelflg.eq.'n'


d	   if(iopt.eq.2)print*,' checkpoint 90 in lst'

!        * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!        WRITE THE  FLOW ERRORS .xyz file FOR EACH WATERSHED:
!        * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         open(unit=66,file=filename(66),status='unknown',iostat=ios)
c         if(ios.ne.0)then    ! added Nov. 10/14  nk
c           print*
c           print*,'Unable to open file',filename(66)(1:40)
c           print*,'Possible causes:'
c           print*,'file in use by another application'
c           print*,'or target directory does not exist'
c           stop 'Program aborted in lst.f @ 974'
c         endif
c!        changed filename and unit number to fln list mar02/04  nk
c      
c!        WRITE ..\simout\error.xyz   
c         if(llflg.ne.'y')then
c           istep=int(astep)
c           do i=1,ycount
c             do j=1,xcount
c               if(nbasin(i,j).gt.0)then
c                 ig=iymin*1000+(i-1)*istep*1000
c                 jg=jxmin*1000+(j-1)*istep*1000
c                 write(66,207)float(jg),float(ig),basinerr(i,j)
c!                 write(*,207)float(jg),float(ig),basinerr(i,j)
c               endif
c             end do
c           end do
c         else
c           istep=int(astep)
c           do i=1,ycount
c             do j=1,xcount
c               if(nbasin(i,j).gt.0)then
c                 write(66,201)(float(jxmin)+float(j-1)*grde)/60.0,
c     *                 (float(iymin)+float(i-1)*grdn)/60.0,
c     *                 basinerr(i,j)
c               endif
c             end do
c           end do
c         endif
c         close(unit=66,status='keep')

        endif   !id=ni
        endif   ! if(writeflg(1))


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     NOW WRITE A COMMA SEPERATED FILE FOR A SPREADSHEET:  spl.csv
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

d      if(iopt.eq.2)then
d        print*,' In lst before reading to unit 60'
d        print*
d      endif

!     Append to existing files
!     FOR THE RESUME CASE, GO TO THE END OF THE output FILEs: 
!     rev. 9.8.27  Sep.  27/12  - NK: changed action on resumflg='s' - keep tbcflg='y'
      if(resumflg.ne.'n'.and.contflg.eq.'y'.and.id.eq.1)then
        ios=0
        do while (ios.eq.0)
          read(51,6040,iostat=ios)tempflg  ! simout\spl.txt
        end do
        ios=0
        do while(ios.eq.0)
          read(60,6040,iostat=ios)tempflg  ! simout\spl.csv
        end do
        ios=0
        do while(ios.eq.0)       
          read(68,6040,iostat=ios)tempflg  ! simout\wetland.csv
        end do
        ios=0
        do while(ios.eq.0)
          read(73,6040,iostat=ios)tempflg  ! simout\resin.csv
        end do
      endif

d      if(iopt.eq.2)then
d        print*,' In lst before writing to unit 60'
d        print*
d      endif

!      write(60,6400)
!     FOR SMALL RECORD LENGTH


!     for grid shifting, we need sequential k's for spaghetti plots
      if(igrdshft.ne.0.and.id.eq.1)then
         kfirst=0
      endif

!     in the future, this section can be used to print out flows only when computed.
!
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     PRINT spl.csv and EVENT MEANS
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

! martin 
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
!     rev. 9.9.66  Apr.  29/15  - NK: Deleted mid_file headers in with tbcflg=y
      if(writeflg(60))then
          if(id.eq.1.and.hdrflg0.and.contflg.eq.'n')then 
!        write the headers for spl.csv
              write(60,60001)start_date,
     *            (gage(l),'_obs',gage(l),'_SIM',l=1,no)
60001         format(a12,<no>(',',a12,a4,',',a12,a4))     
          endif
 
!      if(kt.ge.irdt)then
!     NOTE:
!     qsyn & qhyd were turned into the mean flow for the time step above
!     For WATROUTE  output is wrt.csv and NIT spl.csv

        if(kt.ge.1)then
          do k=kt,nl,kt
            kfirst=kfirst+kt
c                  write(800,*)nl,kt,k,kfirst
            write(60,7401)kfirst,(qhyd(l,k),qsyn(l,k),l=1,no)
          end do
          write(75,7401)year_now,(qhydevent(l),qsynevent(l),l=1,no)
        else
          do k=irdt,nl,irdt
            kfirst=kfirst+irdt*kt
c              write(800,*)nl,kt,k,kfirst
            write(60,7401)kfirst,(qhyd(l,k),qsyn(l,k),l=1,no)
          end do
          write(75,7501)(qhydevent(l),qsynevent(l),l=1,no)
        endif
      endif
      
!     rev. 10.2.15 Feb.  05/18  - NK: Added 'results\monthly_peaks'
      if(iopt99)then
        do j=1,mhtot/732+1
!     rev. 10.2.41 Dec.  10/18  - NK: Added winter monthly peaks
          if(j.le.4.or.j.eq.12)then
            write(21,31000)j,(qp_month_hyd(l,j),
     *                             qp_month_syn(l,j),l=1,no)
          endif
          write(22,31000)j,(qp_month_hyd(l,j),qp_month_syn(l,j),l=1,no)
31000     format(i10,<2*no>f12.3)
        end do
      endif
      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(xmlflg.eq.'y')then
        call write_xml
      endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      

!     Append this year's natural flow to the file:
      if(MRBflg)then
        open(unit=99,file='MRB_WILLISTON_INFLOWS.new',status='old')
!     First go to the end of the file
        n=-23    ! subtract # header lines
        do while(.not.eof(99))
          read(99,*)line
          n=n+1
        end do
!     Then add the nest year's worth of data
        do k=1,nl/24       ! assuming always 24 hour data
          write(99,65001)(MRBdata(i,k+n),i=1,2),qsyn(33,k*24)
65001     format(3f10.3)     
        end do      
        close(unit=99,status='keep')
      endif

!     rev. 9.2.19  Oct.  28/05  - NK: Compute daily & monthly flows
!     rev. 9.7.25  Apr.  28/11  - NK: Fixed daily flows

!     calculate and write the daily flows to ..\simout\spl_dly.csv
c      if(kt.lt.24)then    won't work for other dt's
      if(kt.eq.1.and.writeflg(71))then
        do l=1,no
          do j=1,nl/24*24,24
            qhyd_dly_sum(l)=0.0
            qsyn_dly_sum(l)=0.0
	      i=0
	      ii=0
	      tmpflg=.true.
            do k=j,j+23,kt
              if(qhyd(l,k).gt.0)then
                qhyd_dly_sum(l)=qhyd_dly_sum(l)+qhyd(l,k)
	          ii=ii+1
	        else
	          qhyd_dly_sum(l)=0.0
	          tmpflg=.false.
	        endif
              qsyn_dly_sum(l)=qsyn_dly_sum(l)+qsyn(l,k)
	        i=i+1
            end do
	      if(ii.gt.0)then
              qhyd_dly(l,(j+23)/24)=qhyd_dly_sum(l)/float(ii)
	      else
              qhyd_dly(l,(j+23)/24)=-1.0  ! set -1 if data = missing
	      endif
            qsyn_dly(l,(j+23)/24)=qsyn_dly_sum(l)/float(i)
          end do
        end do

        if(writeflg(71))then
          do j=1,nl/24,kt
            write(71,7404)j,(qhyd_dly(l,j),qsyn_dly(l,j),l=1,no)
          end do
        endif
	else
!       get rid of any hourly file
        inquire(FILE=filename(71),EXIST=exists)  
	  if(exists)close(unit=71,status='delete')      
      endif
      
      
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
      if(frcflg.eq.'y'.and.allocated(iso_18O))then
        if(.not.allocated(iso_sumO_18O))then    !changed by NK Jan 22/16
c        if(firstpass.eq.'y')then
          open(unit=82,file=filename(82),status='unknown')
          open(unit=83,file=filename(83),status='unknown')
          allocate(iso_sumO_18O(n18O),iso_sumO_2H(n18O),
     *             iso_sum_18O(n18O),iso_sum_2H(n18O),
     *             iso_msim2_18O(n18O),iso_msim2_2H(n18O),
     *             iso_obs2_18O(n18O),iso_obs2_2H(n18O),
     *             iso_mobs2_18O(n18O),iso_mobs2_2H(n18O),
     *             iso_obs_18O(n18O),iso_obs_2H(n18O),
     *             iso_csum_18O(n18O),iso_csum_2H(n18O),     
     *             iso_rms_18O(n18O),iso_rms_2H(n18O),
     *             iso_n_18O(n18O),iso_n_2H(n18O),
     *             iso_omobs2_18O(n18O),iso_omobs2_2H(n18O),
     *         stat=iAllocate)
          if (iAllocate.ne.0) STOP 
     *    'Warning: error with allocation of iso error arrays in lst'
          do j=1,n18O
            iso_sumO_18O(j)=0.0
            iso_sumO_2H(j)=0.0
            iso_sum_18O(j)=0.0
            iso_sum_2H(j)=0.0
            iso_msim2_18O(j)=0.0
            iso_msim2_2H(j)=0.0
            iso_obs2_18O(j)=0.0
            iso_obs2_2H(j)=0.0
            iso_mobs2_18O(j)=0.0
            iso_mobs2_2H(j)=0.0
            iso_obs_18O(j)=0.0
            iso_obs_2H(j)=0.0
            iso_csum_18O(j)=0.0
            iso_csum_2H(j)=0.0
            iso_n_18O(j)=0
            iso_n_2H(j)=0
            iso_omobs2_18O(j)=0.0
            iso_omobs2_2H(j)=0.0
          end do
        endif

!     rev. 10.1.88 May   23/17  - NK: Fixed Julian_day problems for iso R/W
c        if(id.eq.ni.or.jul_day_now.ge.365)then
        if(jul_day_now.eq.1.and.hour_now.eq.24)then
            
!         write the isotope data after the last event or on the last day of the year
!         as the recorded isotope file is always for a whole calendar year.        

          if(leapyear)then  
            jd=366
          else
            jd=365
          endif
              
          do i=1,jd
            write(82,82001)i,(iso_18O(i,j),dstr_18O(i,j),j=1,n18O)
82001       format(i10,<2*n18O>(',',f10.2))          
          end do  
          
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
          if(flg2H.eq.2)then
          do i=1,jd
            write(83,82001)i,(iso_2H(i,j),dstr_2h(i,j),j=1,n2H)
          end do 
          endif 
          
 !     rev. 10.1.05 Oct.  11/15  - NK: Iso RMS error 
          write(81,*)
          write(81,*)'Item      Jul_day     location        #',
     *               '   Observed       Computed'     
          do j=1,n18O
            do i=1,jd
 !             select only days with data for error calculation            
              if(iso_18O(i,j).gt.-90000.0)then
!               accumulate sums for average values                                    
                iso_sumO_18O(j)=iso_sumO_18O(j)+iso_18O(i,j)
                iso_sum_18O(j)=iso_sum_18O(j)+dstr_18O(i,j)  
                iso_obs2_18O(j)=iso_obs2_18O(j)+
     *                          (dstr_18O(i,j)-iso_18O(i,j))**2
                iso_obs_18O(j)=iso_obs_18O(j)+
     *                          (dstr_18O(i,j)-iso_18O(i,j))
                iso_n_18O(j)=iso_n_18O(j)+1
                write(81,*)'18O',i,j,
     *                   iso_n_18O(j),iso_18O(i,j),dstr_18O(i,j)
              endif
            end do
          end do
          write(81,*)
          if(flg2H.eq.2)then
          write(81,*)'Item     Jul_day     location        #',
     *               '   Observed       Computed'     
          do j=1,n2H
            do i=1,jd
!             select only days with data for error calculation            
              if(iso_2H(i,j).gt.-90000.0)then
!               accumulate sum of squares              
                iso_sumO_2H(j)=iso_sumO_2H(j)+iso_2H(i,j)
                iso_sum_2H(j)=iso_sum_2H(j)+dstr_2H(i,j)  
                iso_obs2_2H(j)=iso_obs2_2H(j)+
     *                          (dstr_2H(i,j)-iso_2H(i,j))**2
                iso_obs_2H(j)=iso_obs_2H(j)+
     *                          (dstr_2H(i,j)-iso_2H(i,j))  
                iso_n_2H(j)=iso_n_2H(j)+1
                write(81,*)'2H',i,j,
     *                  iso_n_2H(j),iso_2H(i,j),dstr_2H(i,j)
              endif
            end do
          end do
          endif  
!       calcualte the iso error up to this point        
          do j=1,n18O
            iso_rms_18O(j)=sqrt(iso_obs2_18O(j)/iso_n_18O(j))
            iso_rms_2H(j)=sqrt(iso_obs2_2H(j)/iso_n_2H(j))
          end do
          
          write(81,*)'Isotope RMS error 18O'
           write(81,*)'   Location   # readings   RMS error'
          do j=1,n18O
             write(81,*)j,iso_n_18O(j),iso_rms_18O(j)
          end do
          write(81,*)
           write(81,*)'Isotope RMS error 2H'
           write(81,*)'   Location   # readings   RMS error'
          do j=1,n2H
             write(81,*)j,iso_n_2H(j),iso_rms_2H(j)
          end do
           write(81,*),'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
           
c!     rev. 10.1.89 May   25/17  - NK: Added errflg = 11 for isotope DDS
c          if(id.eq.ni)then
c            write(*,*)
c            write(*,*)'Isotope RMS error 18O'
c            write(*,*)'   Location   # readings   RMS error'
c             do j=1,n18O
c             write(*,*)j,iso_n_18O(j),iso_rms_18O(j)
c            end do
c            write(*,*)
c            write(*,*)'Isotope RMS error 2H'
c            write(*,*)'   Location   # readings   RMS error'
c            do j=1,n2H
c              write(*,*)j,iso_n_2H(j),iso_rms_2H(j)
c            end do
c          endif
        endif
      endif   ! frcflg.eq.'y'
 !     END rev. 10.1.05 Oct.  11/15  - NK: Iso RMS error 

!     calculate and write the monthly flows to ..\simout\spl_mly.csv
!     We can do this only if the event is at least one month duration
!     Split the year in to 12 equal lenght events of 730 hour
      if(mhtot.gt.670.and.mhtot.lt.745)then
        nmonths=1
      else
        nmonths=min(mhtot/720,12)     ! approximately
!       so if there are fractions of months, the last month will be ignored        
      endif

!     write the monthly flows
      evt_hrs=mhtot
C      if(nmonths.eq.1)then    ! changed ge to eq   mar 15/06  nk
      if(nmonths.ge.1.and.iopt99)then    ! changed ge to eq   mar 15/06  nk
!       we have month-long events
        do l=1,no
          if(nmonths.eq.1)then
	      hour_no=0
            j=1
            qhyd_mly_sum(l)=0.0
            qsyn_mly_sum(l)=0.0
            do k=1,nl
	        hour_no=hour_no+1    ! count hours
              qhyd_mly_sum(l)=qhyd_mly_sum(l)+qhyd(l,k)
              qsyn_mly_sum(l)=qsyn_mly_sum(l)+qsyn(l,k)
            end do
            qhyd_mly(l,j)=amax1(0.0001,qhyd_mly_sum(l)/float(hour_no))
            qsyn_mly(l,j)=qsyn_mly_sum(l)/float(hour_no)
          else  
!          we have events > 1 month - yearly events assumed
!     rev. 9.5.56  Mar.  26/09  - NK: Fix bug with month in yearly events
	      hour_no=0
            do j=1,nmonths
c            do j=1,12
              qhyd_mly_sum(l)=0.0
              qsyn_mly_sum(l)=0.0
              hour_neg=0
!             run through the months
              do k=1,mohours(j)
	          hour_no=hour_no+1    ! count hours
	          if(qhyd(l,hour_no).lt.0.0)then    !<<<<<<<<<<<<<<<<<<<<<<<
                  hour_neg=hour_neg+1 ! count hours with no data
!                 dont plot data for month with missing data
                endif
                qhyd_mly_sum(l)=amax1(0.0,qhyd_mly_sum(l))+
     *                  			              qhyd(l,hour_no)
                qsyn_mly_sum(l)=qsyn_mly_sum(l)+qsyn(l,hour_no)
              end do
              qhyd_mly(l,j)=qhyd_mly_sum(l)/float(mohours(j))
              qsyn_mly(l,j)=qsyn_mly_sum(l)/float(mohours(j))
c	pause
            end do
          endif
        end do
        do j=1,nmonths
c        do j=1,12
          write(78,7401)j,(qhyd_mly(l,j),qsyn_mly(l,j),l=1,no)
        end do
        do j=1,12
          do l=1,no
            if(qhyd_mly(l,j).gt.0.0)then
               mnl_diff(l)=qsyn_mly(l,j)-qhyd_mly(l,j)
            else
               mnl_diff(l)=-1.0e10
            endif
          end do
            write(960,7401)j,(mnl_diff(l),l=1,no)

!     rev. 9.9.24  Aug.  20/14  - NK: Added monthly mean flow csv file spl_mly_nn.csv
c         write(960+j,7401)id,(qhyd_mly(l,j)*kt,qsyn_mly(l,j),l=1,no)
        end do
          write(961,7401)id,(qhyd_mly(l,1)*kt,qsyn_mly(l,1),l=1,no)
          write(962,7401)id,(qhyd_mly(l,2)*kt,qsyn_mly(l,2),l=1,no)
          write(963,7401)id,(qhyd_mly(l,3)*kt,qsyn_mly(l,3),l=1,no)
          write(964,7401)id,(qhyd_mly(l,4)*kt,qsyn_mly(l,4),l=1,no)
          write(965,7401)id,(qhyd_mly(l,5)*kt,qsyn_mly(l,5),l=1,no)
          write(966,7401)id,(qhyd_mly(l,6)*kt,qsyn_mly(l,6),l=1,no)
          write(967,7401)id,(qhyd_mly(l,7)*kt,qsyn_mly(l,7),l=1,no)
          write(968,7401)id,(qhyd_mly(l,8)*kt,qsyn_mly(l,8),l=1,no)
          write(969,7401)id,(qhyd_mly(l,9)*kt,qsyn_mly(l,9),l=1,no)
          write(970,7401)id,(qhyd_mly(l,10)*kt,qsyn_mly(l,10),l=1,no)
          write(971,7401)id,(qhyd_mly(l,11)*kt,qsyn_mly(l,11),l=1,no)
          write(972,7401)id,(qhyd_mly(l,12)*kt,qsyn_mly(l,12),l=1,no)
      else
        if(iopt99)write(78,*)'Events shorter than one month. No data'
      endif

d	if(iopt.eq.2)print*,' checkpoint 100 in lst'

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     NEW - REV. 7.41 - OUTPUT TO INPUT FILES . . . . .
!     Write the background streamflow file        strout.1
!     This file can be used to compare future runs to previous runs
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c********************************************************************
!     ver. 9.1 - sediment component
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
d      if(iopt.eq.2)then
d        print*,' In lst before writing to unit 69'
d        print*
d      endif

      if(sedflg.eq.'y')then
!       output the sediment output:

        if(id.eq.1)then
          col_name(1)=
     *'SC kg/m3,   SM ton,  NC kg/m3,    NM ton,  PC kg/m3,  PM ton'
          do l=1,no
            col_name(l)=col_name(1)
          end do
           WRITE(69,6678)(l,col_name(l),l=1,no)
6678       format('hour,',99(i3,a60,'  ,'))
       endif


        do 393 k = kt,nl,kt
           write(69,6404)K,(sedsyn(L,K),sedmss(L,K),
     *       nitsyn(L,K),nitmss(L,K),
     *       phssyn(L,K),phsmss(L,K),L=1,No)
  393   continue
      endif

d      if(iopt.eq.2)then
d        print*,' In lst before writing to unit 59'
d        print*
d      endif

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     CREATE THE STAGE HYDROGRAPHS FOR THE DAMAGE LOCATIONS:
!     THE QHYD VARIABLE IS NOW REUSED AS THE STAGE VARIABLE
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     output for stage plots      
      if(writeFlg(59))then
          open(unit=59,file=filename(59),status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
              print*
              print*,'Unable to open file',filename(59)(1:40)
              print*,'Possible causes:'
              print*,'file in use by another application'
              print*,'or target directory does not exist'
              stop 'Program aborted in lst.f @ 1285'
          endif
          write(59,1000)fln(10)
          write(59,1015)no,nl,mhtot,kt
          write(59,1004)ndam
          do l=1,ndam
               write(59,1100)iys(l),jxs(l),damage(l),(frcs(l,i),i=1,4),
     *                 datum(l)
          end do

!         WRITE THE NOTES:
          write(59,1098)(note(i),i=1,nnnote)
          write(59,6400)
          if(ndam.gt.0)then
              do k=kt,nl,kt
                  do l=1,ndam
                      i=iys(l)
                      j=jxs(l)
!       TS - SET QLOC TO MIN VALUE B/C OF ERRORS IN DEBUGGING (OCT.30/03)
	                if(qloc(l,k).lt.0.0)qloc(l,k)=0.0001
	                if(s(i,j).ne.0)then
c                         qhyd(l,k)=frcs(l,1)*qloc(l,k)**frcs(l,2)
                      else
	                    qloc(l,k)=0.0001
	                    qhyd(l,k)=0.0001
	                endif
                  end do
                  write(59,6403,iostat=ios)
     *                  k,(qloc(l,k),qhyd(l,k),l=1,ndam)
              end do
          endif
!      close(59,status='keep')
      endif   
      
	If(ios.ne.0)then
	  write(99,99001)
99001   format(' WARNING: problems writing to unit 59 stg.plt')
	endif

d      if(iopt.eq.2)then
d        print*,' In lst before writing to unit 70'
d        print*
d      endif


! fix fix fix fix fix fix fix fix fix fix fix fix fix fix fix fix fix fix fix fix fix 
! this section needs to be done by write_flow1d
! populate the array etc.

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     WRITE THE FILE FOR external DW routing in tb0 format
!     WRITE THE FILE FOR external DW routing in tb0 format
!     WRITE THE FILE FOR external DW routing in tb0 format
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


!     REV. 7.2 Sept. 19/94 - ADDED IREACH(N) FOR DWOPER or Flow1D INPUT 
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     REV. 10.1.15 Jan.  08/16  - NK: Custom coding for Mackenzie River Basin Hydraulic Model
!     REV. 10.1.28 Apr.  26/16  - NK: Fixed first day of output for master_inflows file 
!     maxr = no of reservoirs - i.e. don't write if there are no reservoirs/reaches
      if(maxr.gt.0.and.routeflg.eq.'q'.and.dds_flag.eq.0)then
!       averaging the flows for the day      
        do mh=24,mhtot,24
  	    do i=1,maxr
	      q24(i)=0.0
            do j=mh-23,mh
              q24(i)=q24(i)+qdwpr(i,j)

c      if(i.eq.1)print*,j,q24(i),qdwpr(i,j)

	      end do
            q24(i)=q24(i)/24.0

c      if(i.eq.1)print*,j,q24(i)

c      if(i.eq.1)pause 'lst.f'


	    end do
!         Note:the last 22 lakes are not to be written to this file but 
!         passed through rerout for natural lake routing
!         Note: write_flow1d_tb0 is not called here as we need to write
!         the file continuously for the whole run
!         tHE HEADER IS WRITTEN IN sub.F
!         this should be changed to call write_flow1d_tb0 if possible
          write(70,5001)(q24(i),i=1,126)  ! just the # of MRBHM nodes
        end do
	endif

!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     WRITE THE RESERVOIR INFLOW FILES 
!     DWOPER FLOWS IN COLUMN FORMAT FOR PLOTTING MAYBE 
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
d        if(iopt.eq.2)then
d          print*,' In lst before writing to unit 70'
d          print*
d        endif
      if(id.eq.1)then
        write(53,*)'In lst'
        write(53,*)'Reservoir inflows this event:'
	endif

      do mh=kt,mhtot,kt
 	  do i=1,maxr
	    q24(i)=0.0
          do j=mh-kt+1,mh
            q24(i)=q24(i)+qdwpr(i,j)
	    end do
          q24(i)=q24(i)/float(kt)
	  end do
        write(53,201)(q24(i),i=1,maxr)
      end do

!     rev. 9.6.03  Mar.  31/10  - NK: replaced leakage.dat by nbs.tb0
!     write the net basin supply for the great lakes
c      if(resname(1).eq.'Superior     ')then
c        do mh=kt,mhtot,kt
c          write(79,5001)qdwpr(1,mh),
c     *		qdwpr(2,mh)-qdwpr(1,mh),
c     *		qdwpr(3,mh)-qdwpr(2,mh),
c     *		qdwpr(4,mh)-qdwpr(3,mh),
c     *		qdwpr(5,mh)-qdwpr(4,mh)
c	  end do
c	endif


d	if(iopt.eq.2)print*,' checkpoint 110 in lst'
 
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     WRITE THE RESERVOIR FILE (RESIN.csv) FOR PLOTTING
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

!     MOVED FROM SUB  NOV.19/97
!     REV. 8.95 - Mar.  15/99 -  COMPUTED MEAN FLOWS FOR TIME INCREMENT
!     CHECKED BY SPREADSHEET APRIL 12/99 NK.

!     SET THE PEAKS EQUAL TO ZERO FOR EACH EVENT:
      do l=no+1,no+noresvi
c         dpeakh(l)=0.1e-10
c         dpeaks(l)=0.1e-10
         volhyd(l)=0.0
         volsyn(l)=0.0
      end do

!     rev. 9.1.79  Mar.  30/05  - NK: ktri to area2 for reservoir inflow dt
!       CALCULATE THE MEAN FLOW for the reservoir inflow timestep
	if(noresvi.gt.0)then    ! added Apr. 10/05 nk
!       only when there are reservoir inflows
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
        if(id.eq.1.and.hdrflg0)then
          write(73,60002)start_date,
     *    (',obs_',resnamei(i),',SIM_',resnamei(i),i=1,noresvi)
60002     format(a12,<noresvi>(a5,a12,a5,a12))     
        endif
        if(id.eq.1)then
          allocate(qdwpr_sum(noresvi),stat=iAllocate)
          if (iAllocate.ne.0) STOP 
     *    'Error: error w/allocation of lake elv min/max in lst @ 1392' 
        endif

!     rev. 9.9.57  Feb.  08/15  - NK: Fixed resv inflow output resin & lake_sd
        ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
c        do k=ktri,mhtot,ktri
        do k=1,mhtot,ktri
          do i=1,noresvi
!     rev. 9.5.81  Jan.  16/09  - NK: allow reservoirs outside watershed in resv file
            if(inbsnflg(no+i).eq.1)then
!     rev. 9.5.78  Nov.  04/09  - NK: matched resvin locations to reach numbers
	        l=resin_reach(i)   ! determined in read_resvin
!             calculate mean for time step
              qdwpr_sum(i)=0.0     ! qdwpr mean daoly
              do kk=k,k+ktri-1
                qdwpr_sum(i)=qdwpr_sum(i)+qdwpr(l,kk)
c      if(l.eq.1)write(835,*)kk,qdwpr(l,kk),qdwpr_sum(i)
              end do
              kk=kk-1
              
              qdwprmd(i)=qdwpr_sum(i)/float(ktri)
            endif
!           for pre-emption
            ashnum(no+i)=ashnum(no+i)+(qinfl(i,k)-qdwprmd(i))**2
          end do
!         write obs & comp inflow values to resin.csv
!     rev. 9.8.32  Oct.  19/12  - NK: Fixed format for resin.csv in lst.f
!     rev. 9.2.41  Jun.  15/06  - NK: changed the resin.txt file to resin.csv
c      write(835,*)kk,qdwprmd(1),qdwpr_sum(1)
c      write(835,*)
c          write(73,7405)kk,(qinfl(i,k+ktri-1),qdwprmd(i),i=1,noresvi)
          write(73,7405)int(totaltime)-mhtot+kk,
     *                   (qinfl(i,k+ktri-1),qdwprmd(i),i=1,noresvi)
        end do


!       add the diversion flow to the receiving grid for the continuity check
        do k=ktri,mhtot,ktri
          do i=1,noresvi
	      if(inbsnflg(no+i).eq.1)then
!     rev. 9.5.78  Nov.  04/09  - NK: matched resvin locations to reach numbers
	        l=resin_reach(i)   ! determined in read_resvin
              if(qinfl(l,k).gt.0.01)then
                volhyd(no+i)=volhyd(no+i)+qinfl(l,k)
                volsyn(no+i)=volsyn(no+i)+qdwpr(l,k)
              endif
	      endif
          end do
        end do
        do i=1,noresvi
          if(inbsnflg(no+i).eq.1)then
            volhyd(i)=volhyd(i)*float(ktri)*3.6e-03
            volsyn(i)=volsyn(i)*float(ktri)*3.6e-03
	    endif
	  end do
                                            
!     rev. 10.1.61 Jan.  03/17  - NK: Changed results\peaks.txt to write peak flows
c        write(76,7601)(dpeakh(l),dpeaks(l),l=1,no+noresvi)
        write(76,7601)id,(qpeakh(l),qpeaks(l),l=1,no+noresvi)
c        write(77,7602)(volhyd(l),volsyn(l),l=1,no+noresvi)
      endif

      if(id.eq.ni) write(51,6030)optim

d	if(iopt.eq.2)print*,' checkpoint 120 in lst'

!     rev. 9.1.67  Oct.  21/04  - NK; added unit 80 for lake_stor & lake_outflow

!     write the header for lake_sd.csv
      if(firstpass.eq.'y')then
          col_name(1)=
     *'  lake_elv,  lake_stor,     lake_inflow,        NBS    ,'
          col_name(2)=
     *'lake_outflow, net_lake_outflow, delta_stor,'
        do l=1,noresv*2,2
          col_name(l)=col_name(1)
	    col_name(l+1)=col_name(2)
        end do
c        if(.not.contflg.eq.'y'.and.iopt99)then
        if(.not.contflg.eq.'y'.and..not.netCDFflg.or.iopt99)then
          WRITE(80,6680)'     time,',((l+1)/2,
     *             col_name(l),col_name(l+1),l=1,noresv*2,2)
6680      format(a10,999(i5,a56,a44))
c        else
c            print*,'Set iopt>0 to get lake_sd.csv file'
        endif
!     rev. 9.9.12  Apr.  04/14  - NK: Added min & max lake_level output file
        allocate(lake_elv_min(noresv),lake_elv_max(noresv),
     *         lake_elv_min_av(noresv),lake_elv_max_av(noresv),   
     *         sum_min(noresv),sum_max(noresv),   
     *         lake_elv_min_id(noresv,ni),lake_elv_max_id(noresv,ni),
     *         stat=iAllocate)
        if (iAllocate.ne.0) STOP 
     *    'Error: error w/allocation of lake elv min/max in lst @ 1392' 
        do l=1,noresv
          lake_elv_min(l)=1.0e+31
          lake_elv_max(l)=-1.0E+31
          sum_min(l)=0.0
          sum_max(l)=0.0
          do i=1,ni   ! can not use id here as it would stop the run after 1 event
            lake_elv_min_id(l,i)=1.0e+31
            lake_elv_max_id(l,i)=-1.0E+31
          end do
        end do
!     rev. 9.9.58  Feb.  13/15  - NK: Added time column to levels.txt
!     rev. 10.1.53 Nov.  09/16  - NK: Changed levels.txt to levels.csv
       if(lvlflg.and.writeFlg(953))then
          write(953,95301)'total_time,  ',
     *               (gname_lvl(n),',',gname_lvl(n),',',n=1,nolvl)
95301     format(a13,<2*nolvl>(a13,a1)) 
        endif               
      endif

!     rev. 9.8.24  Aug.  07/12  - NK: Added reading yyyymmdd_lvl.tb0 for lake levels
!     rev. 9.8.41  Jan.  28/13  - NK: fixed bug in lst for level print statement
!     rev. 9.8.45  Jan.  31/13  - NK: disabled some writes for iopt = 99
c	if(dds_flag.ne.1.and.iopt.lt.99)then
c	if(iopt.lt.99)then  !Changed May 6/13 nk
	if(writeFlg(1))then  !Changed May 6/13 nk
        if(lvlflg)then
          do k=ktlvl,mhtot,ktlvl
            do l=1,nolvl
              i=lvl_reach(l)
              if(i.eq.0)then
!                give this value so file will be readable by grapher etc.              
c                lvl_calc(l,k)=lake_elv(i,k)
                lvl_calc(l,k)=-999.0
!               if the level station is not in a reach, the value of lvl_reach(l)
!               will be 0 and will not work as a subscript below. So we set a dummy 
!               reach number at noresv+1
                lvl_reach(l)=noresv+1
              endif
            end do
          end do

!     rev. 9.9.58  Feb.  13/15  - NK: Added time column to levels.txt
!     REV. 10.1.22 Jan.  25/16  - NK: Fixed flowinit for partial basins
!     rev. 10.1.53 Nov.  09/16  - NK: Changed levels.txt to levels.csv
!         don't write if any station is outside the sub-basin model  
          if(dds_flag.eq.0)then
            do k=ktlvl,mhtot,ktlvl
              write(953,95300)int(totaltime)-mhtot+k,',',
     *         (lvl_obs(l,k),',',lake_elv(lvl_reach(l),k),',',l=1,nolvl)
95300         format(i12,a1,<2*nolvl>(f13.3,a1))         
            end do   
          endif  
        endif
      endif   

!     rev. 9.9.12  Apr.  04/14  - NK: Added min & max lake_level output file
      if(ktlvl.eq.0)ktlvl=1
      if(ktlvl.gt.0)then
        do l=1,noresv
c          do k=ktlvl,mhtot,ktlvl
          do k=1,mhtot
!           find the overall min & max        
            lake_elv_min(l)=amin1(lake_elv(l,k),lake_elv_min(l))
            lake_elv_max(l)=amax1(lake_elv(l,k),lake_elv_max(l))
!           find the event min & max          
            lake_elv_min_id(l,id)=
     *              amin1(lake_elv(l,k),lake_elv_min_id(l,id))
            lake_elv_max_id(l,id)=
     *                amax1(lake_elv(l,k),lake_elv_max_id(l,id))
          end do
        end do
      endif

!     write the data for lake_sd.csv & NBS
      if(noresv.gt.0)then
        
c      if(dds_flag.ne.1.and.iopt99)then
      if(.not.netCDFflg.and.dds_flag.ne.1.or.iopt99)then           
c	    if(trcflg.eq.'y')then
	    if(trcflg.eq.' ')then
            do k=1,mhtot
              write(80,8001)(lake_elv(l,k),lake_stor(l,k),
     *          lake_inflow(l,k),net_lake_inflow(l,k),lake_outflow(l,k),
     *          net_lake_outflow(l,k),del_stor(l,k),l=1,noresv)
c     *         ,(isolakeGW(l,k),l=1,noresv)
8001          format(g14.6,9999(',',g14.6))
	      end do
          else   !   trcflg.eq.'n'

!     rev. 9.5.60  Sep.  01/09  - NK: added deltat_report for lake_sd.csv file
!     rev. 9.7.21. Mar.  07/11  - NK: Fixed delta_report for longer periods in lst
!     rev. 10.2.27 Jul.  08/18  - NK: replaced lake_inflow_sum with temp_flow_sum
            if(deltat_report.ne.1)then
!             calculate the average for the reporting period:            
              do l=1,noresv
                do k=1,mhtot,deltat_report
 	            temp_flow_sum(l)=0.0
                  i=0
                  if(k.le.mhtot)then
	              do j=k,k+deltat_report-1
	                if(j.le.mhtot)then
	                  i=i+1
                        temp_flow_sum(l)=
     *				       temp_flow_sum(l)+lake_inflow(l,j)
                      endif
                    end do
                    j=j-1     
                  endif
c                  lake_inflow(l,j)=temp_flow_sum(l)/deltat_report
		      end do
                lake_inflow(l,j)=temp_flow_sum(l)/float(i)
              end do

!             Get the average for the net-inflow:              
              do l=1,noresv
                do k=1,mhtot,deltat_report
 	            temp_flow_sum(l)=0.0
                  i=0
                  if(k.le.mhtot)then
	              do j=k,k+deltat_report-1
	                if(j.le.mhtot)then
	                  i=i+1
                        temp_flow_sum(l)=
     *				      temp_flow_sum(l)+net_lake_inflow(l,j)
c      if(l.eq.1)print*,k,j,net_lake_inflow(l,j),temp_flow_sum(l)                        
                      endif
                    end do     
                    j=j-1     
                    net_lake_inflow(l,j)=temp_flow_sum(l)/float(i)
c      if(l.eq.1)then
c          print*,k,j,net_lake_inflow(l,j),temp_flow_sum(l)/float(i) 
c          print*,float(i)
c      endif
                  endif
		      end do
              end do
              
!             calculate the outflow average for the reporting period:            
              do l=1,noresv
                do k=1,mhtot,deltat_report
 	            temp_flow_sum(l)=0.0
                  i=0
                  if(k.le.mhtot)then
	              do j=k,k+deltat_report-1
	                if(j.le.mhtot)then
	                  i=i+1
                        temp_flow_sum(l)=
     *				      temp_flow_sum(l)+lake_outflow(l,j)
c      if(l.eq.1)print*,k,j,i,lake_outflow(l,j),temp_flow_sum(l)                        
                      endif
                    end do 
                    j=j-1     
!     rev. 10.2.67 Nov.  03/19  - NK Fixed flow averaging in lst
                lake_outflow(l,j)=temp_flow_sum(l)/float(i)
c      if(l.eq.1)print*,k,j,i,lake_outflow(l,j)
c      if(l.eq.1)print*
                  endif    
 		      end do
              end do

!             calculate the net_outflow average for the reporting period:            
              do l=1,noresv
                do k=1,mhtot,deltat_report
  	            temp_flow_sum(l)=0.0
                  i=0
                  if(k.le.mhtot)then
	              do j=k,k+deltat_report-1
	                if(j.le.mhtot)then
	                  i=i+1
                        temp_flow_sum(l)=
     *				      temp_flow_sum(l)+net_lake_outflow(l,j)
                      endif
                    end do   
                    j=j-1     
!     rev. 10.2.67 Nov.  03/19  - NK Fixed flow averaging in lst
                net_lake_outflow(l,j)=temp_flow_sum(l)/float(i)
                  endif  
		      end do
              end do
	      endif     ! deltat_report.ne.1

c      print*,'deltat_report=',deltat_report
c      pause

!           write the data to lake_sd.csv
c            do k=1,mhtot,deltat_report
            do k=deltat_report,mhtot,deltat_report
	        time_sum=time_sum+deltat_report
!     rev. 10.2.03 Oct   28/17  - NK: Revert to old G format for lakeSD.csv
              write(80,8002)time_sum,(lake_elv(l,k),lake_stor(l,k),
     *          lake_inflow(l,k),net_lake_inflow(l,k),lake_outflow(l,k),
     *          net_lake_outflow(l,k),del_stor(l,k),l=1,noresv)
c8002          format(i10,<noresv>(',',g14.6))
8002          format(i10,999(',',g14.6))
c8002          format(i10,<noresv>(',',f14.3,',',f14.0,
c     *            ',',f14.3, ',',f14.3, ',',f14.3, ',',f14.3,',',f14.0))
c8002          format(i10,<noresv>(',',g14.3,',',g14.0,
c     *            ',',g14.3, ',',g14.3, ',',g14.3, ',',g14.3,',',g14.0))
            end do
            
!     REV. 10.1.36 Jul   12/16  - NK: Added results\LakeName.tb0
!     REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files

!           write the data to the lake_Resv_name.tb0 file
            if(noresv.gt.0.and.tb0flg.eq.'y')then
              do l=1,noresv            
                do k=deltat_report,mhtot,deltat_report
                  write(2000+l,3030)lake_elv(l,k),lake_inflow(l,k),
     *                              lake_outflow(l,k)
3030              format(14x,3f12.3)     
                end do
              end do
            endif
           
          endif   !trcflg

!     rev. 9.5.37  Oct.  14/08  - NK: added deltat_report to lake_sd.csv file write
!     rev. 9.5.37  UNDONE

        endif    !  dds_flag.ne.1
            
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
c        if(writeFlg(957))then
!       Write NBS data      
        if(nbsflg.eq.'y')then
          do k=deltat_report,mhtot,deltat_report
              write(957,95703)(net_lake_inflow(l,k),l=1,noresv)
95703         format(<noresv>f8.1)                  
          end do
          
!     rev. 10.2.65 Sep.  30/26  - NK Changed format of net_basin_supply output
!     This is intended to work nly for regl & glb 2 & 8 day forecasts  
!     When iether of the last 2 evets is larger than 100 days = 2400 hours,
!     we're probably not in a forcast event. So don't write NBS          
          if(id.eq.ni-1.and.mhtot.le.2400)then   ! regl forecast
              i=0
              do k=deltat_report,mhtot,deltat_report
                  i=i+1
                  do l=1,noresv
                      net_basin_supply(l,i)=net_lake_inflow(l,k)
                  end do
              end do
          elseif(id.eq.ni.and.mhtot.le.2400)then  ! glb forecast
              i=2
              do k=deltat_report,mhtot,deltat_report
                  i=i+1
                  do l=1,noresv
                      net_basin_supply(l,i)=net_lake_inflow(l,k)
                  end do
              end do
              do l=1,noresv
                  write(957,95704)l,resname(l),
     *                            (net_basin_supply(l,j),j=1,i)
c                  write(957,*)
95704             format(i5,5x,a12,999f8.1)                  
              end do
          endif
          if(mhtot.gt.2400)then
            print*
            print*,'Warning:'
            print*,'NBS file not written as last 1 or 2 events appear'
            print*,'larger than the regl or glb events.'
            print*
          endif
        endif
      endif    ! if(noresv.gt.0)then

!     write data to nbs.tb0
c     do k=24,mhtot,24
c        write(70,5001)(net_lake_inflow(l,k),l=1,noresv)
c      end do

      close(unit=27,status='keep')

      firstpass='n'

!     rev. 9.9.12  Apr.  04/14  - NK: Added min & max lake_level output file
c      if(id.eq.ni)then
      if(modelflg.eq.'n'.and.writeFlg(1))then
        open(unit=99,file='results\min_max_lake_elevations.txt',
     *          status='unknown',iostat=ios)
        if(ios.ne.0)then
          print*,'Unable to open the file'
          print*,'results\min_max_lake_elevations.txt'
          print*,'so written in the working directory instead'
          open(unit=99,file='min_max_lake_elevations.txt',
     *          status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file min_max_lake_elevations.txt'
            print*,'Possible causes:'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in lst.f @ 1693'
          endif
        endif
!       write the column headers        
        write(99,*)
     *   'Lake#    Min       MinAv       MaxAv         Max     ',
     *                       'RangeAv    RangeAll'
        do l=1,noresv
          if(id.eq.ni)then ! last event
            do i=1,ni
              sum_min(l)=sum_min(l)+lake_elv_min_id(l,i)
              sum_max(l)=sum_max(l)+lake_elv_max_id(l,i)
            end do
            lake_elv_min_av(l)=sum_min(l)/float(ni)
            lake_elv_max_av(l)=sum_max(l)/float(ni)
          else
            lake_elv_min_av(l)=-1.0
            lake_elv_max_av(l)=-2.0
          endif
          write(99,99002)l,lake_elv_min(l),lake_elv_min_av(l),
c          write(99,*)l,lake_elv_min(l),lake_elv_min_av(l),
     *          lake_elv_max_av(l),lake_elv_max(l),
     *          lake_elv_max_av(l)-lake_elv_min_av(l),
     *          lake_elv_max(l)-lake_elv_min(l)  
99002     format(i5,6g12.3)          
        end do
        close(unit=99,status='keep',iostat=ios)
      endif

! FORMATS

  201 format(999g12.3)
  205 format(3i5,f10.2)
  206 format('    l    i    j      error ycount,xcount:',2i5)   
  207 format(3(' ',f12.3))
  503 format(30f8.2)
 1000 format(' ',a20)
 1004 format(' ',i5,' damage sites:')
 1011 format(a12,20x,999f8.0)
 1012 format(' ',a12,f8.0,4f7.0,f7.0,f7.1,2g12.0,3f8.0,99i5)
 1013 format('  no  location      area  precip  o/ro <->c/ro c/ro(t)',
     *'   Dv%   Nash_E    qp/m     qp/c      WS_A    spl_A   %Diff',
     *'  Class fractions')
 1014 format(' location',26x,'precip   o        c       qpm     qpc')
 1015 format(5i5)
 1017 format(5x,'peak flow measured =',f8.1,' cms, computed =',f8.1)
 1018 format(' ',i5,1x,a12,f8.0,4f7.0,2f8.2,5f9.0,99f6.2)
 1019 format(' ',i5,1x,a12,f8.0,11f7.0,99f5.2)
 1098 format(a80)
 1100 format(' ',2i5,1x,a12,7x,4e10.3,f8.3)
 2001 format(i8,' ',3(f8.2,' '),i8,a12,20x)
 2002 format(8f10.3)
 5000 format(' reach # ',2i5)
 5001 format(12x,999f12.3)
 5002 format(/' reach no ',i5,' of ',i5)
 6000 format(5x,'grid length =',f6.0,' meters')
 6002 format('0',10x,'* * * time=',f5.2,' hours * * *',/)
 6011 format(5x,'summary',i2,'   scm(1-5)=',999f5.2/)

 6013 format(' ',2f12.3,i5,1x,a12,99f5.2)
 6014 format(' ',2i5,99f5.2)

 6015 format(' runtime  ',a2,':',a2,':',a2,2x,a4,'-',a2,'-',a2)
 6016 format(/' soil moistutes are ',5f5.2/)
 6020 format(' ','id/l/area(l)/',2i5,f10.3,
     *'    area(l) can not be less than 0.0',
     *'    please check stream gauge location for gauge l'/)
 6030 format(5x,'total average error this run is',e15.6,' cms')
 6040 format(8a)
 6070 format(5x,'the initial unit flow for storm no. =',i5,' is',f10
     *.4,' m**3/sec/km**2')
 6080 format(5x,'lst: the maximum calculated flows are:',/)
 6081 format(9x,'n',5x, 'yyy(n)',5x,'xxx(n)',5x,'da(n)',3x,
     *'qmax(n)',3x,'sump(n)',3x,'sumrff(n)',3x,'eloss(n)',/)
 6089 format(999g7.1)
 6090 format(999f7.0)
 6091 format(999f5.2)
 6092 format(999f5.0)
 6093 format(5(i5,','),5(f10.1,','))
 6094 format(3(i5,','),2(f12.3,','),5(f10.1,','))
 6095 format(2(f10.3,','),g10.2)
 6400 format('  hour    then pairs of measured and computed flows')
 6409 format('  hour    then measured flows')
 6410 format('  hour    then computed flows')
 6401 format(i8,18g12.2)
 6402 format(i8,18g12.1)
 6404 format(i8,<no*2>(',',e10.3))
 7401 format(i8,<no*2>(',',g12.4))
c 7402 format(i8,512(',',f12.4))
c 7403 format(i8,512(',',f12.4))
 7404 format(i8,<no*2>(',',g12.4))
 7405 format(i11,<noresvi*2>(',',g16.4))
 7501 format(f9.4,<no*2>(',',g12.4))
c 7502 format(f9.3,512(',',f9.3))
c 7503 format(f9.2,512(',',f9.2))
c 7504 format(f9.1,512(',',f9.1))
c7402 format(i4,',',60(f6.1,','))
 6403 format(i4,60(f8.2,f8.3))
 6667 format('  snow on each element in mm, scaled by ',f5.2)
 6668 format('  precip. on each element in mm, scaled by ',f5.2)
 6666 format(3i10,3f10.1) 
 6669 format('  final soil moisture for each element is:')
 6670 format(i5,' Errors in %.Runtime ',
     *    a2,':',a2,':',a2,2x,a4,'-',a2,'-',a2)
c6671 format('  l,next,subarea(l),suberr(l)/',2i5,4f9.1)
 6671 format(' table for surfer: ',2(i2,':'),i2,2x,2(i2,'/'),i4,/
     *'nbasin   i    j      N            E     error    precip',
     *'    runoff      evap  runoff coeff.')
 6672 format('  runoff from each grid in mm')
 6673 format('  losses from each grid in mm')
 6674 format('  runoff coefficient ')
 6675 format('nbasin,   i,    j,         N,            E,    precip,'
     *,'    runoff,      evap,  runoff coeff.')
 6676 format(' table for surfer: ',2(i2,':'),i2,2x,2(i2,'/'),i4,/
     *'nbasin    i     j    N    E      error')
 6677 format('  recharge in each grid in mm')
 6679 format(a60)
 7601 format(i8,80f8.1)
 7602 format(80E11.3)
! 7777 format(10f10.1)

      RETURN

      END SUBROUTINE lst
