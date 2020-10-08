      subroutine rd_wfo_spec()

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
     
!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      REAL(4)		::     wfo_spec_version_number
      INTEGER           :: nratio,nj,iallcnt1,m,i,j,ios,n,ii,
     *                     iallocate
      logical           :: exists

	DATA iallcnt1/0/

!     NOTE: FOR MONTHLY RUNS, DIMENSIONS CAN BE CHANGED FROM 
!           3,8784,500  TO  12,366,3000

!>>>>>>>>>>>>>  AB: STUFF FOR ENSIM
      INTEGER(4) :: wfo_yy,wfo_mm,wfo_dd,wfo_hh,wfo_mi,wfo_ss,
     *                   wfo_ms
      INTEGER(4) :: wfo_seq

      if(iopt.eq.2)print*,' In rdwfospec after definitions'

!     RESET SETS ALL INITIAL VARIABLES

!      wfo_open_flg='n'
                
      INQUIRE(FILE='wfo_spec.txt',EXIST=exists)
      IF(exists)THEN
        open(unit=99,file='wfo_spec.txt',status='old',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file  wfo_spec.txt'
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          stop 'Program aborted in rd_wfo_spec @ 40'
        endif
      else
!     rev. 10.2.06 Oct   28/17  - NK: wfo_spec.txt in working OR basin directory
        INQUIRE(FILE='basin\wfo_spec.txt',EXIST=exists)
        IF(exists)THEN
          open(unit=99,file='basin\wfo_spec.txt',
     *                       status='old',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file  basin\wfo_spec.txt'
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            stop 'Program aborted in rd_wfo_spec @ 40'
          endif
        else
          print*,' WFO_SPEC.TXT file not found in the working OR'
          print*,' basin directories'
          print*,' although ensim flag = y'
          print*,' Program continues with ensimflg set = n'
          print*
          print*,' BSN.EXE can create a wfo_spec.new file'
          print*,' that can be copied to the working directory'
          print*,' as wfo_spec.txt    '
          print*
          print*,' Hit Ctrl C to abort or'
          pause ' Hit enter to continue       (in rd_wfo @ 49)'
        endif  
      endif

!     rev. 9.1.50  Jan.  14/04  - NK: version number added to the wfo_spec.txt file
        read(99,99000,iostat=ios)wfo_spec_version_number
        if(ios.ne.0.or.wfo_spec_version_number.ne.6.0)then    ! see below <<<<<<!!!!
          print*,' Problem reading wfo_spec.txt'
          print*,' You may need to update the wfo_spec.txt file'
          print*,' File version number =',wfo_spec_version_number
          print*,' Required version number = 6.0'
	    print*,' As of Feb. 15/18'
          print*,' Run the current bsn.exe to create a new'
          print*,' basin\wfo_spec.new file and copy it to the'
          print*,' working directory as wfo_spec.txt'
          print*
          stop ' Program aborted in spl @ 1195'
        endif
        read(99,99004,iostat=ios)nj
!       This number comes from the bsn.for program
!       See "no of attributes"

        allocate(attname(nj),attunits(nj),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *   'Error with allocation of arrays in write_both_headers'

! Note:  number below needs to be updated if items are added <<<<<<<<<<<<<<<<<!!!!!
        
        if(nj.ne.15+15*(classcount))then
          print*,' The wfo_spec.txt file needs to be updated'
          print*,' additional attributes have been added'
          print*,' Re-run the updated bsn.exe on the map file ' 
          print*,'     and create a wfo_spec.new file'
          print*,' classcount=',classcount
          print*,' No entrees found =',nj
          print*,' Looking for', 15+15*(classcount)
          print*,' Copy wfo_spec.new to ..\wfo_spec.txt' 
          print*,' Revision 9.1.29'
          print*
          stop ' Program aborted in rd_wfo_spec @ 406'
        endif

!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
        if(.not.allocated(wfo_pick))then
          allocate(wfo_pick(nj),wfo_attributes(nj),
     *          wfo_sum_p(na),wfo_cum_p(na),
     *          wfo_sum_qsyn(na),wfo_sum_qhyd(na),stat=iAllocate)
          if(iAllocate.ne.0)then
	      print*,'nj=',nj,'  na=',na
            STOP
     *     'Error with allocation of area2 arrays in rd_wfo_spec @ 74'
          endif
	  endif

!     Cumulative precip is set to 0.0 only the first time that the wfo stuff
!     is initialized. So if you want the total precip over a period, you have
!     to have the wfo flag on for just that period.

c       read(99,99004,iostat=ios)ireport
        read(99,*,iostat=ios)ireport
        if(ios.ne.0)then
          print*,' Problem reading line two of the'
          print*,' wfo_spec.txt file'
          print*
          stop ' Program aborted in spl @ 1195'
        endif
        if(ireport.eq.0)ireport=1

!       ireport_start and ireport_end mark the beginning and ens of the 
!       watflood.wfo file It is based on the total time from the beginning 
!       of a simulation - not the event time.

c       read(99,99004,iostat=ios)ireport_start
        read(99,*,iostat=ios)ireport_start
        if(ios.ne.0)then
          print*,' Problem reading  line two of the'
          print*,' wfo_spec.txt file'
          print*
          stop ' Program aborted in spl @ 1195'
        endif

c       read(99,99004,iostat=ios)ireport_end
        read(99,*,iostat=ios)ireport_end
        if(ios.ne.0)then
          print*,' Problem reading  line three of the'
          print*,' wfo_spec.txt file'
          print*
          stop ' Program aborted in spl @ 1195'
        endif
        if(ireport_end.eq.0)ireport_end=8760000         ! 1000 years

        do j=1,nj
          read(99,99003,iostat=ios)wfo_pick(j),wfo_attributes(j)
          if(ios.ne.0)then
            print*,'iostat code =',ios
            print*,' Read to line ',j,' in wfo_pick.txt'
            print*,' then a problem was found'
            STOP 'program aborted in sub.for @ 1190'
          endif
        end do
        close(unit=99,status='keep')

!     this check added 24/10/05  nk
      if(trcflg.ne.'y'.and.wfo_pick(5).eq.1)then
	  wfo_pick(5)=0
	  print*
	  print*,'Warning:'
	  print*,'You can not use the groundwater output function'
	  print*,'for ENSIM if trcflg = `n` '
	  print*,'wfo_pick(5) is set to 0'
	  print*,'I.e. the ground water contribution for ENSIM is'
	  print*,'turned off.'
	  print*
	  print*,'If you want the GW contribution for ENSIM, please'
	  print*,'set the trcflg = `y` in the event file'
	  print*
	  pause 'Hit enter to continue with GW turned off'
	endif

      print*,'Finished reading the wfo_spec.txt file'
      print*

99000 format(f5.1)
99001 format(f25.0)
99003 format(i1,5x,a50)
99004 format(i5)

      return
      end subroutine rd_wfo_spec