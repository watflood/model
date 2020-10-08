      SUBROUTINE options(ix,e1,conv,scale,smc5,nhr,nhf)

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
     
!**********************************************************************
!   REV. 9.00   March 2000  - TS: CONVERSION TO FORTRAN 90
!   REV. 9.03   Nov   2000  - TS: ADDED WATFLOOD SWAMP ROUTING 
!   rev. 9.2.08  Jul.  29/05  - NK: opt work-around in options 
!	
!   REV		Nov/2005 - BT (Bryan Tolson).  Added DDS optimization option                      
!   rev. 9.2.24  Dec.  07/05  - BT: DDS optimization 
!     rev. 9.5.18  Mar.  03/08  - NK: added conv to options & sub argument list
!     rev. 9.5.44  Oct.  27/08  - NK: optlow= moved from options.f to compute_error.f
!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!**********************************************************************

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(128):: qstr
      CHARACTER(72) :: junk
!      CHARACTER(14) :: date
      CHARACTER(1)  :: smok 
      CHARACTER*14    :: date
	INTEGER    :: iallcnt,icnt(5),ndir(5),nchr,ix,ios,icase,
     *              iallocate,igrdshft,iyshiftmin,iyshiftmax,
     *              jxshiftmin,jxshiftmax,ishift,jshift,inum,jnum,
     *              l,iw1,iw2,iv,iflg,i,n,ii,j,nhr,nhf,jan
      REAL(4)    ::   smc5(16),errold(5),err(5),chng(5),best(5)
      REAL(4)    :: optlow,e1,scale,ddtemp,cc1,cc2,crit,conv,best1
	real(4)    :: optlast
      integer*2  :: result1,ntest

!      DATA ntest/-1358/qstr/'options'/nchr/7/
      DATA ntest/43814/qstr/'optionsss'/nchr/9/
	DATA iallcnt/0/


!     *** MHTOT MUST BE LESS THAN 360 ***

! OPENING UP THE EVAPORATION FILES MOVED HERE SO THAT VER AND FLGEVP2
! HAVE VALUES - FRANK S: OCT/99
! VER 7.9 IF VER>7.9 OPEN EVAPORATION FILES
      if(ver.ge.7.9.and.flgevp2.eq.3.0)then
        open(unit=49,file=fln(19),status='unknown',iostat=ios)
        if(ios.ne.0)then
          write(*,99132)fln(19)
          write(98,99132)fln(19)
99132     format(' Warning: Error opening or reading fln:',a30/
     *'   Probable cause: missing yymmdd.flx - radiation input file'/
     *'   This file is optional - used to compute evaporation using '/
     *'   the priestly taylor routine'/
     *'   OR: in config.sys have you set files=100 & buffers=50?'/)
!          call errormsg(ios)
          print*,'iostat code =',ios
          STOP 'program aborted in spl.for @ 380'
        endif

        open(unit=50,file=fln(20),status='unknown',iostat=ios)
        if(ios.ne.0)then
          write(*,99142)fln(20)
          write(98,99142)fln(20)
99142     format(' Warning: Error opening or reading fln:',a30/
     *'   Probable cause: missing yymmdd.flx - point radiation file'/
     *'   This file is optional - used to compute evaporation using '/
     *'   the priestly taylor routine'/
     *'   needed if resinflg=y in yymmdd.evt'/
     *'   OR: in config.sys have you set files=100 & buffers=50?'/)
!          call errrormsg(ios)
          print*,'iostat code =',ios
          STOP 'program aborted in spl.for @ 381'
        endif
      endif


d      if(iopt.eq.2)print*, ' In options: 3 - before call rdsdc'

!      if(snwflg.eq.'y') call rdsdc()

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     rev. 9.5.44  Oct.  27/08  - NK: removed code & obj modules for hasp & rainbow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!     remove for unix
c      call keychk(qstr,nchr,result1)
c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c      call userchk(ntest,result1)
c!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      totaltime=0.0
      mo=mo1
      icase=numa
      nnn=0


!! ALLOCATIONS MOVED HERE FROM BELOW SO AS TO MISS 7777 CONTINUE BYPASS
!! SUPPOSED TO BE ALLOCATED WITH INUM/JNUM, BUT NOT ALWAYS DEF'ND -> SO 500    
      if(iallcnt.eq.0)then

        if(numa.gt.50)then     ! added Apr. 5/06  nk.
          print*,'No of parameters chosen for optimization > 50)'
	    print*,'This too many and is not allowed'
	    print*
	    stop 'Program aborted in options @ 111'
	  endif

        if(numa.gt.0)then
!         TS - ALLOCATION OF AREA8A ARRAYS
          allocate(a(50),b(50),ddelta(50),checkl(50),checkh(50),
     *       ssave(50),les(50),ba(50),iclosl(50),iclosh(50),
     *       nsign(50),odelta(50),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *      'Error with allocation of area8a arrays in optionsa'
        else
          allocate(a(1),b(1),ddelta(1),checkl(1),checkh(1),
     *         ssave(1),les(1),
     *    ba(1),iclosl(1),iclosh(1),nsign(1),odelta(1),stat=iAllocate)
          if(iAllocate.ne.0) STOP
     *      'Error with allocation of area8a(a) arrays in optionsa'
        endif
	  iallcnt=1
	endif

!     icase = numa  -1 for smc opt; 0 for one run; >1 for par opt

!     rev. 9.2.24  Dec.  07/05  - BT: DDS optimization 
c	if(dds_flag.eq.1)go to 8888

      if(icase)9999,7777,8888
!     8888 -> pattern search
!     7777 -> regular run
!     9999 -> operation: precip adjustment factor      


! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     SUBROUTINE SUBA DOES ALL RUNOFF AND FLOW CALCULATIONS  

 7777 CONTINUE

!     REV. 8.98   July 15/99 - MET GRID SHIFTING FOR WEATHER MODELS 
!     READ THE GRID SHIFTING PARAMETERS:
      open(unit=99,file='grdshift.txt',status='old',iostat=ios)
      if(ios.eq.0)then
        read(99,9901,iostat=ios)igrdshft
        if(ios.ne.0)then
          print*,' problems reading the gridshift.txt file'
          print*,' ios = ',ios
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          STOP ' program aborted in options @ 148'
        endif
        write(51,9902)igrdshft
        read(99,9903)junk
        write(98,9903)junk
9903    format(a72)
        read(99,9901)
     *     iyshiftmin,iyshiftmax,jxshiftmin,jxshiftmax,ishift,jshift
        write(51,9901)
     *     iyshiftmin,iyshiftmax,jxshiftmin,jxshiftmax,ishift,jshift
        close(unit=99)
        optlow=0.0
      else
!        these added Nov. 25/04  nk
         iyshiftmax=0
         jxshiftmax=0
         iyshiftmin=0
         jxshiftmin=0
         ishift=0
         jshift=0
         igrdshft=0
      endif

      inum=ycount+abs(iyshiftmin)+abs(iyshiftmax)
      jnum=xcount+abs(jxshiftmin)+abs(jxshiftmax)

d	if(iopt.eq.2)print*,'in options, inum.jnum/',inum,jnum

!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
	sensitivityflg=.false.
      if(dds_flag.eq.-1)then
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  call sensitivity(jan,smc5,conv,scale,icase,smok,optlow,igrdshft)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return
      endif


!     WARNING: DON'T COMBINE NEXT SEGMENT WITH PREVIOUS SEGMENT 
!     BECAUSE IGRDSHFT CAN BE SET = 0 IN THE FILE (i.e. NO IOS ERROR)
      if(igrdshft.le.0)then
        maxn=1
!        nnn=nnn+1   changed this June 23/02 nk killed allocation in flowinit
d        if(iopt.eq.2)print*, ' In options: 4 - before call sub'

        jan=1
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        nnn=nnn+1
      else
        maxn=((iyshiftmax-iyshiftmin)/ishift+1)*
     *     ((jxshiftmax-jxshiftmin)/jshift+1)
        write(52,5201)
5201    format(' dy_shift, dx_shift,peak_flows')
        do iyshift=iyshiftmin,iyshiftmax,ishift
          do jxshift=jxshiftmin,jxshiftmax,jshift
d            if(iopt.eq.2)print*, ' In options: 5 - before call sub'
            jan=1
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d            if(iopt.eq.2)print*, ' In options: 6 - before call lst'
c        if(.not.netCDFflg.or.iopt99)then
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call lst(scale,igrdshft)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c        endif
            write(52,5202)iyshift,jxshift,(qpeaks(l),l=1,no)
5202        format(2i10,40f10.3)
            nnn=nnn+1
          end do
        end do
      endif

      RETURN

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

!              * * * OPTIMIZATION * * *

!8888  if(iw.eq.0)then
!        iw1=1
!        iw2=5
!      else
!        iw1=iw
!        iw2=iw
!      endif

8888  if(icase.gt.0.and.resumflg.ne.'y')then
        print*,' WARNING: You should use a resume.txt file'
        print*,'          to ensure that each opt run starts'
        print*,'          with the same initial conditions'
        print*,' Of course a spinup period is OK too'
        print*,' '
        write(*,'(A)',advance='no')
     *            '  Be warned! Be warned! Be warned!'
        print*,' '
!        read(*,*)
        write(98,10001)
10001   format(' WARNING: You should use a resume.txt file',
     *         '          to ensure that each opt run starts',
     *         '          with the same initial conditions')
      endif

!     ASSIGN INITIAL VALUES TO PARAMETERS TO BE OPTIMIZED:

      
      if(par_file_version.ge.par_file_version_latest)then
!       par file will be updated from rdpar
!       before user can proceed
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
        call par_init()
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
      endif

c        call write_par(99,1)
        call write_par_10(99,25)


        nnn=0
        nstart=0 
        optim=0.0
        optlow=0.1e+30
        optlast=-1.0

8000    CONTINUE
!       >>>>>>>>>> START OF OPTIMIZATION LOOP <<<<<<<<<

!       ASSIGN INITIAL VALUES TO PARAMETERS TO BE OPTIMIZED:
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
        call par_assign()
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

!       SUB ORGANIZES ALL HYDROLOGIC COMPUTATIONS

        optlast=optim
        optim=0.0
	  totaltime=0.0
        mo=mo1
	  rewind 58  ! spl.plt
	  rewind 60  ! spl.scv
	  rewind 80  ! lake_sd

d       if(iopt.eq.2)print*, ' In options: 7 - before call sub'
        jan=1
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!       A NEW PARAMETER FILE CALLED 'NEW.PAR' WILL BE WRITTEN
!       FOR THE LOWEST ERROR YIELDING PARAMETERS:
        if(optim.lt.optlow)then
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c          call write_par(99,1)
          call write_par_10(99,25)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          optlow=optim
        endif

c          print*,optim,optlow
   
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call opt(*8999,*8000)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

8999    RETURN
      
c	ENDIF ! DDS or PAttern option
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

 9999 CONTINUE

!     REAL-TIME OPTIMIZATION of initial soil moisture:

!     	       * * * * * * * * * * * *
      if(icase.eq.-1)then
!       OPTIMIZATION FOR SOIL MOISTURE  
        if(ni.gt.1)then
          write(6,1001)
          write(6,1002)ni
          RETURN
        endif
!       CALCULATE THE ERROR FOR THE INITIAL PARAMETERS:
        nnn=0
        do ii=1,classcount-2
           best(ii)=smc5(ii)
           err(ii)=0.0
           errold(ii)=0.1e+32
           chng(ii)=0.05
           if(smc5(ii).lt.chng(ii)) smc5(ii)=chng(ii)+0.01
           if(smc5(ii).gt.por-chng(ii)) smc5(ii)=por-chng(ii)-0.01
           icnt(ii)=0
           ndir(ii)=0	
        end do
        write(52,8003)smc5
!       SET DIRECTION FOR SEARCH AND SET NEW SMC:
!       * * * * * * * * * * * * * * * * * * * * *
        crit=0.0
d        if(iopt.eq.2)print*, ' In options: 10 - before call sub'
        jan=1
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 9015   call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       FIND THE ERROR ASSOCIATED WITH EACH PERMEABILITY CLASS
!       FOR EACH TRIAL SMC:
        nnn=nnn+1
!       NNN IS THE NNN'TH ITERATION
!        optim=0.0
        do l=1,no
           i=iy(l)	
           j=jx(l)
           n=s(i,j)
           ii=ibn(n)
c           optim=optim+delta(1,l)
           err(ii)=err(ii)+delta(1,l)
        end do
        write(52,8001)ii,(err(ii),ii=1,5)
!       a problem here with smc5
!       it needs to be somehow related to river classes of basin areas.
        do ii=1,classcount-2
!          IF ERROR IS REDUCED, KEEP GOING IN THE SAME DIRECTION
!          WE ARE ASSUMING A CONCAVE ERROR FUNCTION W/NO BUMPS
           if(err(ii).lt.errold(ii))then
              best(ii)=smc5(ii)
              smc5(ii)=smc5(ii)+chng(ii)
           else
              ndir(ii)=ndir(ii)+1
              if(ndir(ii).le.1)then
!                IF ERROR IS BIGGER, REVERSE DIRECTION
                 chng(ii)=-1.0*chng(ii)
              else 
!                AND CUT THE INTERVAL ON SECOND AND SUBSEQUENT 
!                DIRECTION CHANGES
                 chng(ii)=-0.4*chng(ii)
              endif
              smc5(ii)=smc5(ii)+chng(ii)
           endif
!          CHECK THE LIMITS AND STAY WITHIN
           if(smc5(ii).le.0.0)then
              chng(ii)=0.4* abs(chng(ii))
              smc5(ii)=0.0
           elseif(smc5(ii).ge.por)then
              chng(ii)=0.4*abs(chng(ii))
              smc5(ii)=por
           endif
           errold(ii)=err(ii)
           err(ii)=0.0
        end do
        write(52,8003)smc5
        write(52,8004)chng
        do ii=1,classcount-2
           if(smc5(ii).le.0.9*por)then
              crit=0.03*por
           else
              crit=0.01*por
           endif
           if(abs(chng(ii)).gt.crit)go to 9015
        end do
!         * * * * * * * * * * * * * * * * * * * * *
        do ii=1,classcount-2
           smc5(ii)=best(ii)
        end do
        write(52,6000)smc5
!       ICASEA IS SET AS A FLAG TO MAKE THE FINAL RUN FOR THE FULL 
!       FORECAST PERIOD
        icase=-11
!       UPDATE THE EVENT FILE:
        mo1=mo
        fln(20)='event/event.evt'
d        if(iopt.eq.2)print*, ' In options: 11 - before call outevt'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call outevt(conv,scale,smc5,nhr,nhf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      elseif(icase.eq.-2)then
      
!       * * * * RADAR SCALING (OPTIMIZATION) * * * *
!       ASSIGN INITIAL VALUES TO PARAMETERS TO BE OPTIMIZED:
        a(1)=1.0
        errflg=2
        numa=1
        ddelta(1)=0.25
        checkl(1)=0.1
        checkh(1)=3.
        nstart=0
        nnn=0
        optlow=0.1e+32
        optim=0.0
!       START OF OPTIMIZATION LOOP
 8080   CONTINUE
d        if(iopt.eq.2)print*, ' In options: 12 - before call sub'

        if(nnn.gt.0)then
          fln(99)='event/event.evt'
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call read_evt(date,conv,scale,smc5,nhr,nhf)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ensimflg='n'
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    call read_shed_ef(31,1)	
          call read_par_parser(32,2)    
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  endif
        optim=0.0
        scale=a(1)        
        errflg=2
        numa=1
        jan=1
        dds_flag=-1
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(optim.lt.optlow)then
          best1=scale
          optlow=optim
        endif

!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call opt(*8081,*8080)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

8081	continue

!       Done optimizing - do a final run
        fln(99)='event/event.evt'
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call read_evt(date,conv,scale,smc5,nhr,nhf)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ensimflg='n'
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        scale=best1
        numa=0
        jan=1
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        if(ni.le.1)call lst(scale,igrdshft)

      RETURN



!       NUMA IS SET AS A FLAG TO MAKE THE FINAL RUN FOR THE FULL  
!       FORECAST PERIOD
 8089   icase=-12
        scale=best1
!       UPDATE THE EVENT FILE:
        mo1=mo
c        fln(20)='event/event.evt'
c        call outevt(conv,scale,smc5,nhr,nhf)
      endif
!     	       * * * * * * * * * * * *
!     TO MAKE A FORECAST RUN OF NL HRS SET:
!     ENTRY POINT IF ICASE = -10 IN XXXX.PAR FILE
d      if(iopt.eq.2)print*, ' In options: 15 - before call sub'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call read_evt(date,conv,scale,smc5,nhr,nhf)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ensimflg='n'
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	call read_shed_ef(31,1)	
      call read_par_parser(32,2)    
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      scale=best1
      numa=0
      jan=1
      call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
      if(ni.le.1)call lst(scale,igrdshft)

      RETURN

! FORMATS
 1001 format
     *(' you cannot run smc optimization for more than 1 storm')
 1002 format('        you have',i3,' events')
 6000 format(' optimized smcs:'/5f10.3/)
 6500 format(6e10.3)
c 6501 format(' a(',i2,')too close to lower constraint')
 6502 format(' a(',i2,')too close to upper constraint')
c 6503 format(' ddtemp    cc2       checkl    checkh    a(i)')
 6510 format(' the constraints on a(',i3,') are too close together')
 8001 format(' err=',i5,5e10.3/)
 8003 format('    smc=',5f10.5)
 8004 format('   chng=',5f10.5)
 8123 format(' ','a(1),scale,best1,optim,optlow/',5f10.5)
 9345 format(' classcount,nbsn,numa/',3i5/)
 9901 format(10i5)
 9902 format(/' grid shift index =',i2,'  0=no shift  1=shifting')

      END SUBROUTINE options
