      subroutine rdresume()

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
     
!     rev. 9.2.31  Feb.  09/06  - NK: Added area check to rdresume
!     rev. 9.2.40  Jun.  09/06  - NK: added tto(),ttomin(),ttomax() to resume
!     rev. 9.8.87  Oct.  25/13  - NK: Added error message for mismatched resume file

      use area_watflood
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      real*4, dimension(:),   allocatable :: area_temp

      DIMENSION     :: smc5(16)
      CHARACTER(14) :: date
      CHARACTER(3)  :: eofmark
      CHARACTER(1)  :: lineflg,smok,change_flag
      character(20) :: junk
      character(256) :: line
      REAL(4)		:: optlow,time,tot1,qwert,conv,scale,&
                  smc5,tj1,clock,t,thr,dtmin,dtmax,div,aintvl,sintvl,&
                  tot2,tdum,qtemp,diff2,sdlz,&
                  wfo_spec_version_number
      INTEGER       :: rbin,inta,block,no1,classcount1,na1,&
                      ycount1,xcount1,&
                      iallcnt1,n1,ii1,jan,m,ios,n,iallocate,&
                      l,ii,juold,jj,lun,nhr,nhf,&
                      i,j,icase,iz,jz,nchr,mz,ju,mon,&
                      iwest,ieast,isouth,inorth,nocolumns,nno,k,&
                      nu,igrdshft,minsnwflg,oldjan,flgevp22,&
                      noresv1,nj,npick,n_version
      CHARACTER(128):: qstr
      character*15  :: junk1
      logical       :: exists

      INTEGER(kind=2) :: result1,ntest

      CHARACTER(10) :: coordsys
!      INTEGER      :: xcount,ycount
!      REAL         :: xorigin,yorigin,xdelta,ydelta,
      real         :: a66

      real :: ha,fpw,kdn,nratio

!     rev. 9.9.33  Oct.  16/14  - NK: Added checks for files existing for a resume'
      INQUIRE(FILE='resume.txt',EXIST=exists)
      IF(.not.exists)THEN
        print*
        print*,'Expecting to find a file called'
        print*,'resume.txt'
        print*,'in the working directory but it is not found'
        print*
        stop 'Program aborted in rdresume @ 54'
      endif
      open(unit=99,file='resume.txt',status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,' Problems opening the resume.txt file'
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        print*
        stop ' Program aborted in rdresume @ ~625'
      endif

      read(99,*,iostat=ios)junk1,n_version
      if(ios.ne.0.or.n_version.ne.7)then
        print*,'*********************************************'
        print*
        print*,'Resume file format has changed'
        print*,'You must create a new resume file to be '
        print*,'compatible with this program version'
        print*,'Version found =',n_version,' Ver. 5 required'
        print*
        print*,'Resume now required the use of the' 
        print*,'soil_init.r2c and flow_init.r2c files'
        print*,'These are automatically created with'
        print*,'resumflg=`y'
        print*
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
        print*,'NEW:      Jan. 21, 2019'
        print*,'A new array has been added to the flow_init.r2c file'
        print*,'to make a smoother transition for a resume'
        print*,'Please rerun your spinup period to create new '
        print*,'resume files'
        print*
        print*,'*********************************************'
        stop 'Program aborted in rd_resume @ 80'
      endif
      
      read(99,*)          ! read start date  ver. 5
      read(99,*)          ! read end date    ver. 5

      read(99,*,iostat=ios)classcount1,na1,ycount1,xcount1,no1,noresv1
      if(ios.ne.0)then
        write(98,*)'Error'
        write(98,*)'Resume file format has changed'
        write(98,*)'You must create a new resume file to be '
        write(98,*)'compatible with this program version'
        write(98,*)'Program aborted in rd_resume @ 90'
        print*,'Resume file format has changed'
        print*,'You must create a new resume file to be '
        print*,'compatible with this program version'
        print*
        stop 'Program aborted in rd_resume @ 90'
      endif
      
!     rev. 9.8.87  Oct.  25/13  - NK: Added error message for mismatched resume file
      if(no1.ne.no)then
        write(98,*)'Error'
        write(98,*)'No of flow stations in the str file is different'
        write(98,*)'before & after the resume'
        write(98,*)'Program aborted in rdresume @ line 91'
        print*
        print*,'No of flow stations in the str file is different'
        print*,'before & after the resume'
        print*
        stop 'Program aborted in rdresume @ line 132'
      endif
      if(noresv1.ne.noresv)then
        print*
        print*,'No of flow stations in the rin file is different'
        print*,'before & after the resume'
        print*
        stop 'Program aborted in rdresume @ line 98'
      endif

      if(.NOT.allocated(area_temp))then         ! added Dec. 29/07 nk
        allocate(area_temp(no1),stat=iAllocate)
        if(iAllocate.ne.0) STOP 'Error with allocation of area_temp() in rdresume'
      endif


!     rev. 10.2.61 Aug.  28/19  - NK removed inbsnflg from read_resume
      read(99,*)(junk1,n=1,no1+noresv1)
!c      read(99,*,iostat=ios)(inbsnflg(n),n=1,no1+noresv1)
!      print*,inbsnflg
!    pause 55555

    
    
      if(ios.ne.0)then
        print*,'problems reading the first block'
        print*,' in the resume.txt file'
        print*
        stop ' Program aborted in sub @ 787'
      endif
      if(classcount1.ne.classcount.or.na1.ne.na.or.ycount1.ne.ycount.or.xcount1.ne.xcount)then
        print*,' resume file array size not compatible'
        print*,' with present watershed (.shd) file'
        print*,' '
        write(*,'(A)',advance='no') 'hit enter to continue - at your own risk'
        read(*,*)
      endif
      
!d      print*,'meta data read'
!d      print*,'Reading block 1'
      read(99,*,iostat=ios)block,junk    !block 1
!     print*,block
!     rev. 9.2.40  Jun.  09/06  - NK: added tto(),ttomin(),ttomax() to resume
      do n=1,na1
        read(99,*,iostat=ios)n1,sump(n),qmax(n),sumrff(n)
        if(ios.ne.0)then
          print*,' Problems reading ',block
          print*,' at grid no. ',n,' out of ',na1
          print*,n1,sump(n),qmax(n),sumrff(n)
          print*
          stop 'Program aborted in rdresume @ 833'
        endif
      end do
!d      print*,'Part 1 of block 1 read'
      do n=1,na1
        read(99,*,iostat=ios)n1,(effpor(n,ii),ii=1,classcount)
        if(ios.ne.0)then
          print*,' Problems reading ',block
          print*,' at grid no. ',n,' out of ',na1
          print*,'n,ii,effpor/',n,ii,effpor(n,ii)
          print*
          stop 'Program aborted in rdresume @ 148'
        endif
      end do
!d      print*,'Part 2 of block 1 read'

      read(99,*,iostat=ios)block,junk    !block 2
!d      print*,'Reading block 2'
      do n=1,na1-1
        read(99,*,iostat=ios)tto(n),ttomin(n),ttomax(n)
        if(ios.ne.0)then
          print*,' Problems reading ',block
          print*,' at grid no. ',n,' out of ',na1-1
          print*
          stop 'Program aborted in rdresume @ 146'
        endif
      end do
!     took out qsn(n)  Nov. 7/01  NK
!d      print*,' block ',block,'  read'

!     taken out at one point but needed to make continuous rff files
!     there's no harm at all in having it
      read(99,*,iostat=ios)block,junk    !block 3
!d      print*,'Reading block 3'
      read(99,*,iostat=ios)(sr(ii),x4(ii),sexcess(ii),&
           sq1(ii),sq1fs(ii),sqint(ii),sqintfs(ii),&
           sdrng(ii),sdrngfs(ii),ii=1,classcount1)
!c           sdrng(ii),sdrngfs(ii),ii=1,classcount1+1)
      if(ios.ne.0)then
        write(*,98001)block
        write(98,98001)block
        print*
        STOP ' Program aborted in rdresume @ 169'
      endif
!c      if(iopt.eq.2)print*,' block ',block,'  read'
!d      print*,' block ',block,'  read'

      read(99,*,iostat=ios)block,junk    !block 4
!d      print*,'Reading block 4'
      do n=1,na
        do ii=1,classcount
          read(99,*,iostat=ios)n1,ii1,ssumr(n,ii),sumf(n,ii),sumffs(n,ii)
          if(ios.ne.0)then
            write(*,98001)block
            write(98,98001)block
            write(*,*)'grid no.=',n,' class no.=',ii
            print*
             stop 'Program aborted at 187'
          endif
        end do
      end do

!c      if(iopt.eq.2)print*,' block ',block,'  read'
!d      print*,' block ',block,'  read'
      read(99,*,iostat=ios)block,junk    !block 5
!d      print*,'Reading block 5'
      do n=1,na
        i=yyy(n)
        j=xxx(n)
        read(99,*,iostat=ios)sn1(i,j),basinerr(i,j),nhyd(i,j),nbasin(i,j)
      end do
      if(ios.ne.0)then
        print*,' error in block =',block
             stop 'Program aborted at 201'
      endif

!     rev. 9.8.34  Oct.  23/12  - NK: Added sums to the resume.txt file
!d      print*,'Reading block 7'
      read(99,*,iostat=ios)block,junk    !block 7
      read(99,*)intevt
      read(99,*)evt
      read(99,*)sump
      read(99,*)ssumr
      read(99,*)sumf
      read(99,*)sumffs
      read(99,*)sumrff
      read(99,*)slzinflw
      read(99,*)sdlz
      read(99,*)sum_sublim
      read(99,*)sum_pet
      read(99,*)sum_et



!d      print*,'Reading block 8'
      read(99,*,iostat=ios)block,junk    !block 8
      read(99,*,iostat=ios)&
           (aa(i),bb(i),cc(i),ashnum(i),ashden(i),qpeakh(i),&
           qpeaks(i),nq(i),nqc(i),nxtbasin(i)&
                ,area_temp(i),i=1,no1)

      if(ios.ne.0)then
        write(*,98001)block
        write(98,98001)block
        print*
             stop 'Program aborted at 237'
      endif
!c      if(iopt.eq.2)print*,' block ',block,'  read'
!d      print*,' block ',block,'  read'
      
!     rev. 10.3.03 Jan.  21/20  = NK added dd_ice & dd_thaw to the resume.txt file      
!     For additional dd_ice & dd_thaw data fields 
      read(99,99002)line
99002 format(a256)      
      if(line(1:3).ne.'eof')then
            read(line,*)i,junk
            if(i.eq.9)then
                read(99,*)(dd_ice(n),n=1,na)
                read(99,99002)line
                read(line,*)i,junk
                if(i.eq.10)then
                    read(99,*)(dd_thaw(n),n=1,na)
                endif
                read(line,*)i,junk
                if(i.eq.11)then
                    read(99,*)(ice_fctr(n),n=1,na)
                endif
            else
                read(line,9701,iostat=ios)eofmark
 9705           if(eofmark.ne.'eof')then
                    write(6,9702)
                    STOP 'Program stopped in sub @ eof check'
                endif

                close(unit=99,status='keep')
          
            endif
      else      
          close(unit=99,status='keep')
      endif
      
      



!     rev. 9.2.31  Feb.  09/06  - NK: Added area chaeck to rdresume
      open(unit=99,file='changed_areas.txt',status='unknown',iostat=ios)
      if(ios.ne.0)then    ! added Nov. 10/14  nk
        print*
        print*,'Unable to open file  changed_areas.txt'
        print*,'Possible cause(s):'
        print*,'file in use by another application'
        stop 'Program aborted in rdresume @ 276'
      endif
      change_flag='n'
      do i=1,no1
        if(area_temp(i).ne.area(i))then
          if(change_flag.eq.'n')then
            write(98,98000)'Warning: The drainage area has been changed since the ',&
           'resume file was written likely due to ',&
           'due to re-ordering the flow stations ',&
           'New values are used in this run ',&
           'Program continues but initial flows may take ',&
           'time to settle. ',&
           'Changed_areas have been written ',&
           'in debug\warnings.txt'
98000       format(a9,a45,a38,a37,a32,a45,a16,a32,a21)            
     
            write(98,*)'Warning: Flow Station #, area from resume file, new area computed in flowinit:'

            change_flag='y'
          endif
          write(98,*)'Warning: ',i,'     ',area_temp(i),'     ',area(i)
        endif
      end do
      close(unit=99,status='keep')
      write(98,*)'Info: FINISHED READING THE RESUME.TXT FILE'

      !     REV. 8.61 - Dec.12/97 - ADDED CONTFLG FOR STATS CONTN
      if(contflg.ne.'y')then    
!       RESTART CALC'S FOR THE STATISTICS
        do l=1,no  
          nq(l)=0
          nqc(l)=0
          aa(l)=0.
          bb(l)=0.
          cc(l)=0.
          ashnum(l)=0.0
          ashden(l)=0.0
          qpeakh(l)=0.1e-10
          dpeakh(l)=0.1e-10
          qpeaks(l)=0.0e-10
          dpeaks(l)=0.0e-10
        end do
        do n=1,naa  
          do ii=1,classcount       
            intevt(n,ii)=0.0
            evt(n,ii)=0.0
            sump(n)=0.0
            ssumr(n,ii)=0.0
            sumf(n,ii)=0.0
            sumffs(n,ii)=0.0
            sumrff(n)=0.0
!c        slzinflw(n)=0.0
!c        sdlz(n)=0.0
            sum_sublim(n,ii)=0.0
            sum_pet(n,ii)=0.0
            sum_et(n,ii)=0.0
          end do
        end do
      endif

      if(iopt.ge.1)print*,'Finished reading the resume.txt file' 
      if(iopt.ge.2)pause 'Hit enter to continue'
      
      return

 9701 format(a3)
 9702 format(' resume.txt file length does not match current'/&
             ' memory array spec. Regenerate the file'/)

99000 format(f5.1)
99001 format(f25.0)
99003 format(i1,5x,a50)
99004 format(i5)
!c99182   format(' Warning: Error opening or reading fln:',a30/
!c     *  ' Probable cause: missing strfw/yymmdd.str input file'/
!c     *  ' OR: in config.sys have you set files=100 & buffers=50?'/)
98001 format(' Warning: Error reading resume.txt'/&
     ' Probable cause: wrong format or end of file reached'/&
     ' This could be due to a change in grid characteristics'/&
     ' since creating the resume file'/&
     ' Solution: create a new resume file with current executable'//&
     ' OR, look for *** or NaN in the resume file and change to '/&
     ' 0.0E+00     Look in:'//&
     ' block= ',i5)


      end subroutine rdresume

