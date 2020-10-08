      subroutine write_resume(jz)

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
     
!     rev. 9.2.40  Jun.  09/06  - NK: added tto(),ttomin(),ttomax() to resume
!     rev. 10.2.50 Mar.  22/19  - NK: New resume files written to \resume\*.*

      use area_watflood
      use area_debug
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      DIMENSION     :: smc5(16)
      CHARACTER(14) :: date
      CHARACTER(3)  :: eofmark
      CHARACTER(1)  :: lineflg,smok
      character(20) :: junk
      REAL(4)		:: optlow,time,tot1,qwert,conv,scale,
     *               smc5,tj1,clock,t,thr,dtmin,dtmax,div,aintvl,sintvl,
     *               tot2,e1,tdum,qtemp,diff2,sdlz,
     *               wfo_spec_version_number
      INTEGER       :: rbin,inta,block,no1,classcount1,na1,
     *                   iallcnt1,n1,ii1,jan,m,ios,n,iallocate,
     *                   l,ii,juold,jj,lun,nhr,nhf,
     *                   i,j,icase,iz,jz,nchr,mz,ju,mon,
     *                   iwest,ieast,isouth,inorth,nocolumns,nno,k,
     *                   nu,igrdshft,minsnwflg,oldjan,flgevp22,
     *                   ycount1,xcount1,noresv1,nj,npick
      CHARACTER(128):: qstr,line,tmpLine

      INTEGER(kind=2) :: result1,ntest

      CHARACTER(10) :: coordsys
      real         :: a66
      real :: ha,fpw,kdn,nratio
      logical :: exists

      print*,'The to_be_continued flag tbcflg= ',tbcflg
      print*,' Writing a new resume.txt file - in write_resume'

      
!       THIS SECTION IS USED TO WRITE ALL STATE VARIABLES
!       SO PRGRAM CAN BE CONTINUED WHERE IT LEFT OFF
!     rev. 10.2.50 Mar.  22/19  - NK: New resume files written to \resume\*.*
        open(unit=99,file='resume\resume.txt',status='unknown',
     *                              iostat=ios)
        if(ios.ne.0)then
          print*
          print*,'Unable to open file  resume.txt'
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          stop 'Program aborted in write_resume.f @ 46'
        endif
!         we're assuming that noresv>noresvi
        write(99,*)'resume_version= 7'
	  write(99,*)'start y/m/d/h:',
     *              year_start,mo_start,day_start,hour_start
	  write(99,*)'end   y/m/d/h:',
     *              year_now,month_now,day_now,hour_now
        write(99,*,iostat=ios)classcount,na,ycount,xcount,no,noresv,
     *        (inbsnflg(n),n=1,no+noresv)
        block=1
        write(99,*)block,' block'   
!     rev. 9.2.40  Jun.  09/06  - NK: added tto(),ttomin(),ttomax() to resume
        do n=1,na
          write(99,*)n,sump(n),qmax(n),sumrff(n)
	  end do
!       took qsn(n) out  Nov. 5/01  NK

        do n=1,na
          write(99,*)n,(effpor(n,ii),ii=1,classcount)
	  end do

        block=2
        write(99,*)block,' block'   
!     rev. 9.2.40  Jun.  09/06  - NK: added tto(),ttomin(),ttomax() to resume
        do n=1,na-1
          write(99,*)tto(n),ttomin(n),ttomax(n)
        end do


        block=3           !  sums
        write(99,*)block,' block'   
        write(99,*)(sr(ii),x4(ii),sexcess(ii),
     *     sq1(ii),sq1fs(ii),sqint(ii),sqintfs(ii),
     *     sdrng(ii),sdrngfs(ii),ii=1,classcount)

        block=4           !  sums
        write(99,*)block,' block' 
        do n=1,na
          do ii=1,classcount  
            write(99,*)n,ii,ssumr(n,ii),sumf(n,ii),sumffs(n,ii)
          end do
        end do

        block=5
!       write(*,*)block,' block' 
        write(99,*)block,' block' 
        do n=1,na
          i=yyy(n)
          j=xxx(n)
!          write(*,*)sn1(i,j),basinerr(i,j),ssmc(i,j),nhyd(i,j),
!     *                 nbasin(i,j)
          write(99,*)sn1(i,j),basinerr(i,j),nhyd(i,j),
     *               nbasin(i,j)
        end do

        block=6

        do i=1,no
          aa(i)=max(-1.0e+32,aa(i))
          bb(i)=max(-1.0e+32,bb(i))
          cc(i)=max(-1.0e+32,cc(i))
          ashnum(i)=max(-1.0e+32,ashnum(i))
          ashden(i)=max(-1.0e+32,ashden(i))
          qpeakh(i)=max(-1.0e+32,qpeakh(i))
          qpeaks(i)=max(-1.0e+32,qpeaks(i))
          nq(i)=max(-99999,nq(i))
          nqc(i)=max(-1.0e+32,aa(i))
          nxtbasin(i)=max(-99999,nxtbasin(i))
          area(i)=max(-1.0e+32,area(i))
          aa(i)=min(1.0e+32,aa(i))
          bb(i)=min(1.0e+32,bb(i))
          cc(i)=min(1.0e+32,cc(i))
          ashnum(i)=min(1.0e+32,ashnum(i))
          ashden(i)=min(1.0e+32,ashden(i))
          qpeakh(i)=min(1.0e+32,qpeakh(i))
          qpeaks(i)=min(1.0e+32,qpeaks(i))
          nq(i)=min(99999,nq(i))
          nqc(i)=min(1.0e+32,aa(i))
          nxtbasin(i)=min(99999,nxtbasin(i))
          area(i)=min(1.0e+32,area(i))
        end do

!     rev. 9.8.34  Oct.  23/12  - NK: Added sums to the resume.txt file
        block=7      ! additional sums
        write(99,*)block,' block   SV sums'   
        write(99,*)intevt
        write(99,*)evt
        write(99,*)sump
        write(99,*)ssumr
        write(99,*)sumf
        write(99,*)sumffs
        write(99,*)sumrff
        write(99,*)slzinflw
        write(99,*)sdlz
        write(99,*)sum_sublim
        write(99,*)sum_pet
        write(99,*)sum_et

        block=8      ! carry statistics
        write(99,*)block,' block    stats'   
        write(99,*)(aa(i),bb(i),cc(i),ashnum(i),ashden(i),
     *    qpeakh(i),qpeaks(i),nq(i),nqc(i),nxtbasin(i),area(i),i=1,no)
        
!     rev. 10.3.03 Jan.  21/20  = NK added dd_ice & dd_thaw to the resume.txt file      
        block=9      ! dd_ice
        write(99,*)block,' block    dd_ice'   
        write(99,*)(dd_ice(n),n=1,na)
            
        block=10      ! dd_thaw
        write(99,*)block,' block    dd_thaw'   
        write(99,*)(dd_thaw(n),n=1,na)

        block=11      ! dd_thaw
        write(99,*)block,' block    dd_thaw'   
        write(99,*)(ice_fctr(n),n=1,na)

        eofmark='eof'
        write(99,9701)eofmark
        close(unit=99,status='keep')

        author='write_resume      '

!     rev. 10.2.50 Mar.  22/19  - NK: New resume files written to \resume\*.*
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call write_flowinit()
        call write_soilinit()
!     rev. 9.8.91  Oct.  30/13  - NK: Got rid of lzs_init.r2c - data is in flow_init.r2c already
!       this data is in flow_init.r2c so this is redundant.
c        call write_lzsinit()
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!     rev. 9.9.11  Mar.  20/14  - NK: Added lake_level_init.pt2 file for a resume
!     rev. 10.2.50 Mar.  22/19  - NK: New resume files written to \resume\*.*

        open(unit=99,file='resume\lake_level_init.pt2',status='unknown')
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file  lake_level_init.pt2'
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          stop 'Program aborted in sub.f @ 4451'
        endif

        write(99,99011)'########################################   '
        write(99,99011)':FileType pt2  ASCII  EnSim 1.0            '
        write(99,99011)'#                                          '
        write(99,99011)'# DataType               EnSim PT2 Set     '  
        write(99,99011)'#                                          '
        write(99,99011)':Application             WATFLOOD          '
        write(99,99011)':Version                 2.1.23            '   
        write(99,99011)':WrittenBy          spl                    ' 
        write(99,99011)':CreationDate       2014-03-25  07:27      '
        write(99,99011)'#                                          '
        write(99,99011)'#--------------------------------------    '
        write(99,99011)':SourceFile                   flow_data    '
        write(99,99011)'#                                          '
        write(99,99011)':Name                 InitLakeLevels       '
        write(99,99011)'#                                          '
	  write(99,3004)':Projection         ',coordsys1
	  if(coordsys1.eq.'UTM       ')then
          write(99,3004)':Zone               ',zone1
	    write(99,3004)':Ellipsoid          ',datum1
	  elseif(coordsys1(1:7).eq.'LatLong')then
	    write(99,3004)':Ellipsoid          ',datum1
        endif
        write(99,99011)'#                                          '
        write(99,99011)':SampleTime        0000/00/00  00:00:00.0  '
        write(99,99011)'#                                          '    
        write(99,99011)':AttributeName 1     StationName           '  
        write(99,99011)':AttributeType 1     text                  '    
        write(99,99011)':AttributeName 2     InitialElevation      '   
        write(99,99011)':AttributeType 2     float                 '    
        write(99,99011)':AttributeName 3     Datum                 '   
        write(99,99011)':AttributeType 3     float                 '   
        write(99,99011)':AttributeName 4     Depth                 '   
        write(99,99011)':AttributeType 4     float                 '   
        write(99,99011)':AttributeName 5     SafeMax               '   
        write(99,99011)':AttributeType 5     float                 '   
        write(99,99011)':AttributeName 6     Qmin                  '   
        write(99,99011)':AttributeType 6     float                 '   
        write(99,99011)':AttributeName 6     DecayT                  '   
        write(99,99011)':AttributeType 6     float                 '   
        write(99,99011)':endHeader                                 '    
99011   format(a43)                               
        do l=1,noresv    
          write(99,99010)xres(l),yres(l),resname(l)(1:7),
     *      lake_elv(l,jz),b7(l),b7(l)-LKinvert(l),safe_max(l),qmin(l),
     *      DecayT(l)
99010     format(2f20.6,5x,a7,6f12.3)        
        end do     
      
        close(unit=99,status='keep')
 
      go to 77777
        
!     rev. 10.2.42 Jan.  16/19  - NK: Added resume\yyyy-mm-dd folder with resume files
!       Create the resume file directory and file name 
        line='mkdir resume'
        CALL execute_command_line(line(1:12))
!       Create the sub dir for this date      
        if(month_now.le.9.and.day_now.le.9)then
          write(tmpLine,10001)year_now,month_now,day_now
10001     format('mkdir resume\',i4,'-0',i1,'-0',i1)   
        elseif(month_now.le.9.and.day_now.gt.9)then
          write(tmpLine,10002)year_now,month_now,day_now
10002     format('mkdir resume\',i4,'-0',i1,'-',i2)   
        elseif(month_now.gt.9.and.day_now.le.9)then
          write(tmpLine,10003)year_now,month_now,day_now
10003     format('mkdir resume\',i4,'-',i2,'-0',i1)   
        else
          write(tmpLine,10004)year_now,month_now,day_now
10004     format('mkdir resume\',i4,'-',i2,'-',i2) 
        endif
        CALL execute_command_line(tmpLine(1:23))
      
!       make the copy commands 
        write(line,10005)tmpLine(7:23)
10005   format('copy resume.txt  ',a23,'\resume.txt')
        if(debug_output)write(63,*)line(1:51)
        CALL execute_command_line(line(1:51))
        
        write(line,10006)tmpLine(7:23)
10006   format('copy lake_level_init.pt2  ',a23,'\lake_level_init.pt2')
        if(debug_output)write(63,*)line(1:90)
        CALL execute_command_line(line(1:90))
        
        write(line,10007)tmpLine(7:23)
10007   format('copy flow_init.r2c  ',a23,'\flow_init.r2c')
        if(debug_output)write(63,*)line(1:72)
        CALL execute_command_line(line(1:72))
        
        write(line,10008)tmpLine(7:23)
10008   format('copy soil_init.r2c  ',a23,'\soil_init.r2c')
        if(debug_output)write(63,*)line(1:72)
        CALL execute_command_line(line(1:72))
        
77777   continue   
        
        write(98,*)
     *   'Info: A new  resume\resume.txt          file has been written'
        write(98,*)
     *   'Info: A new  resume\flow_init.r2c       file has been written'
        write(98,*)
     *   'Info: A new  resume\soil_init.r2c       file has been written'
        write(98,*)
     *   'Info: A new  resume\lake_level_init.pt2 file has been written'
        write(98,*)
     *   'Info: in the resume directory overwriting the old file'
        write(98,*)'Info: To NOT overwrite, set the' 
        write(98,*)'Info: tbcflg = n in the LAST(!) event file.'
        
      RETURN

3004  format(a20,a10,2x,a10)
9701  format(a3)
 9702 format(' resume.txt file length does not match current'/
     *       ' memory array spec. Regenerate the file'/)

99000 format(f5.1)
99001 format(f25.0)
99003 format(i1,5x,a50)
99004 format(i5)
98001 format(' Warning: Error reading resume.txt'/
     *' Probable cause: wrong format or end of file reached'/
     *' This could be due to a change in grid characteristics'/
     *' since creating the resume file'/
     *' Solution: create a new resume file with current executable'/
     *' block= ',i5)
      end subroutine write_resume

