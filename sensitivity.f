      subroutine sensitivity
     *    (jan,smc5,conv,scale,icase,smok,optlow,igrdshft)

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

!     This subroutine carries out a sensitivity analysis of all
!     optimizable parameters


      use area_watflood
      USE EF_module
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      CHARACTER(128):: qstr
      CHARACTER(72) :: junk
!      CHARACTER(14) :: date
      CHARACTER(1)  :: smok,rvrflg,cvrflg,debug_flag 
      character(20) :: rvrname(12),classname(24)
      CHARACTER*14    :: date
      INTEGER    :: iallcnt,icnt(5),ndir(5),nchr,ix,ios,icase,
     *              iallocate,igrdshft,iyshiftmin,iyshiftmax,
     *              jxshiftmin,jxshiftmax,ishift,jshift,inum,jnum,
     *              l,iw1,iw2,iv,iflg,i,n,ii,j,nhr,nhf,jan
      REAL(4)    ::   smc5(16),errold(5),err(5),chng(5),best(5),hmax
      REAL(4)    :: optlow,e1,scale,ddtemp,cc1,cc2,crit,conv,best1
      real(4)    :: opt_base,up_factor,down_factor,delta_factor
      integer*2  :: result1,ntest

!      DATA ntest/-1358/qstr/'options'/nchr/7/
      DATA ntest/43814/qstr/'optionsss'/nchr/9/
      DATA iallcnt/0/


      allocate(nrvr_array(nrvr,12),class_array(classcount,24),
     *            stat=iAllocate)
      if(iAllocate.ne.0) STOP
     *      'Error with allocation of area8a arrays in optionsa'

!     initialize array values:
      do ii=1,nrvr
        do i=1,12
          nrvr_array(ii,i)=0.0
        end do
      end do
      do ii=1,classcount
        do i=1,24
          class_array(ii,i)=0.0
        end do
      end do
      optim=0.0

      ensimflg='n'
      sensitivityflg=.true.  ! default = false - set in options

      print*
      print*,'Do you want sensitivities on the routing parameters? y/n'
      read*,rvrflg
      print*
      print*,'Do you wnat sensitivities on the hydrol. parameters? y/n'
      read*,cvrflg
      print*
	Print*,'Please enter the % delta you would like to use:'
	print*,'10% is not a bad value'
      read*,delta_factor
	print*

	if(rvrflg.eq.'y'.or.rvrflg.eq.'n')then
	  continue
	else
	  print*,'rvrflg=',rvrflg
	  print*,'Improper response given for routing question'
	  stop 'Program aborted in sensitivities @ 73'
	endif
	if(cvrflg.eq.'y'.or.cvrflg.eq.'n')then
	  continue
	else
	  print*,'cvrflg=',cvrflg
	  print*,'Improper response given for hydrol. question'
	  stop 'Program aborted in sensitivities @ 80'
	endif
	if(delta_factor.lt.0.0.or.delta_factor.gt.50)then
	  print*,'WARNING: delta_factor not between 0 & 50'
	  print*,'Enter your domain at your peril'
	  pause 'Hit enter to continue ^C to quit'
	endif
	
	debug_flag='n'
	if(iopt.ge.1)then
	  print*,'iopt = ',iopt,' is debug mode'
	  print*,'Do you want to run in this way y/n ?'
	  read*,debug_flag
	  if(debug_flag.eq.'n')then
	    debugflg=.false.
	    print*,'iopt set to 0'
	  else
	    print*,'OK. It will be messy!'
	  endif
	else
	  debugflg=.false.
	endif

      print*,'OK, thank you'

      rvrname(1)='flz'
      rvrname(3)='pwr'
      rvrname(5)='r2n'
      rvrname(7)='theta'
      rvrname(9)='kcond'
      rvrname(11)='rlake'
      classname(1)='rec'
      classname(3)='ak'
      classname(5)='akfs'
      classname(7)='retn'
      classname(9)='ak2'
      classname(11)='ak2fs'
      classname(13)='r3'
      classname(15)='mf'
      classname(17)='base'
      classname(19)='fratio'
      classname(21)='not_used'
      classname(23)='sublim_rate'

       
      id=1
      fln(99)='event/event.evt'
      call read_evt(date,conv,scale,smc5,nhr,nhf)
      ensimflg='n'
	call read_shed_ef(31,1)	
      call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
      opt_base=optim
      optim=0.0
      print*,'base value =',opt_base
      print*,'errflg=',errflg
      print*,'-----------------------'

c      id=1
c      fln(99)='event/event.evt'
c      call read_evt(date,conv,scale,smc5,nhr,nhf)
c      ensimflg='n'
c	call read_shed_ef(31,1)	
c      call read_par_parser(32,2)    
c      call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
c      opt_base=optim
c      optim=0.0
c      print*,'base value =',opt_base
c      print*,'-----------------------'
      
      up_factor=1.0+(delta_factor/100.)
      down_factor=1.0-(delta_factor/100.)

      if(rvrflg.eq.'y')then
      print*,'flz:'
      do ii=1,nrvr
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        flz_o(ii)=flz_o(ii)*down_factor
        do n=1,naa
          flz(n)=flz_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
90111   format(a13,i3,a2,i3,a2,2f15.5)     
        nrvr_array(ii,1)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
	  if(debug_flag.eq.'n')iopt=0
        do n=1,naa
          flz(n)=flz_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,2)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'pwr:'
      do ii=1,nrvr
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        pwr_o(ii)=pwr_o(ii)*down_factor
         do n=1,naa
          pwr(n)=pwr_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
     
        nrvr_array(ii,3)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        pwr_o(ii)=pwr_o(ii)*up_factor
         do n=1,naa
          pwr(n)=pwr_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,4)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'r2n:'
      do ii=1,nrvr
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        r2n_o(ii)=r2n_o(ii)*down_factor
        do n=1,naa
          r2n(n)=r2n_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,5)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        r2n_o(ii)=r2n_o(ii)*up_factor
         do n=1,naa
          r2n(n)=r2n_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,6)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'theta:'
      do ii=1,nrvr
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        theta_o(ii)=theta_o(ii)*down_factor
        do n=1,naa
          theta(n)=theta_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,7)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        theta_o(ii)=theta_o(ii)*up_factor
        do n=1,naa
          theta(n)=theta_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,8)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'kcond:'
      do ii=1,nrvr
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        kcond_o(ii)=kcond_o(ii)*down_factor
        do n=1,naa
          kcond(n)=kcond_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,9)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        kcond_o(ii)=kcond_o(ii)*up_factor
        do n=1,naa
          kcond(n)=kcond_o(ibn(n))
        end do
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        nrvr_array(ii,10)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'rlake:'
      do ii=1,nrvr
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
	  if(rlake_o(ii).gt.0.0)then
          rlake_o(ii)=rlake_o(ii)*down_factor
          do n=1,naa
            rlake(n)=rlake_o(ibn(n))
          end do
          call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
          write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
          nrvr_array(ii,11)=(opt_base-optim)/opt_base*100.0
        else
          nrvr_array(ii,11)=-9.999
        endif
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
        ensimflg='n'
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)   
	  if(rlake_o(ii).gt.0.0)then
          rlake_o(ii)=rlake_o(ii)*up_factor
          do n=1,naa
            rlake(n)=rlake_o(ibn(n))
          end do
          call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
          write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
          nrvr_array(ii,12)=(opt_base-optim)/opt_base*100.0
        else
          nrvr_array(ii,12)=-9.999
        endif
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)
      endif


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(cvrflg.eq.'y')then

      print*,'rec:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        rec(ii)=rec(ii)*down_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,1)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        rec(ii)=rec(ii)*up_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,2)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'ak:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        ak(ii)=ak(ii)*down_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,3)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        ak(ii)=ak(ii)*up_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,4)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'akfs:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        akfs(ii)=akfs(ii)*down_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,5)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        akfs(ii)=akfs(ii)*up_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,6)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'retn:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        retn(ii)=retn(ii)*down_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,7)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        retn(ii)=retn(ii)*up_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,8)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'ak2:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        ak2(ii)=ak2(ii)*down_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,9)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        ak2(ii)=ak2(ii)*up_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,10)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'ak2fs:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        ak2fs(ii)=ak2fs(ii)*down_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,11)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        ak2fs(ii)=ak2fs(ii)*up_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,12)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

      print*,'r3:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        r3(ii)=r3(ii)*down_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,13)=(opt_base-optim)/opt_base*100.0
        optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        r3(ii)=r3(ii)*up_factor
        call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *      (opt_base-optim)/opt_base*100.0,optim
        class_array(ii,14)=(opt_base-optim)/opt_base*100.0
        optim=0.0
      end do
      call write_sens_result(rvrname,classname,delta_factor)

	print*,'fm:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        fm(ii)=fm(ii)*down_factor
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	  class_array(ii,15)=(opt_base-optim)/opt_base*100.0
	  optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        fm(ii)=fm(ii)*up_factor
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	  class_array(ii,16)=(opt_base-optim)/opt_base*100.0
	  optim=0.0
	end do
	call write_sens_result(rvrname,classname,delta_factor)

	print*,'base:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        base(ii)=base(ii)-1
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	  class_array(ii,17)=(opt_base-optim)/opt_base*100.0
	  optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        base(ii)=base(ii)+1
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	  class_array(ii,18)=(opt_base-optim)/opt_base*100.0
	  optim=0.0
	end do
	call write_sens_result(rvrname,classname,delta_factor)

	print*,'fratio:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        fratio(ii)=fratio(ii)*down_factor
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	  class_array(ii,19)=(opt_base-optim)/opt_base*100.0
	  optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
        fratio(ii)=fratio(ii)*up_factor
	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	  class_array(ii,20)=(opt_base-optim)/opt_base*100.0
	  optim=0.0
	end do
	call write_sens_result(rvrname,classname,delta_factor)

	print*,'not_used:'
c      do ii=1,classcount
c        id=1
c        fln(99)='event/event.evt'
c        call read_evt(date,conv,scale,smc5,nhr,nhf)
c	  call read_shed_ef(31,1)	
c        call read_par_parser(32,2)    
c        ensimflg='n'
c!       fpetmo is calculated in read_par and is prop'l to ftall
c!       so we can't operated on ftall here
c        do i=1,12
c          fpetmo(i,ii)=fpetmo(i,ii)*down_factor
c        end do
c	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
c        write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
c     *	  (opt_base-optim)/opt_base*100.0,optim
c	  class_array(ii,21)=(opt_base-optim)/opt_base*100.0
c	  optim=0.0
c
c        id=1
c        fln(99)='event/event.evt'
c        call read_evt(date,conv,scale,smc5,nhr,nhf)
c	  call read_shed_ef(31,1)	
c        call read_par_parser(32,2)    
c        ensimflg='n'
c!       fpetmo is calculated in read_par and is prop'l to ftall
c!       so we can't operated on ftall here
c        do i=1,12
c          fpetmo(i,ii)=fpetmo(i,ii)*up_factor
c        end do
c	  call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
c        write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
c     *	  (opt_base-optim)/opt_base*100.0,optim
c	  class_array(ii,22)=(opt_base-optim)/opt_base*100.0
c	  optim=0.0
c	end do
	call write_sens_result(rvrname,classname,delta_factor)

	print*,'sublim_rate:'
      do ii=1,classcount
        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
	  if(sublim_rate(ii).gt.0.0)then
          sublim_rate(ii)=sublim_rate(ii)*down_factor
	    call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
          write(*,90111)'sensitivity -',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	    class_array(ii,23)=(opt_base-optim)/opt_base*100.0
        else
          class_array(ii,23)=-9.999
	  endif
	  optim=0.0

        id=1
        fln(99)='event/event.evt'
        call read_evt(date,conv,scale,smc5,nhr,nhf)
	  call read_shed_ef(31,1)	
        call read_par_parser(32,2)    
        ensimflg='n'
	  if(sublim_rate(ii).gt.0.0)then
          sublim_rate(ii)=sublim_rate(ii)*up_factor
	    call sub(jan,1.0,smc5,conv,scale,icase,smok,optlow,igrdshft)
          write(*,90111)'sensitivity +',int(delta_factor),'%(',ii,')=',
     *	  (opt_base-optim)/opt_base*100.0,optim
	    class_array(ii,24)=(opt_base-optim)/opt_base*100.0
        else
          class_array(ii,24)=-9.999
        endif
	  optim=0.0
	end do
	call write_sens_result(rvrname,classname,delta_factor)
	endif

	return

	end subroutine sensitivity

!******************************************************************

