      SUBROUTINE rules_sl(n,div,thr,l,jz,at,dtmin,date,time,firstpass)

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
!   THIS S/R  reads target levels at fixed dates and interpolates the 
!   target levels for each day of the year.
!
!     rev. 9.9.65  Apr.  `3/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
!     rev. 9.9.73  Aug.  31/15  - NK: Finshed rules s/r - ready for beta testing
!                                     Fixed for leapyears Jan. 16 NK
!     rev. 10.2.13 Jan.  31/18  - NK: Re-wrote rules.f to mimic stop log operations -> rules_sl.f
!     rev. 10.4.25 Jul.  25/20  = NK Changed rules_sl to use hourly rule intervales
!***********************************************************************

      use area_watflood
	implicit none
      save
	
c !DEC$ ATTRIBUTES DLLIMPORT :: rules_MH

      Integer  :: ios,nnu,j,k,nrr,i,n,l,ic,jm,jz,iz,ll,nnext,ju_next
	integer  :: newrel,newrin,iDeallocate,iAllocate
      integer  :: dayrad(12),last_month,ndiv_max,lvl_sta_no
      real*4   :: hold,wt,dtmin,at,div,thr,time   !,q_divert,q_fixed
      real*4   :: sup,mhu,stc,eri,ont,mean_elv,delta_elv,temp_elv(100)
      real*4   :: sup_init,mhu_init,stc_init,eri_init,ont_init
      real*4   :: retard_factor(12,5)  ! for great lakes ice-weed retardation
      real*4   :: monthly_evap(12,5),hourly_evap(12,5)    ! for great lakes evap
      real*4   :: qbear,qgold,qkettle,elv,hgold,hbear
      real*4   :: qtake,outflow,last_elv,dlth,ddd,dddlast,ddm,ddmlast
      real*4   :: lif,lif_min_lm,lif_min_lw,log_adj,last_obs,qsum
      real*4   :: datum_LMan,flow,spill,sill,s_factor,store_live
      real*4   :: midrange,midrange1,weight
      real*4   :: drawdown,drawdown_rate,fill_rate
      real*4   :: xrule(999),yrule(999),dif
      real*4   :: Qhourly(50,168),Qoutflow(50,169)
      real*4   :: datum_local,store_temp

	integer  :: lcount,day_last,lastres,nlast,rank,hour_next
	character(10)  :: yyyymmdd(366),hhmmss(8784)
	character(10)  :: yyyymmddhh(8784)
	
	real*4   :: elvlast,qraw,qtemp
c!     rev. 9.9.16  Jun.  06/14  - NK: Added location file for Root R. diversion
c	integer  :: divX(3),divY(3),divGridNo(3)
c	real*4   :: divlon(3),divlat(3),qtweak
c	character*20 :: cjunk(20)

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      character(20) :: junk
      character(30) :: newfilename
      character(10) :: fileformat
	character(14) :: date
      character*1   :: firstpass
      character*256 :: line
      LOGICAL exists,firstpass_local(99),cedar_lvl,warningflg
      logical foundEndHeader
      
      real*4, dimension(:),   allocatable :: old,raise,lower,inflow
      real*4, dimension(:,:),   allocatable :: qqold
      logical, dimension(:),  allocatable :: full

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

	DATA firstpass_local/99*.true./
	DATA day_last/0/
	DATA warningflg/.false./

      nnu=0

!     index = 1 for first pass each new chained event
!     index = 2 for subsequent passes. set in sub

!     rev. 9.1.11  Feb.  07/02  - fixed bug in reservoir routing 

      store1(n)=store2(n)   !  moved from below 'if'  09/11/04 nk

   
!~~~~~~~~~~~~~~~~~~~~~~~~~~START FIRST PASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(firstpass_local(99))then
        allocate(raise(noresv),lower(noresv),inflow(noresv),
     *                       qqold(noresv,9999),stat=iAllocate)
          if(iAllocate.ne.0)STOP
     *      'Error with allocation of raise/lower in rerout @ 121'
        allocate(lower_range(366*24,noresv),upper_range(366*24,noresv),
     *           old(noresv),
     *           full(noresv),max_range(noresv),
     *           ruleNo(noresv),resvNo(noresv),stat=iAllocate)
        
        INQUIRE(FILE='resrl\rules.ts5',EXIST=exists)      

        if(exists)then  
          open(unit=99,file='resrl\rules.ts5',
     *                            status='old')
c        read(99,*)
c        do j=1,365
c          read(99,*)(lower_range(j,i),upper_range(j,i),i=1,noresv)
c        end do
          ll=0
          foundEndHeader=.false.
c          print*,'NEW<<<<<'
c          print*,'Coordinates for operating rules:'
          do while(.not.foundEndHeader)
            read(99,99000)line
c            print*,line(1:72)
99000       format(a256)       
            if(line(1:10).eq.':endHeader')foundEndHeader=.true.   
            if(line(1:10).eq.':EndHeader')foundEndHeader=.true.   
            if(line(1:6).eq.':Point')then
              ll=ll+1
              read(line,*)junk,xrule(ll),yrule(ll)
c              print*,junk,ll,xrule(ll),yrule(ll)
            endif
          end do
        endif  
        nrules=ll/2   ! # of reservoirs/lakes with rules
c        print*,'# reservoirs with rules =',nrules
c        print*,'       rule#        row         col   rank        resv#'
          
!       check the locations in the resrl\rules.ts5 file
!       coincide with the reservoir locations in the rel file
!       and mark the lake/resv as having a rule: 
        do ll=1,noresv
!         These were initially set to .true. in read_resv_ef.f   
!         to be sure rules.f is called at least once to look for the rules 
!         if they exist.     
          resRuleFlg(ll)=.false.
          ruleNo(ll)=-1
          resvNo(ll)=-1
        end do     
        
c      print*,
c     *'         ll       rank      res(rank)    resvNo(res(rank))'
        write(53,*)
        write(53,*)'Reservoirs with rules:'
        do ll=1,nrules
!         convert to local coordinate unit system for new .res file
          j=int((xrule(ll*2)-xorigin)/xdelta)+1
          i=int((yrule(ll*2)-yorigin)/ydelta)+1
          rank=s(i,j)
!         find the resv # for this rule curve
          if(rank.gt.0)then
!           this reservoir has a rule          
            resRuleFlg(res(rank))=.true.
            ruleNo(res(rank))=ll  !res(rank)
            resvNo(ll)=res(rank)
            write(53,*)'Rule #',ll,' grid # ',rank,' res # ',res(rank)
          else  
            warningflg=.true.
            print*,'Rule coordinates do not match a reservoir location'
            print*,'for rule No.',ll 
          endif
        end do
      print*  
      
        if(warningflg)then
          print*
          print*,'Please check that coordinates in the ts5 file'
          print*,'match those in the rel files'
          stop 'Program aborted in rule @ 165'
        endif         
        
        do j=1,8784
            do i=1,nrules
                lower_range(j,i)=-1.
                upper_range(j,i)=-1.
            end do
        end do
c        do j=1,365
        do j=1,8737,24
          read(99,*)yyyymmddhh(j),hhmmss(j),
     *                    (lower_range(j,i),upper_range(j,i),i=1,nrules)
c          print*,j,yyyymmddhh(j),lower_range(j,14),upper_range(j,14)
        end do
        close(unit=99,status='keep')
!       Add a day for leap years:
        do i=1,nrules
              lower_range(8784,i)=lower_range(8737,i)
              upper_range(8784,i)=upper_range(8737,i)
        end do
          
!       fill in missing values
!       lower range
       
        Do ll=1,nrules
          nlast=1
!         find 1st >0 value
          j=2
          do while(nlast.lt.8784)
            do while(lower_range(j,ll).lt.0.0)                
              j=j+1            
            end do
            nnext=j
!           found a +ve value so interpolate missing for days - hours now          
            dif=lower_range(nnext,ll)-lower_range(nlast,ll)
            do j=nlast,nnext
              lower_range(j,ll)=lower_range(nlast,ll)+
     *                         dif*float(j-nlast)/(nnext-nlast)
            end do
            nlast=nnext
          end do
        end do
!       end lower range        

!       upper range
        Do ll=1,nrules
          nlast=1
!         find 1st >0 value
          j=2
c          do while(nlast.lt.365)
          do while(nlast.lt.8784)
            do while(upper_range(j,ll).lt.0.0)                
              j=j+1            
            end do
            nnext=j
!           found a +ve value so interpolate missing for days          
            dif=upper_range(nnext,ll)-upper_range(nlast,ll)
            do j=nlast,nnext
              upper_range(j,ll)=upper_range(nlast,ll)+
     *                         dif*float(j-nlast)/(nnext-nlast)
            end do
            nlast=nnext
          end do
        end do
!       end upper range        
        
!       write the rules to a file - hourly delta t        
        open(unit=99,file='debug\rules_echo.txt',
     *                            status='unknown')
        do i=1,ni
          do j=1,8766,deltat_report
            write(99,99001)j,
     *            (lower_range(j,ll),upper_range(j,ll),ll=1,nrules)
99001       format(i5,<nrules*2>f10.5)     
          end do
        end do
        close(unit=99,status='keep')

        do ll=1,nrules
          max_range(l)=-999.0
        end do
        do ll=1,nrules
c          do j=1,365
          do j=1,8784
            max_range(ll)=max(max_range(ll),upper_range(j,ll))
          end do
        end do
        close(unit=99,status='keep')
        
        outflow=0.0

        firstpass_local(99)=.false.
!       open the output file:        
        open(unit=987,file='results\res_levels.csv',status='unknown')
        
        
c        if(iopt99)then        
c          write(840+l,84001)'jul_day_now','l','res(n)','ruleNo(l)',
c    *          'temp_elv(l)','upper_range','lower_range',
c    *          'midrange','datum_local','store2(n)','qo2(n)',
c    *          'lake_elv(l,jz)','datum_local'  
c84001     format(4a10,8a15)          
c        endif
        
      endif   !  end of local first pass to set up the rule table
!~~~~~~~~~~~~~~~~~~~~~~~~~~END FIRST PASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~END FIRST PASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~END FIRST PASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
!     note that this is the global first pass
      if(firstpass.eq.'y')then
!         initialize storage	   
!         storage = live storage
        if(ruleNo(l).gt.0)then 
          temp_elv(l)=b6(l)
          lake_elv(l,jz)=b6(l)    ! initial lake level
          if(hour_now.eq.24)then
             j=jul_day_now*24-48+hour_now
          else
             j=jul_day_now*24-24+hour_now
          endif
          midrange=(upper_range(j,ruleNo(l))
     *              +lower_range(j,ruleNo(l)))/2.0
          
          store_temp=(temp_elv(l)-midrange)*lake_area(l)
          if(store_temp.gt.0.0)then
            qo2(n)=b1(l)*(store_temp)**b2(l)
          else
            qo2(n)=0.0
          endif
          qo1(n)=qo2(n)
          old(l)=qo2(n)
          weight=.1      ! initial weight on last flow
          
          weight=1.0/decayT(l)      ! initial weight on last flow
          
          if(decayT(l).gt.100)then
              write(98,*)'Warning: decayT for lake # :',l
              write(98,*)'Warning: seems too large'
          endif    
          
        endif
        
            do j=1,decayT(l)
                qqold(l,j)=qo1(n)
            end do
            qo2(n)=qo1(n)
c            if(l.eq.60)then
c                write(333,33300)l,n,qo2(n),(qqold(l,j),
c     *                                       j=1,decayT(l))
33300           format(2i5,999f10.3)   
c                write(333,*)l,decayT(l)
c                write(*,*)l,decayT(l),qo1(n)
c            endif
      endif
!     ~~~~~~~~~~~~~~~~~~~END FIRST PASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      inflow(l)=qi2(n)
c      if(jul_day_now.eq.1)full(l)=.false.
      if(hour_now.eq.1)full(l)=.false.

      if(ruleNo(l).gt.0)then
!       use only the upper range values to decide if drawdown rate is fast enough
!       drawdown is calculated using the downward slope of the upper range  
        ju_next=jul_day_now+1
        hour_next=hour_now+1
c        if(ju_next.gt.366)ju_next=1

!       OK to use 365 as target for 366 - target for 365        
        if(mod(year_now,4).eq.0)then
          if(hour_next.gt.8784)hour_next=1
        else
          if(ju_next.gt.8760)hour_next=1
        endif
      
        if(hour_now.eq.24)then
           j=jul_day_now*24-48+hour_now
        else
           j=jul_day_now*24-24+hour_now
        endif
        J=MAX(1,J)  ! this is a cludge - j = 0 for last hour of the year
        
!       Between mid & lower range is used as the target
        midrange=lower_range(j,ruleNo(l))+
     *   (upper_range(j,ruleNo(l))
     *              -lower_range(j,ruleNo(l)))/4

        datum_local=midrange
        store_temp=(temp_elv(l)-datum_local)*lake_area(l)
        
c        old(l)=qo2(n)
!       move each flow back one iteration        
        do j=int(decayT(l)),1,-1
            qqold(l,j+1)=qqold(l,j)
        end do
        
        if(store_temp.gt.0.0)then
              qtemp=b1(l)*(store_temp)**b2(l)
        else
            qtemp=0.0
        endif
        
c        qo2(n)=weight*qo2(n)+(1.0-weight)*old(l)
!       Use area_debug moving average:        
        qqold(l,1)=qtemp
        qo2(n)=0.0
        do j=1,int(decayT(l))
          qo2(n)=qo2(n)+qqold(l,j)
        end do
c        if(l.eq.50)write(333,*)qo2(n),decayT(l)
        qo2(n)=qo2(n)/decayT(l)
c        if(l.eq.50)write(333,33300)l,n,time,qtemp,qo2(n),
c     *              (qqold(l,j),j=1,20)
        
        if(qo2(n).lt.0.)qo2(n)=0.0
        
        store2(n)=store1(n)+(qi1(n)+qi2(n)     ! moved down
     *                -qo1(n)-qo2(n))*div
	  lake_inflow(l,jz)=qi2(n)
        net_lake_inflow(l,jz)=qi2(n)-outflow
        lake_elv(l,jz)=b7(l)+(store2(n)-store_dead(l))/lake_area(l)
c        lake_elv(l,jz)=datum_local+store2(n)/lake_area(l)
        temp_elv(l)=lake_elv(l,jz)
      
        if(temp_elv(l).ge.max_range(l))full(l)=.true.

      endif  
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
!     rev. 9.5.11  Feb.  12/08  - NK: added -ve storage check for reservoirs
!     Fixed for -ve outflow Mar. 4/08 -nk-
!     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso mosed
        resvstore2(n)=store2(n)

	last_month=month_now

  999 RETURN


      END SUBROUTINE rules_sl

