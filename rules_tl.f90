      SUBROUTINE rules_tl(n,div,thr,l,jz,at,dtmin,date,time,firstpass)

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
!     rev. 10.2.73 Dec.  14/19  - NK Convert rules_tl.f90 to Fortran 90 & fix drawdown comps.
!***********************************************************************

      use area_watflood
      use area_debug
	  implicit none
      save
	

      Integer  :: ios,nnu,j,k,nrr,i,n,l,ic,jm,jz,ll,nnext,ju_next
	  integer  :: newrel,newrin,iDeallocate,iAllocate
      integer  :: dayrad(12),last_month,ndiv_max,lvl_sta_no
      integer  :: yyyy,mmmm,dddd
      real*4   :: hold,wt,dtmin,at,div,thr,time   !,q_divert,q_fixed
      real*4   :: sup,mhu,stc,eri,ont,mean_elv,delta_elv,temp_elv(100)
      real*4   :: sup_init,mhu_init,stc_init,eri_init,ont_init
      real*4   :: retard_factor(12,5)  ! for great lakes ice-weed retardation
      real*4   :: monthly_evap(12,5),hourly_evap(12,5)    ! for great lakes evap
      real*4   :: qbear,qgold,qkettle,elv,hgold,hbear
      real*4   :: qtake,outflow,last_elv,dlth,ddd,dddlast,ddm,ddmlast
      real*4   :: lif,lif_min_lm,lif_min_lw,log_adj,last_obs,qsum
      real*4   :: datum_LMan,flow,spill,sill,s_factor,store_live
      real*8   :: midrange,midrange1,weight
      real*8   :: drawdown,drawdown_rate,fill_rate
      real*4   :: xrule(999),yrule(999),dif,ua,ub,c,d

      real*4,dimension(:,:),allocatable :: Qhourly,Qoutflow

      
	integer  :: lcount,day_last,lastres,jz_last,nlast,rank,zone
	character(10)  :: yyyymmdd(366),hhmmss(366)
	
	real*4   :: elvlast,qraw

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      character(20) :: junk
      character(30) :: newfilename
      character(10) :: fileformat
	character(14) :: date
      character*1   :: firstpass,chtr(10)
      character*256 :: line
      LOGICAL exists,firstpass_local(99),cedar_lvl,warningflg
      logical foundEndHeader
      
      real*4, dimension(:),   allocatable :: old,raise,lower,inflow
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
        allocate(raise(noresv),lower(noresv),inflow(noresv),stat=iAllocate)
          if(iAllocate.ne.0)STOP  'Error with allocation of raise/lower in rerout @ 121'
        allocate(lower_range(366,noresv),upper_range(366,noresv),&
                old(noresv),&
                full(noresv),max_range(noresv),&
                ruleNo(noresv),resvNo(noresv),stat=iAllocate)
        INQUIRE(FILE='resrl\rules.ts5',EXIST=exists)      

        if(exists)then  
          open(unit=54,file=filename(54),status='unknown',iostat=ios)
          k=1
          do ll=1,noresv
              k=max(k,int(DecayT(ll)))   ! note: k used in a different way later
          end do
          allocate(Qhourly(noresv,k+1), Qoutflow(noresv,k+1),stat=iAllocate)
        
          open(unit=99,file='resrl\rules.ts5', status='old')
          ll=0
          foundEndHeader=.false.
          print*,'NEW<<<<<'
          print*,'Using Coordinates for operating rules:'
          print*
          do while(.not.foundEndHeader)
            read(99,99000)line
99000       format(a256)       
            if(line(1:10).eq.':endHeader')foundEndHeader=.true.   
            if(line(1:10).eq.':EndHeader')foundEndHeader=.true.   
            if(line(1:6).eq.':Point')then
              ll=ll+1
              read(line,*)junk,xrule(ll),yrule(ll)
!c              print*,junk,ll,xrule(ll),yrule(ll)
            endif
          end do
        endif  
        nrules=ll/2   ! # of reservoirs/lakes with rules
        print*,'# reservoirs with rules =',nrules
!c        print*,'       rule#        row         col   rank        resv#'
          
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
          if(DecayT(ll).lt.1.or.DecayT(ll).gt.10000)then
              write(98,*)'Error: - non fatal:'
              write(98,*)'Error: It looks like your level\yyyymmdd_ill.pt2 file'
              write(98,*)'Error: needs to be upgraded to have attribute # 7 '
              write(98,*)'Error: DcayT  as the last column'
              write(*,*)'Error: - non fatal:'
              write(*,*)'Error: It looks like your level\yyyymmdd_ill.pt2 file'
              write(*,*)'Error: needs to be upgraded to have attribute # 7 '
              write(*,*)'Error: DcayT  as the last column'
          endif
        end do     
        
        if(debug_output)write(63,*)'        x             y         ll       rank      res(rank)'
        do ll=1,nrules
!         convert to local coordinate unit system for new .res file
          j=int((xrule(ll*2)-xorigin)/xdelta)+1
          i=int((yrule(ll*2)-yorigin)/ydelta)+1
          rank=s(i,j)
          if(debug_output)write(63,*)xrule(ll*2),yrule(ll*2),ll,rank,res(rank)
!         find the resv # for this rule curve
          if(rank.gt.0)then
!           this reservoir has a rule          
            resRuleFlg(res(rank))=.true.
            ruleNo(res(rank))=ll  !res(rank)
            resvNo(ll)=res(rank)
          else  
            warningflg=.true.
            print*,'Rule coordinates do not match a reservoir location'
            print*,'for rule No.',ll 
          endif
        end do
        
        open(unit=987,file='rule_locations.csv',status='unknown')
        do ll=1,noresv
            rank=s(ires(ll),jres(ll))
            if(resruleFlg(res(rank)))then
                write(987,*)xres(ll),yres(ll),ll,resname(ll)
            endif
        end do
        close(unit=987,status='keep',iostat=ios)
        open(unit=987,file='rule_no_locations.csv',status='unknown')
        do ll=1,noresv
            rank=s(ires(ll),jres(ll))
            if(.not.resruleFlg(res(rank)))then
                write(987,*)xres(ll),yres(ll),ll,resname(ll)
            endif
        end do
        close(unit=987,status='keep',iostat=ios)
      
        if(warningflg)then
          print*
          print*,'Please check that coordinates in the ts5 file'
          print*,'match those in the rel files'
!c          stop 'Program aborted in rule @ 165'
        endif         
        
        do j=1,365
          read(99,*,iostat=ios)yyyymmdd(j),hhmmss(j),(lower_range(j,i),upper_range(j,i),i=1,nrules)
          if(ios.ne.0)then
              print*,'Failed reading the resrl\rules.ts5 file'
              print*,'last read data line j'
              print*,yyyymmdd(j)
          print*,j,yyyymmdd(j),(lower_range(j,i),upper_range(j,i),i=1,noresv)
              stop 'Program aborted in read_tl.f @ 224'
          endif
        end do

!c          read(yyyymmdd(1),*)(chtr(j),j=1,10)  
          print*,yyyymmdd(1)
          if(yyyymmdd(1)(5:5).ne.'-')then
              print*,'Wrong date format'
              print*,'You can edit the rules.ts5 file in Exel but it '
              print*,'must be opened as a tab & space delimited file'
              Print*,'Before saving, change the date format to'
              print*,'2015-01-01'
              print*,'The program will read / as a space and give'
              print*,'an error later as the data is not read'
              stop 'Program aborted in rules.f @ 234'
          endif
!       fill in missing values
!       lower range
       
        Do ll=1,nrules
          nlast=1
!         find 1st >0 value
          j=2
          do while(nlast.lt.365)
            do while(lower_range(j,ll).lt.0.0)                
              j=j+1            
            end do
            nnext=j
!           found a +ve value so interpolate missing for days          
            dif=lower_range(nnext,ll)-lower_range(nlast,ll)
            do j=nlast,nnext
              lower_range(j,ll)=lower_range(nlast,ll)+dif*float(j-nlast)/(nnext-nlast)
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
          do while(nlast.lt.365)
            do while(upper_range(j,ll).lt.0.0)                
              j=j+1            
            end do
            nnext=j
!           found a +ve value so interpolate missing for days          
            dif=upper_range(nnext,ll)-upper_range(nlast,ll)
            do j=nlast,nnext
              upper_range(j,ll)=upper_range(nlast,ll)+dif*float(j-nlast)/(nnext-nlast)
!c              write(769,*)j,nlast+1,n-1,dif,upper_range(j,ll)
            end do
            nlast=nnext
          end do
        end do
!       end upper range        
        close(unit=99,status='keep')
        
!       write the rules to a file - daily delta t        
        open(unit=99,file='results\rules_echo.txt',status='unknown')
        do i=1,ni
          do j=1,365
            write(99,99001)j,(lower_range(j,ll),upper_range(j,ll),ll=1,nrules)
99001       format(i5,<nrules*2>f10.3)     
          end do
!     rev. 10.2.43 Jan.  17/19  - NK: Fixed bug in leapyear extra day in rules_echo.txt
          if(mod(year_now+i-1,4).eq.0)then
!               write and extra line for leap years  (approx)              
                write(99,99001)366,(lower_range(365,ll),upper_range(365,ll),ll=1,nrules)
          endif  
        end do
        close(unit=99,status='keep')

        do ll=1,nrules
          max_range(ll)=-999.0
        end do
        do ll=1,nrules
          do j=1,366
            max_range(ll)=max(max_range(ll),upper_range(j,ll))
          end do
        end do
        close(unit=99,status='keep')
        
        outflow=0.0

        firstpass_local(99)=.false.
        jz_last=-999
!        do l=1,noresv
!            old(l)=qo2(n)
!        end do
      endif   !  end of local first pass to set up the rule table
!~~~~~~~~~~~~~~~~~~~~~~~~~~END FIRST PASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
!     note that this is the global first pass
      if(firstpass.eq.'y')then
!         initialize storage	   
!         storage = live storage
        if(ruleNo(l).gt.0)then 
          temp_elv(l)=b6(l)
          lake_elv(l,jz)=b6(l)    ! initial lake level
          store2(n)=(lake_elv(l,jz)-b7(l))*lake_area(l)
          lower_range(366,l)=lower_range(365,ruleNo(l))
          upper_range(366,l)=upper_range(365,ruleNo(l))
!         ********************************************************weight
!         ********************************************************weight
!         ********************************************************weight
          raise(l)=0.5   ! was 0.75 but too slow filling
          lower(l)=1.5
          weight=1.0/decayT(l)      ! initial weight on last flow
!         ********************************************************weight
!         ********************************************************weight
!         ********************************************************weight
        endif
        if(Qmin(1).lt.0.0)then
            print*,'Error:'
            print*,'Found *_ill.pt2 file but'
            print*,'Qmin and safe_max values' 
            print*,'not found in yyyymmdd_ill.pt2 file'
            print*,'Please fix'
            print*
            stop 'Program abortd in rules @ 277'
        endif
      endif
      
!c      if(jul_day_now.lt.122.or.jul_day_now.gt.288)then
!c                call rerout(n,div,thr,l,jz,at,dtmin,date,time,firstpass)
!c               return
!      endif
!      if(jul_day_now.eq.123)print*,'Switched to rules'
!      if(jul_day_now.eq.289)print*,'Switched to power functions'
      
    
      inflow(l)=qi2(n)
      if(jul_day_now.eq.1)full(l)=.false.

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(ruleNo(l).gt.0)then
!       use only the upper range values to decide if drawdown rate is fast enough
!       drawdown is calculated using the downward slope of the upper range  
        ju_next=jul_day_now+1
!       OK to use 365 as target for 366 - target for 365        
        if(ju_next.ge.365)ju_next=1
!       drawdawn is for the midrange trendline        
!     rev. 10.2.37 Oct.  18/18  - NK: Changed target levels to real*8 
        if(jul_day_now.le.365)then
!         for 366 just leave as is            
          drawdown=upper_range(jul_day_now,ruleNo(l))&
                     -upper_range(ju_next,ruleNo(l))&
                     +lower_range(jul_day_now,ruleNo(l))&
                     -lower_range(ju_next,ruleNo(l))
          drawdown_rate=drawdown*lake_area(l)/43200.        ! 86400.0)/2.0 = 43200
          fill_rate=-drawdown_rate
        endif
      
!       midrange is used as the target
        midrange=(upper_range(jul_day_now,ruleNo(l))&
                   +lower_range(jul_day_now,ruleNo(l)))/2.0
!          write(840+l,84000)jul_day_now,l,ruleNo(l),
!     *          upper_range(jul_day_now,ruleNo(l)),
!     *          lower_range(jul_day_now,ruleNo(l)),
!     *          midrange,drawdown,qmin(l),safe_max(l)
84000 format(3i10,7f12.3)     

      ua=upper_range(jul_day_now,ruleNo(l))-midrange
      ub=TEMP_ELV(L)-midrange
      c=abs(midrange-TEMP_ELV(L))
      d=upper_range(jul_day_now,ruleNo(l))-midrange
        
      if(TEMP_ELV(L).gt.safe_max(l))then
          qo2(n)=1.05*qi2(n)
          weight=1/decayT(l)
          qo2(n)=weight*qo2(n)+(1.0-weight)*qold(n)
      elseif(abs(fill_rate).le.0.0100)then
!         keep level constant          
          IF(TEMP_ELV(L).gt.upper_range(jul_day_now,ruleNo(l)))then
!             above upper target level so lower regardless 
              if(drawdown.gt.0.0)then
!                  qo2(n)=amax1(1.25*qi2(n),1.25*drawdown_rate)  
                  qo2(n)=max(1.1*qi2(n),1.1*drawdown_rate) 
              else
!                  qo2(n)=1.25*qi2(n)  
                  qo2(n)=1.05*qi2(n)  
              endif
              zone=1
          elseIF(TEMP_ELV(L).gt.midrange)then  

              qo2(n)=(c/d+1)*qi2(n)
!              if(zone.lt.2)qold(n)=qo2(n)   ! kill memory
!              zone=2
        
          elseIF(TEMP_ELV(L).gt.lower_range(jul_day_now,ruleNo(l)))then
        
              qo2(n)=(1.0-c/d)*qi2(n)
!              if(zone.lt.3)qold(n)=qo2(n)   ! kill memory
!              zone=3
              
          elseIF(TEMP_ELV(L).le.lower_range(jul_day_now,ruleNo(l)))then
!             below the lower target level                  
              if(temp_elv(l).le.b7(l))then
!                 below the sill                  
                  qo2(n)=0.0
              else
!                  qo2(n)=0.1*qi2(n)
                  qo2(n)=0.01*qi2(n)
              endif
              zone=4
          endif
          weight=1/decayT(l)
          qo2(n)=weight*qo2(n)+(1.0-weight)*qold(n)
      elseif(fill_rate.gt.0.100)then
!         Filling
          IF(TEMP_ELV(L).gt.upper_range(jul_day_now,ruleNo(l)))then
!             above upper target level so lower regardless 
              qo2(n)=1.1*qi2(n)  
              zone=1
           elseIF(TEMP_ELV(L).gt.midrange)then  
        
              qo2(n)=(c/d+1)*qi2(n)-fill_rate
!               if(zone.lt.2)qold(n)=qo2(n)   ! kill memory
             zone=2
        
           elseIF(TEMP_ELV(L).gt.lower_range(jul_day_now,ruleNo(l)))then
        
              qo2(n)=(1.0-c/d)*qi2(n)-fill_rate
!              if(zone.lt.3)qold(n)=qo2(n)   ! kill memory
!              zone=3
              
        
          elseIF(TEMP_ELV(L).le.lower_range(jul_day_now,ruleNo(l)))then
!             below the lower target level                  
              if(temp_elv(l).le.b7(l))then
!                 below the sill                  
                  qo2(n)=0.0
              else
                  qo2(n)=0.1*qi2(n)
              endif
              zone=4
          endif
          weight=1/decayT(l)
          qo2(n)=weight*qo2(n)+(1.0-weight)*qold(n)
      elseif(drawdown_rate.gt.0.0100)then
!         lowering
        !             above upper target level so lower regardless 
          IF(TEMP_ELV(L).gt.upper_range(jul_day_now,ruleNo(l)))then
              if(drawdown.gt.0.0)then
                  qo2(n)=max(1.25*qi2(n),1.25*drawdown_rate)  
              else
                  qo2(n)=1.25*qi2(n)  
              endif
              zone=1
           elseIF(TEMP_ELV(L).gt.midrange)then  
        
              qo2(n)=(1.0+c/d)*drawdown_rate+qi2(n)    !<<<<<<<<<<<<<<<<<<<<<<<<<
              if(zone.lt.2)qold(n)=qo2(n)   ! kill memory
              zone=2
        
           elseIF(TEMP_ELV(L).gt.lower_range(jul_day_now,ruleNo(l)))then
        
              qo2(n)=(c/d-1.0)*drawdown_rate+qi2(n)
              if(zone.lt.3)qold(n)=qo2(n)   ! kill memory
              zone=3
              
           elseIF(TEMP_ELV(L).le.lower_range(jul_day_now,ruleNo(l)))then
!             below the lower target level                  
              if(temp_elv(l).le.b7(l))then
!                 below the sill                  
                  qo2(n)=0.0
              else
                  qo2(n)=0.1*qi2(n)
              endif
              zone=4
          endif
          weight=1/decayT(l)
          qo2(n)=weight*qo2(n)+(1.0-weight)*qold(n)
      else
          
          qo2(n)=qi2(n)    
          
        weight=1/decayT(l)
        qo2(n)=weight*qo2(n)+(1.0-weight)*qold(n)
      endif

!      weight=1/decayT(l)
!      qo2(n)=weight*qo2(n)+(1.0-weight)*qold(n)

      
!      if(l.eq.1.and.iopt99)then
!          write(666,66600)jul_day_now,jz,&
!          drawdown,drawdown_rate,fill_rate,qi2(n),qo2(n),&
!          upper_range(jul_day_now,ruleNo(l)),&
!                     upper_range(ju_next,ruleNo(l)),&
!                     lower_range(jul_day_now,ruleNo(l)),&
!                     lower_range(ju_next,ruleNo(l)),midrange     
!66600     format(2i5,f15.8,99f10.3)          
!      endif
        
        if(qo2(n).lt.0.)qo2(n)=0.0
        
        if(lake_release(l,day_now).gt.0.0)then
              qo2(n)=lake_release(l,day_now)
        endif
        
        qold(n)=qo2(n)
        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
	    lake_inflow(l,jz)=qi2(n)
        
        
        if(nbsflg.eq.'y')then
          net_lake_inflow(l,jz)=qi2(n) ! in rulestl
        else
          net_lake_inflow(l,jz)=qi2(n)-outflow ! in rulestl
        endif
        
        
        lake_elv(l,jz)=b7(l)+(store2(n)-store_dead(l))/lake_area(l)
        temp_elv(l)=lake_elv(l,jz)
      
        if(temp_elv(l).ge.max_range(l))full(l)=.true.
        
      endif
        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      jz_last=jz    

!     rev. 9.5.11  Feb.  12/08  - NK: added -ve storage check for reservoirs
!     Fixed for -ve outflow Mar. 4/08 -nk-
!     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso mosed
        resvstore2(n)=store2(n)

         if(store2(n).le.0.0.and.iopt.ge.1)then
           Print*,'store2(',n,' ) -ve / needs work in rules @ 345'
           write(54,54000)year_now,time,l,store1(n),store2(n),qold(n),&
                qi1(n),qi2(n),qo1(n),qo2(n),div,&
                inflow(l),raise(l),drawdown_rate,drawdown
54000      format(i5,f8.0,i5,3E15.6,10f10.3)           
!c           print*,'Try increasing lake depth '
!c           print*,'check debug\debug\lake_error.txt'
!c           pause 'in rules @ 365'
         endif

	last_month=month_now

  999 RETURN


      END SUBROUTINE rules_tl

