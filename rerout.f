      SUBROUTINE rerout(n,div,thr,l,jz,at,dtmin,date,time,firstpass)

!DEC$ ATTRIBUTES DLLIMPORT:: rules_MH

!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen 
        
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
     

!     This s/r is meant to be user friendly so operating rules
!     for specific applicatins can be coded.

!***********************************************************************
!   THIS S/R  ROUTES WATER THROUGH A RESERVOIRFq
!   or SIMPLY PASSES ENTERED FLOWS TO THE RIVER
!
!     REV. 7.80   Oct.  29/96 -  SPL7 ADDED YYMMDD.RIN FOR RES INFLOWS
!                             -  UNIT = 39   FLN = 09
!     rev  9.1.03  July  24/01  - added polinomial to reservoir routing
!     rev. 8.99l  Oct.    2001-     fixed reservoir release timing
!     rev. 8.99n  Dec. 31/2001-     fixed nat. res initial flow (JW)
!     rev. 9.1.11  Feb.  07/02  - fixed bug in reservoir routing 
!     rev. 9.1.13  Mar.  23/02  - fixed resv. timing, moved to beginning of dt
!     rev. 9.1.56  Jun.  18/04  - NK: write new rel & rin files to resrl\newfmt folder.
!     rev. 9.1.59  Jul.  15/04  - NK: split rerout into two parts: rdresv & rerout
!     rev. 9.2.39  May.  09/06  - NK: thr added to route & rerout arg list
!     rev. 9.4.11  Jun.  22/07  - NK: reordered rerout for glake 
!     rev. 9.5.11  Feb.  12/08  - NK: added -ve storage check for reservoirs
c     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso mosed
!     rev. 9.5.32  Jun.  04/08  - NK: compute reservoir levels
!     rev. 9.5.38  Oct.  14/08  - NK: added optional coef6 & 7 to rel file for lake levels
!     rev. 9.5.51  Jan.  13/09  - NK: added reading yyyymmdd_ill.pt2 foa all lakes
!     rev. 9.5.67  Oct.  06/09  - NK: fixed bug in rerout
!     rev. 9.7.22. Mar.  07/11  - NK: Changed diversion code: give/route take/rerout
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
!     rev. 9.8.56  Apr.  10/13  - NK: Added check in rerout for -ve storage due to evaporation
!     rev. 9.8.58  Apr.  12/13  - NK: REvised Family Lake (WPEGR) O/R in rerout
!     rev. 9.8.65  May   28/13  - NK: Dimensioned firstpass_local()in REROUT
!     rev. 9.9.16  Jun.  06/14  - NK: Added location file for Root R. diversion
!***********************************************************************

      use area_watflood
	implicit none
	
c !DEC$ ATTRIBUTES DLLIMPORT :: rules_MH
	

      Integer  :: ios,nnu,j,k,nrr,i,n,l,ic,jm,jz
	integer  :: newrel,newrin,iDeallocate,iAllocate
      integer  :: dayrad(12),last_month,ndiv_max,lvl_sta_no
      real*4   :: old,hold,wt,dtmin,at,div,thr,time   !,q_divert,q_fixed
      real*4   :: sup,mhu,stc,eri,ont,mean_elv,delta_elv,temp_elv(100)
      real*4   :: sup_init,mhu_init,stc_init,eri_init,ont_init
      real*4   :: retard_factor(12,5)  ! for great lakes ice-weed retardation
      real*4   :: monthly_evap(12,5),hourly_evap(12,5)    ! for great lakes evap
      real*4   :: qbear,qkettle,elv,hgold,hbear,qgold
      real*4   :: qtake,outflow,last_elv,dlth,ddd,dddlast,ddm,ddmlast
      real*4   :: lif,lif_min_lm,lif_min_lw,log_adj,last_obs,qsum
      real*4   :: datum_LMan,flow,spill,sill,s_factor,store_live
      real*4   :: zp,kkk   ! coefficients for Lake Athabasca routing rule
      real*4   :: fudge,Qfudge,Qlast,WSL   ! coefficients for GSL routing rule
      integer  :: PPx,PPy,PPn    !  needed for Lake Athabasca routing
	integer  :: lcount,day_last
	real*4   :: elvlast,qraw,wollaston_base_store,qtemp
      integer  :: gold_grid_no
c!     rev. 9.9.16  Jun.  06/14  - NK: Added location file for Root R. diversion
c	integer  :: divX(3),divY(3),divGridNo(3)
c	real*4   :: divlon(3),divlat(3),qtweak
c	character*20 :: cjunk(20)

!     for Jenpeg 
      integer  :: jenpeg_grid_no,jenpeg_div_no,jenpeg_gauge_no
!     For the Lake Winnipeg model
      integer  :: FamL_S_grid_no
      integer  :: ivalue
      real(4)  :: xvalue,yvalue,parValue(10)

!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
      character(20) :: junk
      character(30) :: newfilename
      character(10) :: fileformat
	character(14) :: date
      character*1   :: firstpass
      character*256 :: line
      LOGICAL exists,firstpass_local(99),cedar_lvl,warningflg
      logical       :: natural
      
      real*4, dimension(:),   allocatable :: raise,lower,inflow
      real*4, dimension(:),   allocatable :: qo1_temp,qo2_temp

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

	DATA firstpass_local/99*.true./
	DATA day_last/0/
	DATA warningflg/.true./

!     For the great lakes only:
      DATA retard_factor/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.,0.,
     * 590.,480.,110.,110.,0.0,0.0,0.0,0.0,0.0,0.0,0.0,110.,
     * 650.,510.,110.,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,140.,
     * 110.,140.,80.,140.,0.,60.,230.,140.,80.,60.,0.0,0.0,
     * 170.,300.,150.,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,10./
!     initial water elevations - now read from a file
!	data sup_init,mhu_init,stc_init,eri_init,ont_init/
!     *      183.2,175.98,174.8,174.01,74.61/

!     lif = lake ice factor

!     WHEN LOCAL IS -VE, RESERVOIR RELEASES ARE NEED IN AND THESE 
!     BECOME THE OUTFLOWS OF THE SQUARE IN WHICH THE DAM OUTLET 
!     STRUCTURE IS LOCATED.  
!     WHEN LOCAL = 0 , THERE ARE NO RESERVOIRS, REROUT IS NOT CALLED.

!     WHEN LOCAL IS +VE, THE INFLOW TO THE SQUARE IN WHICH THE
!     RESERVOIR IS LOCATED IS PRINTED OUT ON FILE 11 FOR USE IN
!     HEC-5 OR OTHER OPERATING PROGRAM.

!     THE OVERALL OPERATING SEQUENCE IS THEN THE FIRST RUN SIMPLE
!     TO CALCULATE THE RUNOFF PRODUCED IN VARIOUS SUB-BASINS THEN
!     TO DETERMINE RELEASES, AND THEN TO ROUTE THE RELEASES AND LOCAL
!     INFLOWS USING HYMO.  THE RELEASES ARE READ IN FROM FILE 12.

      nnu=0

!     index = 1 for first pass each new chained event
!     index = 2 for subsequent passes. set in sub

!     rev. 9.1.11  Feb.  07/02  - fixed bug in reservoir routing 

      store1(n)=store2(n)   !  moved from below 'if'  09/11/04 nk
      
      if(firstpass_local(98))then
        last_month=month_now
        wrtdiverflg=.false.  !used for writing a new diversion file 
!       if divertflg = 'g' then diversion flows need to be generated 
!       for L. St. Joseph and NOT read in. In this case, qdivert 
!       has not been allocated in read_divert_ef.f so it needs to be 
!       done here for just one location:
        if(.not.allocated(qdivert))then
!         set one diversion location
          nodivert=1
          allocate(qdivert(nodivert,mhtot),stat=iAllocate)
          if(iAllocate.ne.0)STOP
     *      'Error with allocation of aqdivert  in rerout @ 109'
        endif
        if(iopt.ge.1)print*,'qdivert allocated as      1,',mhtot
        firstpass_local(98)=.false.
      endif
      
      if(firstpass_local(99))then
        allocate(raise(noresv),lower(noresv),inflow(noresv),
     *           qo1_temp(noresv),qo2_temp(noresv),stat=iAllocate)
          if(iAllocate.ne.0)STOP
     *      'Error with allocation of raise/lower in rerout @ 121'
     
!       For lake Athabasca ONLY:
        natural=.false.
        INQUIRE(FILE='natural.txt',EXIST=exists)
        if(exists)then
          open(unit=99,file='natural.txt',status='old',iostat=ios)
          read(99,*)natural
          IF(NATURAL)THEN
            print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            PRINT*,'Using naturalized flow for the MRB'
            print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
          endif
        endif
        firstpass_local(99)=.false.
      endif
      
      if(firstpass_local(97))then
!       FIND THE GRID NUMBER FOR pEACE pOINT         
!     REV. 10.1.40 Oct   11/16  - NK: Fixed bug in read_divert for missing u/s DA
        if(resname(l)(1:7).eq.'Athabas'.and.routeflg.ne.'q')then
          PPy=int((59.118-yorigin)/ydelta)+1
          PPx=int((-112.437-xorigin)/xdelta)+1
          PPn=s(ppy,ppx)
          b7(2)=207.0   !  datum
          firstpass_local(97)=.false.
        endif
      endif

!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
!     rev. 9.5.53  Jan.  20/09  - NK: undid rev. 9.5.40
!##########################################################################
!     DIVERSIONS
!##########################################################################
!     rev. 9.5.40  Oct.  21/08  - NK: added diversions to rerout
!     A diversion can only be made to another node providing
!     the receiving node is at a lower elevation - i.e. the grid number is 
!     higher in the order = later in the grid loop n=1,naa
!     Makes sense. Can not divert to a higher grid without pumping

!     If in an outlet node of a lake from which diversion is made,
!     add an outflow from this reach to the inflow in another. 

!     Lake St Joseph:
!     outflow grid = 1803    05QB006
!     receiving grid = 1866
c      if(n.eq.1803.and.na.eq.3039)then
c        q_divert=qhyd(34,(jz-1)/ktr*ktr+ktr)        
c        q_divert=amax1(0.0,q_divert)     ! in case there are missing data (-1.00)
c        qi2(1866)=qi2(1866)+q_divert   ! diversion into receiving grid
c        write(53,*)n,jz,(jz-1)/ktr*ktr+ktr,qi2(1866),q_divert
c      else
c        q_divert=0.000
c      endif
!###########################################################################

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        

!     rev. 9.8.82  Sep.  07/13  - NK: Bypass of hard-coded lake rules when coeff1=0
!     Note:
!     b1(l) can be less or more than 1.00 for the polinimial case. So don't change!
      if(b1(l).ne.0.00000)then  
!      Start of rule-based routing
!      If lake name matches name below, coefficients will be ignored
!      and rules used instead.
          
       if(resname(l).eq.'Superior     ')then

!       Lake Superior
!         this means that:
!         1. we are doing the great lakes
!         2. we are initializing lake superior here
        if(firstpass.eq.'y')then
!         initialize storage	   
!         storage = live storage 
c	    sup=sup_init
	    sup=b6(1)
	    store1(n)=(sup-181.43)*82.1e+09
c	    store1(n)=(sup-b7(1))*82.1e+09
	    store2(n)=store1(n)
          qo2(n)=824.7*(sup-181.43)**1.5-retard_factor(mo1,1)        
          qo1(n)=qo2(n)
	  endif

        if(b1(1).ne.0.0)then
          qo2(n)=824.7*(sup-181.43)**1.5-retard_factor(mo1,1)        
	  else
	     if(qrel(l,jz).gt.0.0)qo2(n)=qrel(l,jz)
	  endif

        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div

        sup=store2(n)/82.1e+09+181.43

        lake_elv(l,jz)=sup
        lake_inflow(l,jz)=qi2(n)
        net_lake_inflow(l,jz)=qi2(n) !Superior
c	print*,sup,mhu,stc,eri,ont
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
      elseif(resname(l)(1:5).eq.'Huron')then

!       Lake Michigan-Huron
        if(firstpass.eq.'y')then
!         initialize storage	   
!         storage = live storage 
c	    mhu=mhu_init
	    mhu=b6(2)
          store1(n)=(mhu-166.98)*117.4e+09
c          store1(n)=(mhu-b7(1))*117.4e+09
	    store2(n)=store1(n)
c	    delta_elv=mhu-stc_init
	    delta_elv=mhu-b6(3)
c	    mean_elv=(mhu_init+stc_init)/2.0
	    mean_elv=(b6(2)+b6(3))/2.0
	  else
!     rev. 9.4.11  Jun.  22/07  - NK: reordered rerout for glake 
	    delta_elv=mhu-stc
	    mean_elv=amax1(166.98+0.1,mean_elv)   ! prevent div by 0
	    delta_elv=amax1(0.001,delta_elv)        ! prevent div by 0
	  endif
!       use stc from the previous time step. Slow change anyway.

        qo2(n)=82.2*(mean_elv-166.98)**1.87*(delta_elv)**0.36
     *           -retard_factor(mo1,2)        
        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
        mhu=store1(n)/117.4e+09+166.98

!        print*,n,resname(l),mhu,qo2(n),store1(n),store2(n)
        lake_elv(l,jz)=mhu
        lake_inflow(l,jz)=qi2(n)
        net_lake_inflow(l,jz)=qi2(n)-lake_outflow(l-1,jz) !Huron
c	print*,sup,mhu,stc,eri,ont
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
      elseif(resname(l)(1:7).eq.'StClair')then

!       Lake St. Clair
        if(firstpass.eq.'y')then
!         initialize storage	   
!         storage = live storage 
c	    stc=stc_init
	    stc=b6(3)
	    store1(n)=(stc-164.91)*1.11e+09
c	    store1(n)=(stc-b7(3))*1.11e+09
	    store2(n)=store1(n)
c          delta_elv=stc_init-eri_init
          delta_elv=b6(3)-b6(4)
	  else
!     rev. 9.4.11  Jun.  22/07  - NK: reordered rerout for glake 
          delta_elv=stc-eri
	    delta_elv=amax1(0.001,delta_elv)        ! prevent div by 0
	  endif
        stc=amax1(0.1,stc)                    ! prevent div by 0
!       use eri from the previous time step. Slow change anyway.

        qo2(n)=28.8*(stc-164.91)**2.28*(delta_elv)**0.305        
     *           -retard_factor(mo1,3)        

        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
        stc=store1(n)/1.11e+09+164.91

!        print*,n,resname(l),stc,qo2(n),store1(n),store2(n)
        lake_elv(l,jz)=stc
        lake_inflow(l,jz)=qi2(n)
        net_lake_inflow(l,jz)=qi2(n)-lake_outflow(l-1,jz) ! StClair
c	print*,sup,mhu,stc,eri,ont
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
      elseif(resname(l)(1:4).eq.'Erie')then

!       Lake Erie
        if(firstpass.eq.'y')then
!         initialize storage	   
!         storage = live storage 
c	    eri=eri_init
	    eri=b6(4)
	    store1(n)=(eri-169.86)*25.7e+09
c	    store1(n)=(eri-b7(4))*25.7e+09
	    store2(n)=store1(n)
        endif       

        qo2(n)=558.3*(eri-169.86)**1.60-retard_factor(mo1,4)
        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
        eri=store1(n)/25.7e+09+169.86  
        lake_elv(l,jz)=eri
        lake_inflow(l,jz)=qi2(n)
        net_lake_inflow(l,jz)=qi2(n)-lake_outflow(l-1,jz) ! Erie

c	print*,sup,mhu,stc,eri,ont
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
      elseif(resname(l)(1:7).eq.'Ontario')then

!       Lake Ontario
        if(firstpass.eq.'y')then
!         initialize storage	   
!         storage = live storage 
c	    ont=ont_init
	    ont=b6(5)
	    store1(n)=(ont-69.474)*18.96e+09
c	    store1(n)=(ont-b7(5))*18.96e+09
	    store2(n)=store1(n)
	  endif

        qo2(n)=555.823*(ont-0.0014*real(2000-1985)-69.474)**1.5
     *           -retard_factor(mo1,5)        

        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
        ont=store1(n)/18.96e+09+69.474         
!        print*,n,resname(l),ont,qo2(n),store1(n),store2(n)
!        print*,'year=',year
        lake_elv(l,jz)=ont
	  lake_inflow(l,jz)=qi2(n)
        net_lake_inflow(l,jz)=qi2(n)-lake_outflow(l-1,jz)  ! Ontario
        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
      include '..\ManHydro\lake_rules.fi'       ! users:  remove this
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        
       else   !  resname(l).eq.
       
!       NATURAL RESERVOIR ROUTING:
!       NATURAL RESERVOIR ROUTING:
!       NATURAL RESERVOIR ROUTING:
!       NATURAL RESERVOIR ROUTING:
!       NATURAL RESERVOIR ROUTING:
!       rev. 9.7.22. Mar.  07/11  - NK: Changed diversion code: give/route take/rerout


!       For the Mackenzir River:
!       REV. 10.1.14 Jan.  05/16  - NK: Added ice rules for Great Slave.
        if(resname(l)(1:7).eq.'Gr_Slav')then
!         works way better without ice_factor!!  NK Jul. 02/16        
          lake_ice_factor(l,month_now)=
     *            amax1(lake_ice_factor(l,month_now),0.7)
c          write(777,*)year_now,month_now,day_now,lake_ice_factor(l,month_now)
        endif

!       get initial value:        
        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
!     rev. 9.8.56  Apr.  10/13  - NK: Added check in rerout for -ve storage due to evaporation
c        if(store2(n).le.0.0)then
c          write(98,9801)time,l,n,store2(n),qi2(n)
c          store2(n)=1.0
c          if(qi2(n).lt.0.0)qi2(n)=0.0 !can't evaporate water that is not there
c        endif
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
c        old=qo1(n)
        qold(n)=qo1(n)
        hold=0.1e+26
        
        if(b3(l).eq.0.0)then
!         using a power function        
!         tried to put this in the iteration loop but got spikes
!     rev. 9.9.44  Nov.  28/14  - NK: Added dead storage to reservoirs 
          if(store2(n)-store_dead(l).gt.0.0)then
            do while(abs(hold-store2(n)).gt.0.003*hold.and.ic.lt.20)
                if(store2(n)-store_dead(l).gt.0.0)then
!                 have to do at least 3 iterations
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
                  qo2(n)=b1(l)*(store2(n)-store_dead(l))**b2(l)
                  qo2(n)=qo2(n)*lake_ice_factor(l,month_now)
                  wt=amax1(0.5,float(ic)/21.0)
!     rev. 10.2.71 Nov.  18/19  - NK Bug fixes in wetland & reservoir routing
c                  qo2(n)=(1.0-wt)*qo2(n)+wt*old
                  qo2(n)=(1.0-wt)*qo2(n)+wt*qold(n)
c                  old=qo2(n)
                  qold(n)=qo2(n)
                  hold=store2(n)
                  store2(n)=store1(n)+
     *                  (qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
                else
                  qo2(n)=0.0
                  store2(n)=store1(n)+
     *                  (qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
                endif
            end do
          else     !(b3(l).gt.0.0)
            qo2(n)=0.0
            store2(n)=store1(n)+
     *                  (qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
!           store2(n) is allowed to go -ve     
          endif
	    lake_inflow(l,jz)=qi2(n)
c        net_lake_inflow(l,jz)=qi2(n)-outflow
        if(nbsflg.eq.'y')then
          net_lake_inflow(l,jz)=qi2(n) ! in rulestl
        else
          net_lake_inflow(l,jz)=qi2(n)-outflow ! in rulestl
        endif

c         lake_elv(l,jz)=-1.0
!     rev. 9.5.32  Jun.  04/08  - NK: compute reservoir levels
!     rev. 9.5.38  Oct.  14/08  - NK: added optional coef6 & 7 to rel file for lake levels
c	    if(b6(l).gt.0.0.and.b7(l).lt.0.000001)then

          lake_elv(l,jz)=b7(l)+(store2(n)-store_dead(l))/lake_area(l)
        
!     rev. 9.8.57  Apr.  12/13  - NK: Added lakeEflg to stop lake evaporation whan levels very low
!         if there is less than 10 mm of water left in a lake, stop the 
!         evaporation so storage can not go -ve & cause routing problems     
!         used in runoff
!         Needs to be done like this so lake evap can be properly considred
!         for te watbal s/r     
          if(lake_elv(l,jz).lt.b7(l))then
            lakeEflg(l)=.false.
          else
            lakeEflg(l)=.true.
          endif
          
c25        continue   ! sorry about that   

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        else       !(b3(l).gt.0.0)
!         using a polinomial        
!     rev. 9.8.56  Apr.  10/13  - NK: Added check in rerout for -ve storage due to evaporation

          if(store2(n).le.0.0)then
            write(98,9801)time,l,n,store2(n),qi2(n)
c            store2(n)=1.0
            if(qi2(n).lt.0.0)qi2(n)=0.0 !can't evaporate water that is not there
          endif

!     rev. 9.9.44  Nov.  28/14  - NK: Added dead storage to reservoirs 
          store_live=store2(n)-store_dead(l) 
          if(store_live.gt.0.0)then
            do while(abs(hold-store2(n)).gt.0.003*hold.and.ic.lt.20)
!             have to do at least 3 iterations
!             rev  9.1.03  July  24/01  - added polinomial
              qo2(n)=store_live*(b1(l)+store_live*(b2(l)+store_live*
     *                    (b3(l)+store_live*(b4(l)+b5(l)*store_live))))
              qo2(n)=qo2(n)*lake_ice_factor(l,month_now)
              qo2(n)=amax1(0.0,qo2(n))
              wt=amax1(0.5,float(ic)/21.0)
              qo2(n)=(1.0-wt)*qo2(n)+wt*old
c              old=qo2(n)
              qold(n)=qo2(n)
              hold=store2(n)
              store2(n)=store1(n)+
     *                  (qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
            end do   ! hold-store
          else     
            qo2(n)=0.0
            store2(n)=store1(n)+
     *                  (qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
!           store2(n) is allowed to go -ve     
          endif

c         if(store2(n).le.0.0)then
c            qo2=0.0
c            store2(n)=1.0
c            write(53,6804)n,l
c          endif
        endif    ! using a polinomial  

	  lake_inflow(l,jz)=qi2(n)
c        net_lake_inflow(l,jz)=qi2(n)-outflow
        if(nbsflg.eq.'y')then
          net_lake_inflow(l,jz)=qi2(n) ! in rulestl
        else
          net_lake_inflow(l,jz)=qi2(n)-outflow ! in rulestl
        endif

c        lake_elv(l,jz)=-1.0
!     rev. 9.5.32  Jun.  04/08  - NK: compute reservoir levels
!     rev. 9.5.38  Oct.  14/08  - NK: added optional coef6 & 7 to rel file for lake levels
c	  if(b6(l).gt.0.0.and.b7(l).lt.0.000001)then

        lake_elv(l,jz)=b7(l)+(store2(n)-store_dead(l))/lake_area(l)
        
!     rev. 9.8.57  Apr.  12/13  - NK: Added lakeEflg to stop lake evaporation whan levels very low
!       if there is less than 10 mm of water left in a lake, stop the 
!       evaporation so storage can not go -ve & cause routing problems     
!       used in runoff
!       Needs to be done like this so lake evap can be properly considred
!       for te watbal s/r     
        if(lake_elv(l,jz).lt.b7(l))then
          lakeEflg(l)=.false.
        else
          lakeEflg(l)=.true.
        endif
          
c        else
c          lake_elv(l,jz)=store2(n)/lake_area(l)
c	  end if

d        if(iopt.ge.2.and.iopt.le.10)then
d          write(53,6004,iostat=ios)n,l,
d    *           qi1(n),qi2(n),store1(n),store2(n),qo1(n),qo2(n)
d          if(ios.ne.0)then
d            print*,'problems for grid #',n,'lake #',l,' @',time
d          endif
d        endif

!       CALCULATE THE DETENTION TIME
!        if(qo1(n).gt.0.001) at=store2(n)/qo1(n)
!         yeah.... fix this:
        if(qo1(n).gt.0.001) at=store2(n)/qo2(n)

!       SELECT MINIMUM TRAVEL TIME FOR THE TIME STEP CALCULATION
c        dtmin=amin1(at,dtmin)
	
!       DTMIN IS THE TIME REQUIRED TO COMPLETELY DRAIN THE FASTEST
!       EMPTYING ELEMENT              


cccc      endif  ! end of rule based lake/reservoir outflow

       endif  !  resname(l).eq.'Superior 
      
!     rev. 9.8.82  Sep.  07/13  - NK: Bypass of hard-coded lake rules when coeff1=0
      else       !   if(b1(l).eq.0.00000)    FROM RELEASE TABLE:

!       From release table
!       From release table
!       From release table
!       From release table

!       rev. 8.99l  Oct.    2001-     fixed reservoir release timing
!         jm=jz+1        old way see JW's e-mail Oct. 23/01

        jm=jz
        if(jz.gt.nrel)jm=nrel

	  if(jm.lt.1)then
          qo2(n)=qrel(l,1)          !+q_divert
	  else
          qo2(n)=qrel(l,jm)          !+q_divert
	  endif
      
        if(qo2(n).lt.0.0)qo2(n)=0.0

!        this line is for the water balance only
!        it doesn't work for releases

c        store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
        store2(n)=
     *     store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
	  lake_inflow(l,jz)=qi2(n)
        if(nbsflg.ne.'y')net_lake_inflow(l,jz)=qi2(n)-qo2(n) ! default Coefficient in rerout

!     rev. 9.5.32  Jun.  04/08  - NK: compute reservoir levels
!     rev. 9.5.38  Oct.  14/08  - NK: added optional coef6 & 7 to rel file for lake levels


!     rev. 9.9.20  Jul.  24/14  - NK: Added dead storage for lakes "stroe_dead"
c          lake_elv(l,jz)=b7(l)+store2(n)/lake_area(l)
          lake_elv(l,jz)=b7(l)+(store2(n)-store_dead(l))/lake_area(l)
          
          

!     rev. 9.5.11  Feb.  12/08  - NK: added -ve storage check for reservoirs
!     Fixed for -ve outflow Mar. 4/08 -nk-
!     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso mosed
        resvstore2(n)=store2(n)


!	  if(frcflg.eq.'y')then
	  if(trcflg.eq.'y')then  ! otherwise upsets tracer code
!     rev. 10.1.80 Apr.  26/17  - NK: Fixed tracer turnoff for -ve resv. storage
c          if(store2(n).lt.cap(n))then
          if(lake_elv(l,jz).le.b7(l)-LKdepth(l))then
c              print*,jz,lake_elv(l,jz),b7(l)-LKdepth(l)
!           leave everything unchanged in the reservoir
!     rev. 9.9.20  Jul.  15/14  - NK: Fix -ve lake storage when release data used
            if(b1(l).eq.0.0)then
!             we're run into a problem where the releases in the rel file are
!             greater than the inflows and the reservoir is already below live storage.
!             We can fix this by having extra depth as specifies in the ill.pt2 file
              trcflg='n'
              print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              print*,'WARNING: tracer turned off!!!!!!!!!!!!!!!!!'
              print*,'Reason:  lake storage has become -ve'
              print*,'in lake # ',l
              if(warningflg)then      ! added Jan. 08/15  NK
                print*,'Possible fix:'
                print*,'Create a initial lake level file with a datum'
                print*,'low enough so storage will not go -ve. '
                print*,'You may have to experiment until you get a '
                print*,'proper datum. The datum is used as the weir '
                print*,'elv for the lake outflow calculation'
                warningflg=.false.
                print*,'Program won`t pause or beep again'
              if(iopt.ge.1)pause 45678
              endif
              print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!             this is fixed by having an yyyymmdd_ill.tb0 file              
            endif

            if(qi2(n).gt.0.0)then
              qo2(n)=(qi1(n)+qi2(n)-qo1(n))/10.0
	        if(qo2(n).lt.0.0) qo2(n)=qi2(n)/2.  ! TS: added to fix -ve outflow problem
	      else
	        qo2(n)=0.0
	      end if
            store2(n)=store1(n)+(qi1(n)+qi2(n)-qo1(n)-qo2(n))*div
	    endif

d	    if(store2(n).le.0.0.and.iopt.ge.1)then
d           Print*,'store2(',n,' ) -ve / needs work in rerout @ 345'
d           print*,store1(n),store2(n),qi1(n),qi2(n),qo1(n),qo2(n)
d         endif
	  endif

c        lake_elv(l,jz)=-1.0
c        lake_elv(l,jz)=store2(n)/lake_area(l)

        if(iopt.ge.2)write(53,6803)l,jm,n,ireach(n),qo2(n)

      endif   !        if(b1(l).ne.0.0)

	last_month=month_now
      
  999 RETURN

! FORMATS

  500 format(256f10.3)
  501 format(3i5,4x,a1)
  502 format(' resv flow data extrapolated ',i5,' hours')
  504 format(' noresv,nrel,ktr/',3i5)
 1011 format(' ',3x,'  i  ires(i) jres(i)    b1(i)     b2(i)',
     *	'    b3(i)     b4(i)')
 1013 format(' ',3x,i3,2i8,5f10.5,a12/)
 4901 format(25i1)
 4902 format(3i5)
 4903 format(a12)
 4904 format(256f10.0)
 4905 format(256f10.3)
 5003 format(2i5,4g10.3,5x,a12)
 3704 format(2i5,5g10.3,5x,a12)

!     rev. 9.1.55  Jun.  12/04  - NK: write new files to resrl\newfmt folder.
 5004 format(a20,a10)
 5005 format(a20,i5)
 5006 format(a20,a1)
 5007 format(a20,256i1)
 5008 format(a20,f12.0)
 5009 format(a20,a2,'-',a2,'-nn',a2)
 5010 format(a20,a2,a4)

 5301 format(' ','Reservoir inflow data echoed:')
 5303 format(6(' ',a12))
 5304 format(' ','Error on unit=99,fln=',a30,'(',i2,')'//)
 5310 format(' -ve flow for reservoir #',i3,'zero flow assumed',i3)
 6004 format('n,l,qi,store,qo/',2i5,2f10.3,2e15.6,2f10.3)
 6005 format(f8.2,2f10.3,2e15.6,2f10.3)
 6801 format('   rerout: reservoir no =',i3,' mhtot =',i5)
 6802 format('   ',256f8.2)
 6803 format(' rerout: l,m,n,ireach(n),qo2(n)/',4i5,f10.2,f12.0)
 6804 format(' warning: store2(',i5,') set = 0.0 for resv no.',i5) 
 9005 format(' iymin,iymax,jxmin,jxmax/',4i5)
9801  format(f10.1,' resv',i3,' grid',i6,
     *         ' store2=',g12.0,'< 0 1.0 assumed')
99182   format(' Warning: Error opening or reading fln:',a30/
     *  ' Probable cause: missing strfw/yymmdd.str input file'/
     *  ' OR: in config.sys have you set files=100 & buffers=50?'/)
9983  format(256(f12.0,x))
9984  format(256(E12.6,x))
9985  format(256I10)
9986  format(256(a12,x))
9987  format(256(f12.7,x))
9988  format(256(a9,x))

      END SUBROUTINE rerout

