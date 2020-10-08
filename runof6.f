      SUBROUTINE runof6(jan,time,t,thr,mon,e1,mz,ju)

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
     
!***********************************************************************
!     runof4 - snowmelt added
!     runof5 - seperate sca and bare areas added
!     runof6 - wetland runoff and routing added
!
!  THIS SUBROUTINE CALCULATES THE SURFACE RUNOFF AND SUBSURFACE
!  FLOW FOR EACH ELEMENT OF THE WATERSHED =qr(n)
!
!     modified for phillips formula ~ 1974
!     modified for multiple classes in each element 1986
!     modified for jd's snow runof4
!     REV. 7.5 seperate snow covered and bare ground
!     modified for separation of snowcovered ground and bare
!              ground by Frank Seglenieks  Feb/1995   runof5
!     REV. 7.75 added ak2fs   - drainage under snow covered ground
!     REV. 7.76   jun.  11/96 - # classes increased to 16 + urban
!     REV. 7.9    dec 18/96   - Added Todd's evaporation
!     REV. 8.1  - Feb.  15/96 - TBC & RSM (to be continued & resume) 
!     REV. 8.5  - Oct.  09/97 - deleted the old interseption stuff
!     REV. 8.51 - Oct.  09/97 - fixed -ve qr() problem in runof5
!     REV. 8.52 - Nov.  14/97 - replaced x4()= in runof
!     REV. 8.60 - Nov.  14/97 - added sl1 to the interflow calculation
!     REV. 8.85 - Oct.  12/98 - fixed rain & snow on water class
!     REV. 8.89 - Nov.  30/98 - simplified uzs parameters
!     REV. 8.99b  Sept. 27/00 - divvy up interflow & drainae
!
!     REV. 9.00    Mar.  2000 - TS: CONVERTED TO FORTRAN 90 
!     REV. 9.03    Nov.  2000 - TS: ADDED WATFLOOD SWAMP ROUTING
!     rev. 9.1     May    7/01  - updated Luis's sed & nutrient stuff
!                                 in runof6 calculate rf(), hsed(), qs()
!     rev  9.1.20  Jun.  25/02  - Added A10 as the power on the UZ discharge function
!     rev. 9.1.37  Mar.  22/03  - Option to turn off leakage by setting LZF < 0.0
!     rev. 9.1.38  Aug.  18/03  - TS: Add qstrm(n) to save value of qstream(n) for Tracer
!     rev. 9.1.48  Dec.  08/03  - NK: sumrechrge() added to get total recharge
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
!     rev. 9.2.07  Jul.  29/05  - NK: soilinit moved from runoff to sub 
!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!     rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
!     rev. 9.2.42  Jun.  20/06  - NK: water class included in the water balance
!     rev. 9.3.10  Jan.  29/07  - NK: routing pars changed to gridded values
!     rev. 9.4.06  May.  09/07  - NK: replaced por with spore(n,ii) in runof6
!     rev. 9.4.08  May.  29/07  - NK: changed baseflow argument list
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
!     rev. 9.5.61  Sep.  03/09  - NK: bug/eloss - added water class for wfo weighted et
!     rev. 9.5.83  Feb.  17/10  - NK: non_basin exclusion for dds_flag=1
!     rev. 9.7.13  Nov.  22/10  - NK: Changed the outfiles.txt for more rff classes
!     rev. 9.8.59  May   14/13  - NK: REmoved psmear & punused from the program
!     rev. 9.8.61  May   22/13  - NK: Introduced flag1 to speed up runof6
!     rev. 9.8.62  May   22/13  - NK: Fixed bug in runof6: (classcount-3) to (classcount-2) 
!     rev. 9.8.67  Jun   06/13  - NK: Added allocation for flag1
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
!     rev. 10.1.51 Nov.  08/16  - NK: removed unused isotope related calculations, merged two 
!     rev. 10.1.51                NK: isotope related calc sections to reduce if statements
!
!     changes made to include c&g model stuff  nk  April. 26/07
!
!          * * * all depths in mm * * *
!  unit 74 is for evt stuff written in aet.for
!  the value for qr (n) is calculeted in this subroutine
!  h is the height of the crop in m.
!  permeability coefficient=ak mm/hr
!  void ratio=vol.voids/vol.solids=e
!  porosity=vol.voids/tot.vol.=por   << no longer use e1

!  aclass(n,classcount) - % of impervious area in element n
!  aclass(n,ii)      - % of class ii in element n
!  p(i,j)            - rainfall in mm.
!  jj : ibn(n)      - basin number (5=max) - used to be iak(ii)
!  intcap(n,ii)      - max. amount of water that can be held as 
!                      intercepted initially v=0.  intcap is satisfied 
!                      before any water reaches the excess then becomes
!                      net precipitation.
!  oldsca(n,ii)      - snow covered area from the previous time step
!                       (calculated)
!  uzsfs(n,ii)       - upper zone storage for snow covered ground
!                       (calculated)
!  d1fs(n,ii)        - d1 for snow covered ground (calculated)
!  r3fs(ii)          - ground roughness for snow covered ground (read
!                      in from parameter file)
!  akfs(ii)          - permeability for snow covered ground (read
!                      in from parameter file)
!  fakefs(ii)        - infiltration capacity for snow covered ground
!                       (calculated)
!  fexcess(n,ii)     - amount of runoff from the snowpack     
!                       (calculated in melt.f)

!  p    - net prec. available for runoff.when we enter this subroutine
!         the rainfall
!  v    - depth of water intercepted
!  thr  - time in hours
!  fake - infiltration in the time step
!  t    - time in seconds
!  rec  - interflow depletion coefficient
!  qlz  - baseflow recession parameter
!  pwr  - is the lz outflow exponent
!  ds   - depression storage maximum value
!  pen  - depth of water penetration into soil
!         see philip paper 1954 #368
!  x2   - a parameter for interception - see linsley
!  step2- is grid length in km squared
!  t    - time step in seconds
!  tdum - converts mm of water depth to cms for the time step
!         so depth in mm is multiplied by tdum, flow is divided by tdum
!  uzs  - depth of water stored below the soil surface (mm)
!  lzs  - depth of water in the saturated zone (mm)
!  fake - infiltration capacity in the current time interval
!  ak   - sat. cond. in mm/sec
!
!***********************************************************************

      use areacg
      use area_watflood
      use area_debug
cc      use omp_lib
      implicit none


!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      character*8  :: tempflg      ! just a dummy variable
      real*4   :: power1,thr,xdum,t,tdum,e1,qlzfrac,a51,zj,time,
     *            dend,duz,demand,supply,fraction,evaptemp,q2,dcheck,
     *            dlz,sdlz,aju,amon,loss,tempv1,tempvmin1,rh1,
     *            eff_bare_area,eff_sc_area,temp_junk

      integer  :: i,jan,classcount_9,ii,ios,n,j,
     *     mon,ju,jj,l,iiwater,mzz,mz,iAllocate

      data tempv1,tempvmin1,rh1/-999.0,-999.0,-999.0/
c      logical  :: firstpass
c      data firstpass/.true./

!     THE GROUNDWATER DEPLETION FUNCTION:
      do n=1,naa
         if(flz(n).le.0.0)then
            print*,'flz for grid river class no',ibn(n),' =',flz(n)
            print*,'In grid #',n
            print*,'in row #',yyy(n)
            print*,'in column #',xxx(n)
            print*,'Value must be larger than 0.0'
            STOP 'Program aborted in runof6 @ 139'
         endif
      end do

!     XDUM CONVERTS 1 mm depth to cu. m. for the centre grid
      xdum=1000.*step2
!     tdum converts 1 mm/hour depth to cms for the centre grid
      tdum=1000.*step2/t

!     note:  jan can't really be replaced by first pass because of various restarts
!            in options: sensitivity, pattern search etc.
      if(jan.eq.1)then
!       firstpass through this routine
        classcount_9=min0(9,classcount)
!       Append to existing rff files
!       FOR THE RESUME CASE, GO TO THE END OF THE output FILEs: 
        if(resumflg.eq.'y'.and.contflg.eq.'y'.and.id.eq.1)then
          if(iopt.ge.1)then
c            do ii=1,classcount_9
            do ii=1,classcount
              ios=0
              do while(ios.eq.0)
                read(80+ii,6040,iostat=ios)tempflg  ! simout\rffii.txt
              end do
            end do
          endif
        endif
!       rev. 9.4.06  May.  09/07  - NK: replaced por with spore(n,ii) in runof6
!       use of por discontinued 
!       rev. 9.8.61  May   22/13  - NK: Introduced flag1 to speed up runof6


c!      rev. 9.8.67  Jun   06/13  - NK: Added allocation for flag1
c	   if(.not.allocated(flag1))then
c	     allocate(flag1(na),stat=iAllocate)
c	     if(iAllocate.ne.0)stop
c     *     'Error with allocation of flag1 array in runof6'
c	    endif
c         do n=1,naa
c           if(wetflg.eq.'y'
c     *         .and.theta(n).gt.0.00001
c     *         .and.nclass(classcount-3)(1:7).ne.'glacier'
c     *         .and.nclass(classcount-4)(1:7).ne.'glacier')then
c!              both of these could be glaciers as we can have 1 or 2 wetland classes     
c             flag1(n)=.true.
c           else
c             flag1(n)=.false.
c           endif
c         end do

c      In flowinit now
c!     REV. 10.1.19 Jan.  15/16  - NK: Fixed initialization of ice_factr - moved from lake_ice > runof6
c!       default ice_fctr        
c        do n=1,naa
c          ice_fctr(n)=1.0
c        end do

!     rev. 10.1.97 Sep   11/17  - NK: Moved hdrflg action in runof6.f
        if(id.eq.1.and.contflg.eq.'n')then 
           if(frcflg.eq.'y')write(95,7402)
!          write the headers in the rff**.txt & resin.csv files           
           if(iopt99)then     ! added Oct. 18/17 NK
             do ii=1,classcount
               if(aclass(nnprint,ii).gt.0.0)then
                open(unit=900+ii,file=filename(900+ii),status='unknown')
                 write(900+ii,6002)
               endif
             end do
           endif  
        endif
      endif

!     * * * * INITIAL VALUES * * * *
!
d     if(iopt.eq.3)print*,' checkpoint 1 in runoff. JAN=',jan

      GO TO(1001,1002,1003)jan

 1001 CONTINUE

!     FIRST TIME THROUGH THE PROGRAM

      sdlz=0.0
      qlzfrac=1.00
      if(numa.eq.0)write(51,5190)qlzfrac

d      if(iopt.eq.3)print*,' checkpoint 1a in runoff - after soilinit'

!     UPDATE THE UNSATURATED ZONE SOIL MOISTURE
!     SO FAR, THERE IS ONLY ONE SOIL MOISTURE FOR ALL CLASSES
!     ALTHOUGH GRANTED, IT'S A WEIRD IDEA FOR SWAMPS, WHERE MAYBE
!     WE OUGHT TO SET IT FOR 100%

! TS: added e-separation file header Apr.5/06
!     rev. 9.8.33  Oct.  23/12  - NK: Deleted header for rff files with resumflg = y
      if(resumflg.eq.'n')then
        if(iopt99)then
c           if(frcflg.eq.'y')write(95,7402)
c           do ii=1,classcount
c             if(aclass(nnprint,ii).gt.0.0)write(900+ii,6002)
c           end do
!       this section added to write state variables at time=0.0
          n=nnprint
            mzz=mz-24*(mz/24)
             classcount_9=min0(9,classcount)
c             amon=float(mon)
             amon=float(mo1)   !  changed jan22/11 nk for sensitivity repeated runs
             aju=float(ju)
!            this section repeated below
             if(snwflg.eq.'y'.or.vapflg.eq.'y')then
!              this is needed because memory is allocated 
!              only for these cases
               tempv1=tempv(n)
               tempvmin1=tempvmin(n)
               rh1=rh(n)
             else
               tempv1=-999.
               tempvmin1=-999.
               rh1=-999.
             endif
c             do ii=1,classcount_9
             do ii=1,classcount
               if(aclass(n,ii).gt.0.0)  
     *                write(900+ii,6000)0.000,intevt(n,ii),evt(n,ii),
     *         p(ipr,jpr),sump(n),ssumr(n,ii),
     *         amin1(fake(ii),1000.),amin1(fakefs(ii),1000.),
!     rev. 10.2.39 Nov.  15/18  - NK: changed snowc(n,ii) to snowc(n,ii)*sca(n,ii) in runof6
     *         sca(n,ii),snowc(n,ii)*sca(n,ii),d1(n,ii),d1fs(n,ii),
     *         sumf(n,ii),sumffs(n,ii),uzs(n,ii),uzsfs(n,ii),lzs(n),
     *         q1(n,ii),q1fs(n,ii),qint(n,ii),qintfs(n,ii),qlz(n),
     *         drng(n,ii),drngfs(n,ii),qr(n),qstream(n),strloss(n),
     *         sumrff(n),fexcess(n,ii),glmelt(n),fmadj(n),
     *         sq1(ii),sq1fs(ii),sqint(ii),sqintfs(ii),
     *         sdrng(ii),sdrngfs(ii),slzinflw,sdlz,amon,aju,def(n,ii),
     *         tempv1,tempvmin1,rh1,      !psmear(n),punused(n),
     *         api(n,ii)*100.0,sublim(n,ii),sum_sublim(n,ii),
     *         v(n,ii),wcl(n,ii),sum_pet(n,ii),sum_et(n,ii),pet(n,ii),
     *         tto(n),ttomin(n)
             end do
        endif
      endif   

      if(numa.eq.0.and.resumflg.eq.'n'.and.iopt99)call watbal(1)

cd      if(iopt.eq.3)print*,' checkpoint 1b in runoff- after WATBAL'

1002  CONTINUE                         ! PRODUCTION MODE

cd      if(iopt.eq.3)print*,' checkpoint 2 in runoff'

      a51=a5**thr
      zj=p(ipr,jpr)

!     CALC RUNOFF FOR EACH SQUARE

      do n=nastart,naend
        i=yyy(n)
        j=xxx(n)
        qr(n)=0.0
        qstream(n)=0.0
        qstrm(n)=0.0
        strloss(n)=0.0
        sump(n)=sump(n)+p(i,j)
        rechrg(n)=0.0
        qdrng(n)=0.0
        qdrngfs(n)=0.0
            
        sumq1(n)=0.0
        sumq1fs(n)=0.0
        sumqint(n)=0.0
        sumqintfs(n)=0.0
!        sumrechrg(n)=0.0
        glmelt(n)=0.0   !added Jun. 3/02  nk
      end do

      if(frcflg.eq.'y')then !   added for iso
        do n=nastart,naend            
          storeGW2(n)=0.0         
        end do
      endif

c????????????????????????????????????????????????????????????
      if(jan.eq.1)then
      do ii=1,classcount
        sr(ii)=0.0
        sqint(ii)=0.0
        sqintfs(ii)=0.0
        sdrng(ii)=0.0
        sdrngfs(ii)=0.0
        sq1(ii)=0.0
        sq1fs(ii)=0.0           
        uzsinit(ii)=0.0
        sexcess(ii)=0.0
      end do
      endif

      leakage=0.0

cd      if(iopt.eq.3)print*,' checkpoint 2a in runoff'

!     FOR FLGEVP2.GE.1.0 EVAPORATION IS FROM PLANTS AND SOIL USING 
!     CLIMATIC DATA USING THE EVAP.DAT TABLE AS THE POTENTIAL RATE

!     FOR FLGEVP2 = 0.0, EVAPORATION IS FROM THE SOIL ONLY AND EQUAL TO 
!     THE VALUES IN THE EVAP.DAT TABLE AS LONG AS MOISTURE IS 
!     AVAILABLE IN THE UZS
cd      if(iopt.eq.2)print*,vapflg,'in runof6'

      if(vapflg.eq.'y'.and.flgevp2.ge.1.0)then
cd         if(iopt.eq.2)print*,' gone to etina'
         call etin(mon,jan,ju)
cd         if(iopt.eq.2)print*,' back from etina & gone to intcept'
         call intcept(mon)
cd         if(iopt.eq.2)print*,' back from intcept & gone to aet'
         call aet(time,t,mon,ju)
cd         if(iopt.eq.2)    print*,' back from aet'
      else
         do n=nastart,naend
            do ii=1,classcount
               i=yyy(n)
               j=xxx(n)
               r(n,ii)=p(i,j)
               ssumr(n,ii)=ssumr(n,ii)+r(n,ii)
            end do
         end do
      endif

cd      if(iopt.eq.3)print*,' checkpoint 2b in runoff'

!     CALL THE SNOWMELT SUBROUTINE IF THERE IS ANY SNOE ON THE GROUND
!     SNOW MELT IS CALCULATED FOR ALL CLASSES FOR ALL ELEMENTS
!           WITH ONE CALL TO THE SNOW MELT S/R.
!     THE S/R IS CALLED AT LEASR ONCE TO SEE IF THERE IS SNOW

      if(snwflg.eq.'y')then
         call melt(time,jan,ju)
      endif

cd      if(iopt.eq.3)print*,' checkpoint 2c in runoff'

!     * * * * * * CALC RUNOFF FOR EACH GRID * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * CALC RUNOFF FOR EACH GRID * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * CALC RUNOFF FOR EACH GRID * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     * * * * * * CALC RUNOFF FOR EACH GRID * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     Grid loop

cc!$OMP PARALLEL DO   NUM_THREADS(2)

      do n=1,naa
      
        jj=ibn(n)

cd     if(iopt.eq.3)print*,'In runof6, passed 101'

      if(slope(n).gt.0.0)then
        i=yyy(n)
        j=xxx(n)
        l=nhyd(i,j)

cd     if(iopt.eq.3)print*,'In runof5, passed 201'

        qlz(n)=0.0
        
        do ii=1,classcount   ! changed Mar. 31/06 nk

cd     if(iopt.eq.3)print*,'In runof5, passed 301'

!         these variables added Jul. 7/06  nk
          eff_bare_area=frac(n)*aclass(n,ii)*(1.0-sca(n,ii))
          eff_sc_area=frac(n)*aclass(n,ii)*sca(n,ii)

!     rev. 9.5.15  Feb.  28/08  - NK: fixed tdum & xdum for proper grid area in lat-long
  
! * * * * * * * *  TS  - WETLAND ROUTING OPTION * * * * * * * * * * * * * 

c!     rev. 9.8.61  May   22/13  - NK: Introduced flag1 to speed up runof6
c          if(ii.eq.classcount-2.and.flag1(n))then


!         rev. 9.8.77  Jul   08/13  - NK: Made universal the use of wetland_flag(n)
          if(ii.eq.classcount-2.and.wetland_flag(n))then
!           rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!                we have to assume that if there is a glacier, 
!                there is no wetland in that grid
!                i.e. we run this code only when NOT in a glacier grid
!                only glacier flow bypasses the wetland
            
            if(sca(n,ii).le.0.001)then
!              NO SCA - ALL PRECIP ON BARE GROUND
               qswrain(n)=r(n,ii)*wetwid(n)*rl(n)/1000.0/t
         
            else
!              THERE IS A FRACTION OF SCA AND BARE GROUND
!              melt is added from the sca and
!              rain is added from the bare area
               qswrain(n)=
     *              (fexcess(n,ii)*sca(n,ii)+r(n,ii)*(1.0-sca(n,ii)))
     *               *wetwid(n)*rl(n)/1000.0/t
            endif

!           TS - NOT USED IN WETLANDS, SET TO ZERO TO AVOID NaN
!           TS - ADDED DRNG(FS) PARAMS (22/08/2006)
!           TS - ADDED DF(FS) PARAMS (02/10/2006)
            d1(n,ii)=0.0
            d1fs(n,ii)=0.0
            uzs(n,ii)=0.0
            uzsfs(n,ii)=0.0
            df(n,ii)=0.0
            dffs(n,ii)=0.0
            drng(n,ii)=0.0
            drngfs(n,ii)=0.0
!           NK - set to zero to prevent write errors fro wetland class.
            q1(n,ii)=0.0
            q1fs(n,ii)=0.0
            qint(n,ii)=0.0
            qintfs(n,ii)=0.0
            qdf(n,ii)=0.0
            qdffs(n,ii)=0.0
            qdrng2(n,ii)=0.0               
            qdrngfs2(n,ii)=0.0
            drng(n,ii)=0.0
            drngfs(n,ii)=0.0
            dprecip(n,ii)=0.0

! * * * * * * * *  TS  - END WETLAND ROUTING OPTION * * * * * * * * * * * * * 

          else   ! non-wetland classes

!     *       PAUSE 'channel runof5 loop for wetland class'

cd     if(iopt.eq.3)print*,'In runof5, passed 401'

            q1(n,ii)=0.0
            q1fs(n,ii)=0.0
            qint(n,ii)=0.0
            qintfs(n,ii)=0.0
            qdf(n,ii)=0.0
            qdffs(n,ii)=0.0
            qdrng2(n,ii)=0.0               
            qdrngfs2(n,ii)=0.0
            drng(n,ii)=0.0
            drngfs(n,ii)=0.0
            dprecip(n,ii)=0.0

          if(frcflg.eq.'y')then
!           TS - ADDED storeSW1(FS) PARAMS (22/08/2006)
!           TS - ADDED storeIF1(FS) PARAMS (02/10/2006)
              storeSW1(n,ii)=storeSW2(n,ii)
              storeSW2(n,ii)=0.0
              storeIF2(n,ii)=0.0
          endif

          if(aclass(n,ii).gt.0.0)then   ! skip if class area = 0.0
          
!           REV. 7.5 SEPERATE SNOW COVERED AND BARE GROUND

!           THE PRECIPITATION NOT INTERCEPTED NOW REACHES THE GROUND
!           WFILE SNOWMELT IS KEPT SEPARATE IN FEXCESS
!           d1(n) IS THE WATER IN SURFACE RETENTION NOW ON PERMEABLE 
!           GROUND
!           sumf(n,ii), sumffs(n,ii) IS THE SUM OF THE INFILTRATION UP 
!           TO THE PRESENT TIME
!           df(n), dffs(n) IS THE NON-CUMMULATIVE INFILTRATION

!   ++++++++++  THIS IS THE MAIN SECTION TO DEAL WITH ++++++++++
!   ++++++++++             BARE GROUND                ++++++++++

cd     if(iopt.eq.3)print*,'In runof5, passed 501'

            if(ak(ii).gt.0.0)then
!             i.e. water is bypassed            
!      +++++++++++++++++++++++++++++++++++++++++++++++++ fs begin
!      ++++++++   RUNOFF FOR BARE GROUND   +++++++++++++
!      +++++++++++++++++++++++++++++++++++++++++++++++++
d     if(iopt.eq.3)print*,'In runof5, passed 600'
              if(sca(n,ii).lt.1.0)then
                d1(n,ii)=d1(n,ii)+r(n,ii)
                dprecip(n,ii)=r(n,ii) 
              endif
              if(nclass(ii)(1:7).eq.'glacier')then
!               NOTE: AUTOMATICALLY ON BARE ICE ONLY
                d1(n,ii)=d1(n,ii)+glmelt(n)
!               ADD GLACIER MELT TO PRECIP FOR WATER BALANCE:
                sump(n)=sump(n)
     *               +glmelt(n)*aclass(n,ii)*(1.0-sca(n,ii))
              endif
cd             if(iopt.eq.3)print*,'In runof5, passed 601'
!             REV. 7.42 CHECK FOR DIVISION BY 0  MAY/95
!              if(uzs(n,ii).gt.0.01)then    !caused discontinuities nk jul 10/06
              if(uzs(n,ii).gt.1.0E-10)then
                fake(ii)=ak(ii)*thr*
     *          (1.+effpor(n,ii)*(d1(n,ii)+pot(ii))/(uzs(n,ii)+0.00001))
!               limit added Dec. 11/00 nk
                fake(ii)=amin1(fake(ii),1000.1)
!               limit added April 19, 2002 AB
                fake(ii)=amax1(fake(ii),0.001)
              else
                fake(ii)=1000.0
              endif
cd             if(iopt.eq.3)print*,'In runof5, passed 701'
c              if(fake(ii).eq.0.0)print*,'ii,fake',ii,fake

!             INFILTRATION:
              if(d1(n,ii).le.fake(ii))then
!               SURFACE STORAGE < INFILTRATION CAPACITY 
!               AND ALL THE WATER IS INFILTRATED-NO RUNOFF

                uzs(n,ii)=uzs(n,ii)+d1(n,ii)
                sumf(n,ii)=sumf(n,ii)+d1(n,ii)*eff_bare_area
                if(frcflg.eq.'y')then
                  qdf(n,ii)=qdf(n,ii)+d1(n,ii)*eff_bare_area*tdum
                endif
!                 prorated for sca Jul. 7/06 nk
                d1(n,ii)=0.0
              else
!               SURFACE STORAGE > INFILTRATION CAPACITY
                uzs(n,ii)=uzs(n,ii)+fake(ii)
!     rev. 10.1.55 Nov.  30/16  - NK: Fixed sumf & sumffs in runof6
c                sumf(n,ii)=sumf(n,ii)+fake(ii)
                sumf(n,ii)=sumf(n,ii)+fake(ii)*eff_bare_area
!               prorated for sca Jul. 7/06 nk
                d1(n,ii)=d1(n,ii)-fake(ii)
                if(frcflg.eq.'y')then
                  qdf(n,ii)=qdf(n,ii)+fake(ii)*eff_bare_area*tdum
                endif
              endif

cd     if(iopt.eq.3)print*,'In runof5, passed 801'

!             OVERLAND FLOW (DIRECT RUNOFF)
              if(over(n)*type1.le.0)then
!               WE HAVE A FLOODPLAIN AND THERE CAN BE RIVER INFLOW
                if(d1(n,ii).le.ds(ii))then
!                 WATER DEPTH IS LESS THAN DEPRESSION STORAGE
!                 AND THERE IS NO SURFACE FLOW
                  q1(n,ii)=0.0
!                 AND D1 REMAINS UNCHANGED BY SURFACE RUNOFF
                else
!                 THE CONVERSE - THE MAX POSSIBLE FLOW IS: 
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
                  q1(n,ii)=(d1(n,ii)-ds(ii))**1.66*sl2(n)*step2/r3(ii)
!                 UNITS: d1,ds IN MM; step2 IN KM^2 
!                 SO R3 HAS DIMENSIONS ACCORDING AS IN MANNING
!                 THE WATER REMAINING AT END OF TIME STEP 
                  dend=d1(n,ii)-q1(n,ii)/tdum
!                 SURFACE FLOW CAN NOT TAKE WATER BELOW DEPRESSION 
!                 STORAGE
d                 if(iopt.eq.3)print*,'In runof5, passed 901'
                  if(dend.ge.ds(ii))then
!                   THERE IS ENOUGH WATER FOR MAX RUNOFF AND        
                    d1(n,ii)=dend
                  else
!                   ALL WATER ABOVE DEP. STOR. IS DRAINED IN THIS DT.
!                   (AND PROBABLY THE SURFACE ROUGHNESS IS TOO LOW!!
                    q1(n,ii)=(d1(n,ii)-ds(ii))*tdum
                    d1(n,ii)=ds(ii)
                  endif
                endif
!               RF() (IN M) IS RUNOFF IN MM FOR THIS DT NEEDED BY SED
                rf(n,ii)=q1(n,ii)/tdum
!               ADJUST FOR PROPER CONTRIBUTING AREA:
                q1(n,ii)=q1(n,ii)*eff_bare_area
cd               if(iopt.eq.3)print*,'In runof5, passed 1001'

!               INTERFLOW AND DRAINAGE:
!               WHERE UZS IS IN MM AND QINT IS CUBIC METER PER SECOND
!     REV. 8.60 - Nov.  14/97 -   ADDED SL1 TO THE INTERFLOW CALCULATION
!     REV. 8.99b  Sept. 27/00 -   DIVVY UP INTERFLOW & DRAINAGE
!               Added check - AB, April 19, 2002
                if(uzs(n,ii).gt.retn(ii))then
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
!                  duz=x4(ii)*((uzs(n,ii)-retn(ii))**A10)*sl1(n)
                  duz=rec(ii)*((uzs(n,ii)-retn(ii))**A10)*sl1(n)
                else
                  duz=0.00000
!                 duz=1.0e-10
                endif
cd     if(iopt.eq.3)print*,'In runof5, passed 1001-2'

!               DRAINAGE OF UZS TO GROUND WATER:
!               THRESHOLD FOR DRAINAGE 
!               IF VOL. WATER CONTENT LESS THAN RETN THEN NO
!               DRAINAGE OCCURS
!               AK2 MUST BE LT 1 - CHECKED IN PARAM.FOR
                if(uzs(n,ii).gt.retn(ii))then
!                 Added check - AB, April 19, 2002
                  drng(n,ii)=ak2(ii)*(uzs(n,ii)-retn(ii))    !**power1
                else
                  drng(n,ii)=0.0000
!                 drng(ii)=1.0e-10
                endif
                drng(n,ii)=amax1(0.000,drng(n,ii))

!               DIVVY UP THE OUTFLOW
                demand=duz+drng(n,ii)
                if(demand.gt.0.0)then
                  if(uzs(n,ii)-retn(ii).ge.demand)then
!                   THERE IS ENOUGH WATER TO SATISFY BOTH
                    uzs(n,ii)=uzs(n,ii)-duz-drng(n,ii)
                  else
!                   DIVVY UP THE WATER
                    supply=uzs(n,ii)-retn(ii)
                    fraction=supply/demand
                    duz=duz*fraction
!     rev. 10.1.84 May   09/17  - NK: Put drng(n,ii)=drng(n,ii)*fraction back into runof6
                    drng(n,ii)=drng(n,ii)*fraction
                    drng(n,ii)=amax1(0.000,drng(n,ii))  ! TS: ADDED ERROR CHECK
                    uzs(n,ii)=retn(ii)
                  endif
                  qint(n,ii)=duz*eff_bare_area*tdum
                  if(n.eq.nnprint)then
                    sq1(ii)=sq1(ii)+rf(n,ii)*eff_bare_area*tdum
                    sqint(ii)=sqint(ii)+duz*eff_bare_area*tdum
                  endif
                  if(n.eq.nnprint)then
                    sdrng(ii)=sdrng(ii)+drng(n,ii)*eff_bare_area  
!                   TS: added *tdum b/c lower zone inflow: Nov 1/07
                    slzinflw=slzinflw  ! FLOW
     *                     +drng(n,ii)*eff_bare_area*tdum 
                  endif
!                 DRAINAGE IS WEIGHTED HERE TO ACOCUNT FOR 
!                 LAND COVER AREA OF EACH CLASS        ! lzs is depth
                  lzs(n)=lzs(n)+drng(n,ii)*aclass(n,ii)*(1-sca(n,ii))
                endif
cd     if(iopt.eq.3)print*,'In runof5, passed 1301'

!               rechrg() is the incremental recharge for output to MODFLOW
!               This could be for instance in 24 hour increments (hard coded now) 
!               The recharge is adjusted for a grid with frac=1.0 !!!
                rechrg(n)=rechrg(n)                  ! adjusted depth
     *                   +drng(n,ii)*eff_bare_area

!     rev. 9.1.48  Dec.  08/03  - NK: sumrechrge() added to get total recharge
!               sumrechrg() is the cummulative recharge for the whole run
!               and is for information wrt doing a water balance AND
!               to see what the recharge is on each grid.
                sumrechrg(n)=sumrechrg(n)+drng(n,ii)*eff_bare_area
                qdrng(n)=qdrng(n)+drng(n,ii)*eff_bare_area*tdum
                if(frcflg.eq.'y') qdrng2(n,ii)=
     *                       qdrng2(n,ii)+drng(n,ii)*eff_bare_area*tdum
!               STREAMFLOW CONTRIBUTION: OVERLAND & INTERFLOW
!               WEIGHTED ACCORDING TO LAND COVER AMOUNT
                qr(n)=qr(n)+qint(n,ii)+q1(n,ii)
              else
!               NO FLOOD PLAIN : TYPE1 = 1.0 IN MAIN PROGRAM
!               WHEN OVER*TYPE1 .GE. 0, CHANNEL INFLOW CALCS ARE
!               BYPASSED 
!               CHANNELS ARE FULL AND THERE IS NO DIRECT RUNOFF OR
!               INTERFLOW
                qr(n)=0.0
                qbase(n)=0.0
              endif
cd     if(iopt.eq.3)print*,'In runof5, passed 1401'

!             EVAPORATION:
!             WATER IS ONLY LOST TO EVAPORATION IF WATER IS AVAILABLE
!             SINCE THE HOURLY AMOUNTS ARE SMALL, WE'RE NOT BOTHERING 
!             TO BRING THE UZS TO 0.0
!       REV 7.9  PROVISION TO ALLOW THE USE OF ORIGINAL AET TABLE DATA
              if(flgevp2.lt.1.0)then
!               alternative method - aet NOT called
                evaptemp=amin1(uzs(n,ii),evap(ii,mon))
                if(evaptemp.gt.0.0)then
!                 WATER AVAILABLE AND PET >0.0, TAKE OFF EVAPORATION:
                  uzs(n,ii)=uzs(n,ii)-evaptemp
                  evt(n,ii)=evt(n,ii)+evaptemp*(1-sca(n,ii))
                  eloss(n)=eloss(n)+evaptemp*aclass(n,ii)*(1-sca(n,ii))
                  if(frcflg.eq.'y')then
                    evloss(n,ii)=evaptemp
                  endif
                endif
!               RETN MAY HAVE TO BE CONSIDERED <<<<<<      ???
                pet(n,ii)=0.0
              endif
!             END OF REV 7.9


c            else
c!             SURFACE WATER:
c!             PERMEABILITY IS -VE FOR WATER
c!             RAIN FALLS ON WATER SURFACE AND IS DIRECTLY ADDED TO 
c!             RIVER FLOW. SET AK(ii) TO 0.0 OR LESS FOR THIS OPTION
c!             REV 7.9 NEW CALCULATION BASED ON NEW EVAPORATION
c!             REV. 8.85 - Oct. 12/98 - FIXED RAIN & SNOW ON WATER CLASS
c              d1(n,ii)=0.0
c              df(n,ii)=0.0
c              uzs(n,ii)=0.0
cd             if(iopt.eq.3)print*,'In runof5, passed 1501'
c!             qstream(n) is added to channel inflow in route
c!             qstream(n) = net precip in mm converted to cms
c              qstream(n)=r(n,ii)*eff_bare_area*tdum        
c!              qstream(n)=amax1(0.1e-10,qstream(n))
c              iiwater=ii
c              if(.not.rd_evp_flg)then
c!               if the water evaporation is read in, skip this part
c!               strloss is computed in sub
c!               EVAPORATION IN THE WATER CLASS IS EQUAL TO
c!               POTENTIAL EVAPORATION AND HAS TO BE TAKEN FROM
c!               RIVER STORAGE AS A FLOW
c!     rev. 9.1.66  Oct.  17/04  - NK; pet*fpet for loss from water instead of pet
c!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!change  nov 25/09 nk
c                if(r(n,ii).gt.0.0)then
c                  ev(n,ii)=0.0    ! ev calculated in AET
c                endif
c                
c!     REV. 10.1.16 Jan.  11/16  - NK: Added subroutine ice_factor.f
c                if(ice_fctr(n).lt.0.0)then
c                  strloss(n)=0.0
c                  ev(n,ii)=0.0
c                else
c!     rev. 9.5.61  Sep.  03/09  - NK: bug/eloss - added water class for wfo weighted et
c!                   for the wfo file:
c                  strloss(n)=ev(n,ii)*eff_bare_area*tdum
c!                 strloss is taken from channel storage in route
c                  if(ireach(n).eq.0)then
c                    if(store2(n).gt.strloss(n)*3600.0)then
c                      evt(n,ii)=evt(n,ii)+ev(n,ii)         ! for rffxx.txt
c!                     for ensim: prorated for sca 
c!                     100% snow cover - no et
c!                       0% snow cover - no reduction
c                      sum_et(n,ii)=sum_et(n,ii)+ev(n,ii)*(1.0-sca(n,ii))   
c                      eloss(n)=eloss(n)+
c     *                              ev(n,ii)*aclass(n,ii)*(1-sca(n,ii))
c
c!     REV. 10.1.17 Jan.  11/16  - NK: Added fpetLakeOverride factor
c!                     Temporary fix for dealing with rogue lakes - i.e. where fpet
c!                     just doesn't fit all the lakes
c!                     Maybve the lake evaporation model will fix this - maybe not.
c!                     It has to be done here so the tracer & isotope stuff 
c!                     is not affected
c                    else
c                      strloss(n)=0.0
c                    endif     !  store2(n).gt.............
c                  endif     !  ireach(n).eq.0
c                endif     !  ice_fctr(n).lt.0.0
c              endif     !  .not.rd_evp_flg
            endif   ! :end SURFACE WATER:  

!           SEDIMENT COMPONENT - REMOVE IF NOT USED
!           EROSION ON BARE GROUND ONLY
            if(sedflg.eq.'y'.or.sedflg.eq.'Y')then
              qs(n,ii)=q1(n,ii)
              hsed(n,ii)=d1(n,ii)-ds(ii)
            endif
d     if(iopt.eq.3)print*,'In runof5, passed 1601'

!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  ++++++++++++++   RUNOFF FOR SNOWCOVERED GROUND   +++++++++++++
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!           BYPASS SNOW COVERED GROUND WHEN THERE IS NO SNOW !!!!
            if(snwflg.eq.'y')then
            if(sca(n,ii).gt.0.001)then
            if(akfs(ii).gt.0.0)then
!             NON-water classes
              d1fs(n,ii)=d1fs(n,ii)+fexcess(n,ii)
                fakefs(ii)=akfs(ii)*thr*(1.+effpor(n,ii)
     *                  *(d1fs(n,ii)+potfs(ii))/(uzsfs(n,ii)+0.00001))
cd               if(iopt.eq.3)print*,'In runof5, passed 1701'
!             INFILTRATION: 
              if(d1fs(n,ii).le.fakefs(ii))then
!               SURFACE STORAGE LT INFILTRATION CAPACITY 
!               AND ALL THE WATER IS INFILTRATED-NO RUNOFF
                uzsfs(n,ii)=uzsfs(n,ii)+d1fs(n,ii)
!     rev. 10.1.55 Nov.  30/16  - NK: Fixed sumf & sumffs in runof6
c                sumffs(n,ii)=sumffs(n,ii)+d1fs(n,ii)
                sumffs(n,ii)=sumffs(n,ii)+d1fs(n,ii)*eff_sc_area
                if(frcflg.eq.'y')then
                  qdffs(n,ii)=qdffs(n,ii)+d1fs(n,ii)*eff_sc_area*tdum  ! fixed nk
                endif
!                prorated for sca Jul. 7/06 nk
                d1fs(n,ii)=0.0
              else
!               SURFACE STORAGE > INFILTRATION CAPACITY
                uzsfs(n,ii)=uzsfs(n,ii)+fakefs(ii)
!     rev. 10.1.55 Nov.  30/16  - NK: Fixed sumf & sumffs in runof6
c                sumffs(n,ii)=sumffs(n,ii)+fakefs(ii)
                sumffs(n,ii)=sumffs(n,ii)+fakefs(ii)*eff_sc_area
!                prorated for sca Jul. 7/06 nk
                d1fs(n,ii)=d1fs(n,ii)-fakefs(ii)
                if(frcflg.eq.'y')then
                  qdffs(n,ii)=qdffs(n,ii)+fakefs(ii)*eff_sc_area*tdum
                endif
              endif

cd     if(iopt.eq.3)print*,'In runof5, passed 1801'

!             OVERLAND FLOW (DIRECT RUNOFF):
              if(over(n)*type1.le.0)then
!               WE HAVE A FLOOD PLAIN AND THERE CAN BE INFLOW
!               TO THE RIVER
                if(d1fs(n,ii).le.dsfs(ii))then
!                 WATER DEPTH IS LESS THAN DEPRESSION STORAGE
!                 AND THERE IS NO SURFACE FLOW
                  q1fs(n,ii)=0.0
!                 AND D1 REMAINS UNCHANGED BY SURFACE RUNOFF
                else
!                 THE CONVERSE - THE MAX POSSIBLE FLOW IS
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
                  q1fs(n,ii)=(d1fs(n,ii)-dsfs(ii))**1.66
     *                         *sl2(n)*step2/r3fs(ii)

!                 AND THE WATER REMAINING AT END OF TIME STEP
                  dend=d1fs(n,ii)-q1fs(n,ii)/tdum
!                 SURFACE FLOW CAN NOT TAKE WATER BELOW DEPRESSION
!                 STORAGE
                  if(dend.ge.dsfs(ii)) then
!                   THERE IS ENOUGH WATER FOR MAX RUNOFF AND
                    d1fs(n,ii)=dend
                  else
!                   ALL WARER ABOVE DEP. STOR. IS DRAINED IN THIS DT.
!                   (AND PROBABLE THE SURFACE ROUGHNESS IS TOO LOW !!)
                    q1fs(n,ii)=(d1fs(n,ii)-dsfs(ii))*tdum
                    d1fs(n,ii)=dsfs(ii)
                  endif
                endif
cd              if(iopt.eq.3)print*,'In runof5, passed 1901'
!               RFFS() IS RUNOFF IN MM FOR THIS DT NEEDED BY SED
                rffs(n,ii)=q1fs(n,ii)/tdum
!               adjust for contributing area
                q1fs(n,ii)=q1fs(n,ii)*eff_sc_area

!               INTERFLOW:
!               WHERE UZS IS IN MM AND QINT IS CUBIC METER PER SECOND
!   REV. 8.60 - Nov.  14/97 -     ADDED SL1 TO THE INTERFLOW CALC
!   REV. 8.99b  Sept. 27/00 -     DIVVY UP INTERFLOW & DRAINAGE
!               Added check - AB, April 19, 2002  added here Apr. 08/03 nk
                if(uzsfs(n,ii).gt.retn(ii))then
!                if(uzs(n,ii)-retn(ii).gt.1.0e-10)then   
                              ! changed Jun. 3/02 nk
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
!                  duz=x4(ii)*((uzsfs(n,ii)-retn(ii))**a10)*sl1(n)
                  duz=rec(ii)*((uzsfs(n,ii)-retn(ii))**a10)*sl1(n)
                else
                  duz=0.00000
!                 duz=1.0e-10
                endif

!               DRAINAGE OF UZS TO GROUND WATER:
!               THRESHOLD FOR DRAINAGE - IF VOL. WATER CONTENT LESS 
!               THAN RETN THEN NO DRAINAGE OCCURS 
!               ak2fs must be lt 1 - checked in param.for
                if(uzsfs(n,ii).gt.retn(ii))then
!                if(uzs(n,ii)-retn(ii).gt.1.0e-10)then   ! changed Jun. 3/02 nk
!               Added check - AB, April 19, 2002
                  drngfs(n,ii)=ak2fs(ii)*(uzsfs(n,ii)-retn(ii))    !**power1
                else
                  drngfs(n,ii)=0.0000
!                 drngfs(ii)=1.0e-10
                endif
                drngfs(n,ii)=amax1(0.000,drngfs(n,ii))  ! TS: ADDED ERROR CHECK

!               DIVVY UP THE UZ OUTFLOW 
cd               if(iopt.eq.3)print*,'In runof5, passed 2001'
                demand=duz+drngfs(n,ii)
                if(demand.gt.0.0)then
                  if(uzsfs(n,ii)-retn(ii).ge.demand)then
!                   THERE IS ENOUGH WATER TO SATISFY BOTH
                    uzsfs(n,ii)=uzsfs(n,ii)-duz-drngfs(n,ii)
                  else
!                   DIVVY UP THE WATER
                    supply=uzsfs(n,ii)-retn(ii)
                    fraction=supply/demand
                    duz=duz*fraction
!     rev. 10.1.84 May   09/17  - NK: Put drng(n,ii)=drng(n,ii)*fraction back into runof6
                    drngfs(n,ii)=drngfs(n,ii)*fraction
                    drngfs(n,ii)=amax1(0.000,drngfs(n,ii))  ! TS: ERROR CHECK
                    uzsfs(n,ii)=retn(ii)
                  endif
                  qintfs(n,ii)=duz*eff_sc_area*tdum
                  if(n.eq.nnprint)then
                    sq1fs(ii)=sq1fs(ii)+rffs(n,ii)*eff_sc_area*tdum
                    sqintfs(ii)=sqintfs(ii)+duz*eff_sc_area*tdum
                  endif
                  if(n.eq.nnprint)then
                    sdrngfs(ii)=sdrngfs(ii)+drngfs(n,ii)*eff_sc_area
                    slzinflw=slzinflw  !FLOW
     *                      +drngfs(n,ii)*eff_sc_area*tdum
                  endif
!               DRAINAGE IS WEIGHTED HERE TO ACCOUNT FOR 
!               LAND COVER AREA OF EACH CLASS
                lzs(n)=lzs(n)+drngfs(n,ii)*aclass(n,ii)*sca(n,ii)  ! depth
              endif
cd             if(iopt.eq.3)print*,'In runof5, passed 2101'
!             The recharge is adjusted for a grid with frac=1.0 !!!
              rechrg(n)=rechrg(n)               ! adjusted depth
     *                 +drngfs(n,ii)*eff_sc_area
!     rev. 9.1.48  Dec.  08/03  - NK: sumrechrge() added to get total recharge
!             see note above
              sumrechrg(n)=sumrechrg(n)+drngfs(n,ii)*eff_sc_area
!                               found mistake here 30/06/00
!             TS: MADE THIS DRNG FROM SNOWMELT FOR TRACER S/R (08/27/03)
              qdrngfs(n)=qdrngfs(n)+drngfs(n,ii)*eff_sc_area*tdum
              if(frcflg.eq.'y') qdrngfs2(n,ii)=
     *                     qdrngfs2(n,ii)+drngfs(n,ii)*eff_sc_area*tdum
!             STREAMFLOW CONTRIBUTION: OVERLAND & INTERFLOW
!             WEIGHTED ACCORDING TO LAND COVER AMOUNT
              qr(n)=qr(n)+qintfs(n,ii)+q1fs(n,ii)
            else
!             NO FLOOD PLAIN: TYPE1 = 1.0  IN MAIN PROGRAM
!             WHEN OVER*TYPE1.GE.0, CHANNEL INFLOW CALCS ARE BYPASSED
!             CHANNELS ARE FULL AND THERE IS NO DIRECT RUNOFF OR
!             INTERFLOW
              qr(n)=0.0
              qbase(n)=0.0
            endif
cd           if(iopt.eq.3)print*,'In runof5, passed 2201'
   

c            else                    !IF(akfs(ii).gt.0.0)THEN
c!             SURFACE WATER:
c!             PERMEABILITY IS -VE FOR WATER
c!             FEXCESS FALLS ON WATER SURFACE AND IS DIRECTLY ADDED TO 
c!             RIVER FLOW. SET AK(ii) TO 0.0 OR LESS FOR THIS OPTION
c!             REV 7.9 UPDATED EVAPORATION    
c!                BUT THERE IS NO EVAP FOR MELTED SNOW
c!             REV. 8.85 - Oct. 12/98 - FIXED RAIN & SNOW ON WATER CLASS
c              d1fs(n,ii)=0.0
c              dffs(n,ii)=0.0
c              uzs(n,ii)=0.0
c              iiwater=ii
c!             qstream(n) is added to channel inflow in route
c!             qstream(n) = net melt+rain in mm converted to cms
c              qstream(n)=qstream(n)+fexcess(n,ii)*eff_sc_area*tdum
c              
c              strloss(n)=0.0
c              
c!             NO EVAPORATION WHEN WATER IS COVERED BY SNOW!
c!                 this is in the last version of SPL8
c               
c             
            endif           !IF(akfs(ii).gt.0.0)THEN
            endif           !IF(sca(n,ii).gt.0.001)THEN
            endif           !IF(snwflg.eq.'y')THEN
!           END OF SNOW COVERED AREA COMPUTATIONS

          endif                         ! IF(aclass().gt.0.0)THEN
cd     if(iopt.eq.3)print*,'In runof5, passed 2301'
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++   fs end 

          endif                     ! if(ii.NE.classcount-2)

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!       TS:ADDED SUMS HERE FOR TRACER S/R (08/27/03)
!       INCLUDES impervious, wetland classes and water class
!          LOOPING FROM ii=1,classcount
          sumq1(n)=sumq1(n)+q1(n,ii)
          sumq1fs(n)=sumq1fs(n)+q1fs(n,ii)
          sumqint(n)=sumqint(n)+qint(n,ii)
          sumqintfs(n)=sumqintfs(n)+qintfs(n,ii)
          if(frcflg.eq.'y')then 
            storeSW2(n,ii)=storeSW2(n,ii)+(d1(n,ii)*eff_bare_area 
     *                 +d1fs(n,ii)*eff_sc_area)*xdum
            storeIF2(n,ii)=storeIF2(n,ii)+(uzs(n,ii)*eff_bare_area
     *                     +uzsfs(n,ii)*eff_sc_area)*xdum   
          endif
          

        end do                ! classcount loop
        
!      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!        if(frcflg.eq.'y')then   ! TS: moved to below "call baseflow" where lzs is defined
!          storeGW2(n)=lzs(n)*xdum*frac(n)
!        endif
cd       if(iopt.eq.3)print*,' checkpoint 5 in runoff'
!     rev. 9.4.08  May.  29/07  - NK: changed baseflow argument list
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call baseflow(n,dlz,sdlz,tdum)
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cd     if(iopt.eq.3)print*,' checkpoint 6 in runoff'
        if(frcflg.eq.'y')then
          storeGW2(n)=lzs(n)*xdum*frac(n)
        endif

! TS - ADDED ROUTING PROVISION FOR WETLANDS - NOV 2000
! * * * * * * * * * * * * * * * * * * * * * * * * * *
!      rev. 9.8.61  May   22/13  - NK: Introduced flag1 to speed up runof6
!      rev. 9.8.62  May   22/13  - NK: Fixed bug in runof6: (classcount-3) to (classcount-2) 
c        if(aclass(n,classcount-2).gt.0.0.and.flag1(n))then
!      rev. 9.8.77  Jul   08/13  - NK: Made universal the use of wetland_flag(n)
!      rev. 9.8.79  Jul   19/13  - NK: Fixed wetland conditional screwed up with rev 9.8.77 in runof6
        if(wetland_flag(n))then
!         for wetlands, qlz(n) is added to qiwet2(n) in route
!      rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!                we have to assume that if there is a glacier, 
!                there is no wetland in that grid
!                i.e. we run this code only when NOT in a glacier grid
!                only glacier flow bypasses the wetland
          qr(n)=amax1(qr(n),0.1e-10)
        else
!         not in a wetland so qlz        
!         water is added directly to the stream in route
          qr(n)=qr(n)+qlz(n)
          qr(n)=amax1(qr(n),0.1e-10)
        endif
! * * * * * * * * * * * * * * * * * * * * * * * * * 
!       HAVE TO DIVIDE BY FRAC TO NORMALIZE TO FULL GRID 
!       THE WATER AREA IS INCLUDED in the water class above
!     rev. 9.2.42  Jun.  20/06  - NK: water class included in the water balance
        sumrff(n)=sumrff(n)+qr(n)/tdum/frac(n)
!       ADD RAIN & SNOWMELT ON WATER AREA TO RUNOFF TOTAL
!       THIS IS NOT INCLUDED IN THE WATERBALANCE CHECK
!       WRITE ALL THE STUFF FOR THE RFF*.TXT FILES
!       FOR DEBUGGING AND PROCESS PLOTS
cd        if(iopt.eq.3)print*,' checkpoint 7 in runoff'
        if(iopt99.or.debug_output)then
c        if()then
          if (n.eq.nnprint)then
            mzz=mz-24*(mz/24)
!           if(mzz.eq.1)then
             classcount_9=min0(9,classcount)
             amon=float(mon)
             aju=float(ju)
!            this section repeated below
             if(snwflg.eq.'y'.or.vapflg.eq.'y')then
!              this is needed because memory is allocated 
!              only for these cases
               tempv1=tempv(n)
               tempvmin1=tempvmin(n)
               rh1=rh(n)
             else
               tempv1=-999.
               tempvmin1=-999.
               rh1=-999.
             endif

!     rev. 9.7.13  Nov.  22/10  - NK: Changed the outfiles.txt for more 50 rff classes
!     rev. 9.9.69  Jun.  10/15  - NK: prevent write ro rff if there is no class area
             do ii=1,classcount
               if(aclass(n,ii).gt.0.0)then
                 write(900+ii,6000)time,intevt(n,ii),evt(n,ii),
     *           p(ipr,jpr),sump(n),ssumr(n,ii),
     *           amin1(fake(ii),1000.),amin1(fakefs(ii),1000.),
!     rev. 10.2.39 Nov.  15/18  - NK: changed snowc(n,ii) to snowc(n,ii)*sca(n,ii) in runof6
     *           sca(n,ii),snowc(n,ii)*sca(n,ii),d1(n,ii),d1fs(n,ii),
     *           sumf(n,ii),sumffs(n,ii),uzs(n,ii),uzsfs(n,ii),lzs(n),
     *           q1(n,ii),q1fs(n,ii),qint(n,ii),qintfs(n,ii),qlz(n),
     *           drng(n,ii),drngfs(n,ii),qr(n),qstream(n),strloss(n),
     *           sumrff(n),fexcess(n,ii),glmelt(n),fmadj(n),
     *           sq1(ii),sq1fs(ii),sqint(ii),sqintfs(ii),
     *           sdrng(ii),sdrngfs(ii),slzinflw,sdlz,amon,aju,def(n,ii),
     *           tempv1,tempvmin1,rh1,           !psmear(n),punused(n),
     *           api(n,ii)*100.0,sublim(n,ii),sum_sublim(n,ii),
     *           v(n,ii),wcl(n,ii),sum_pet(n,ii),sum_et(n,ii),pet(n,ii),
     *           tto(n),ttomin(n)
               endif
             end do
!           endif
          endif
        endif
        
        
        
        
        qstrm(n)=qstream(n)
!     TRISH:  this means that qstream never got added to the lakes ?????????????
cccccccccccccccccccccccccc        qstream(n)=0.0                   
cd        if(iopt.eq.3)print*,' checkpoint 8 in runoff'

!       MODIFY THE UNSATURATED ZONE SOIL MOISTURE
!       THIS USED TO BE AT THE START OF THE LOOP - NO GOOD
!       ADJUST SOIL MOISTURE AFTER ALL CALCS 
!       THIS SHOULD BE DONE SEPERATELY FOR SCA - USE FEXCESS
!       INSTEAD OF RAIN
!       REV 7.9 EVAPORATION    
!       MODIFY THE UNSATURATES ZONE SOIL MOISTURE

!       SOIL MOISTURE IS ADJUSTED FOR DRYING OR NEW PRECIP
!       IF IT IS BELOW FREEZING, NO ADJUSTMENTS ARE MADE
!       THIS IS A KLUDGE
!       AND may HAVE TO BE FIXED SOMEDAY  <<<<<<<
        do ii=1,classcount
          if(snwflg.eq.'y'.or.vapflg.eq.'y')then
            if(ttemp(i,j).gt.0.0)then
              api(n,ii)=api(n,ii)*a51+p(i,j)/100.0          
!              effpor(n,ii)=por-api(n,ii)
!     rev. 9.4.06  May.  09/07  - NK: replaced por with spore(n,ii) in runof6
              effpor(n,ii)=spore(ii)-api(n,ii)
              effpor(n,ii)=amax1(0.0001,effpor(n,ii)) 
              effpor(n,ii)=amin1(por   ,effpor(n,ii))
            endif
          else
!           this bit of code has to be repeated because ttemp() 
!           has not been allocated memory for this case so
!           can't be used in a conditional
            api(n,ii)=api(n,ii)*a51+p(i,j)/100.0          
!            effpor(n,ii)=por-api(n,ii)
!     rev. 9.4.06  May.  09/07  - NK: replaced por with spore(n,ii) in runof6
            effpor(n,ii)=spore(ii)-api(n,ii)
            effpor(n,ii)=amax1(0.0001,effpor(n,ii)) 
            effpor(n,ii)=amin1(por   ,effpor(n,ii))
          endif
        end do     !  ii=
      endif            ! IF(slope(n).gt.0.0)THEN
cd      if(iopt.eq.3)print*,n,ii,time,' checkpoint 9 in runoff'
        end do                !  n=1,naa loop
        
        
cc!$OMP END DO
      

!     * * * * * * *   END OF MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!     * * * * * * *   END OF MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!     * * * * * * *   END OF MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!     * * * * * * *   END OF MAIN LOOP   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

!     QR(n)=TOTAL STREAM INFLOW FROM grid(N) AND IS ROUTED
!     INTO THE DOWNSTREAM ELEMENT.
d      if(iopt.eq.3)print*,' checkpoint 10 in runoff'
      RETURN

1003  CONTINUE
!     END OF RUN
!     WRITE THE WATER BALANCE INFORMATION
d      if(iopt.eq.3)print*,' checkpoint 11 in runoff'

      if(numa.eq.0.and.resumflg.eq.'n'.and.iopt99)call watbal(3)
d      if(iopt.ge.3)write(*,5101)time
d      if(iopt.eq.3)print*,' checkpoint 12 in runoff'

      RETURN

! FORMATS
 5101 format(' at RETURN 6007 time=',f10.2)
 5102 format(/' por = ',f10.3/)
 5190 format(/,' qlzfrac =',f10.2,'  in runof5 <<<<<<<'/)

 6007 format(' n,ipr,jpr,aclass,frac/',3i5,17f8.2)
c 6000 format(' ',16(g8.2,','),1(g8.0,','),50(g8.3,','))
c 6000 format(' ',f8.0,999(',',g10.3))
      
 6000 format(' ',f8.0,3(',',g10.3),f10.1,999(',',g10.3))
     
 6001 format('     time',7x,'v',2x,'intcap',7x,'p',4x,'sump',7x,'r',
     *6x,'d2',6x,'q2',6x,'qr')
 6002 format(    '       time,    intevt,      evt,        p,',
     *  '      sump,      ssumr,     fake,    fakefs,',
     *  '       sca,     snowc,        d1,      d1fs,',
     *  '         sumf,    sumffs,',
     *  '     uzs,     uzsfs,       lzs,        q1,      q1fs,',
     *  '      qint,    qintfs,       qlz,      drng,    drngfs,',
     *  '       qr,    qstream,   strloss,    sumrff,',
     *  '   fexcess,    glmelt,  fmadjust,       sq1,     sq1fs,',
     *  '      sqint,  sqintfs,     sdrng,   sdrngfs,   slzinfl,',
     *  '      sdlz,     month,   jul_day,  heat_def,     tempv,',
     *  '     tempvmin,     rh,    API,',
     *  '     sublm,   sumsublm,       v,        wcl,',
     *  '     sumpet,     sumet,     pet,        tto,     ttomin')
 6003 format(9i8,f8.3)
 6004 format(5f12.5)
 6005 format(' n=',i5,' ii=',i5,' uzs=',f10.5) 
 6006 format(' ipr,jpr/=',2i5)
 6027 format(' ii=',i2,' smc=',f5.2,' effpor=',f5.2,
     * ' -> lower limit is set <-')
 6040 format(8a)
 6100 format
     *('0','    n    i    j        da       lzs     qbase in runof5'/)
 6201 format(i10,7f10.1)
 7010 format(' sump       ',f8.1/)
 7011 format(' ssumr       ',16f8.1/)
 7012 format(' d1         ',16f8.1/)
 7013 format(' sq1        ',16f8.1/)
 7014 format(' sumf       ',16f8.1/)
 7015 format(' usz        ',16f8.1/)
 7016 format(' sqint      ',16f8.1/)
 7017 format(' sdrng      ',16f8.1/)
 7021 format(' sexcess    ',16f8.1/)
 7022 format(' d1fs       ',16f8.1/)
 7023 format(' sq1fs      ',16f8.1/)
 7024 format(' sumffs     ',16f8.1/)
 7025 format(' uzsfs      ',16f8.1/)
 7026 format(' sqintfs    ',16f8.1/)
 7027 format(' sdrngfs    ',16f8.1/)
 7030 format(' aclass     ',16f8.2/)
 7031 format(' slzinflw   ',16f8.1/)
 7032 format(' lzs        ',16f8.1/)
 7033 format(' sdlz       ',f8.1/)
 7034 format(' sumrff     ',16f8.1/)
 7050 format('     ')
! 7401 format('    time fpet2(n)   uzsi'
!     * ,' tto(n)   intev ev(n.ii) pet(n.ii)'
!     * ,'uzs(n.ii)')
 7402 format('    time',<classcount>('   ET   ','   E    '))

! 8429 format(//' in runof5,flz/',5f10.6/
!     * ' wrong parameter value - change flz(',i2,') to +ve value'/)
59999 format(i5,4f12.5,'runof6 //')


      RETURN

      END SUBROUTINE runof6
