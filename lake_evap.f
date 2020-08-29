       SUBROUTINE lake_evap
!***********************************************************************
!    Copyright (C) 1987 by Tegan Holmes and Nicholas Kouwen 
        
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
!  Written by: Tegan Holmes 
!  With initial work by: Phil Slota, March 2013
!       internal revisions:
!          1.0 running
!          2.0 unit conversions, heat flux, net rad v2
!          2.1 solar rad, evap shut down in rain, water T alpha, ice cover
!          2.2 lake depth read-in created
!          3.0 Cleaned out fossil code
!     	   4.0 Removed wind and Granger
!
!  alamb - lambda - latent heat of vapourization
!  delt - delta (slope of saturation vap press curve)
!  gam - psychrometric constant
!  strloss(n) - Evaporation for lakes in cms
!  tempv(n) vectorized air temperature (deg C) for grid n
!  tempw(n) vectorized water temperature (deg C) for grid n
!  alpha - prieslty-taylor alpha coefficient
!  heatflux - heat storage flux in grid
!  net_solar_radv - net solar radiation flux in grid
!  ket_hr - extra-terrestrial solar radiation at middle of hour
!  ice_cover - positive is ice-on, based on water temperature
!
! NOTE: petn is adjusted based on the time increment
! i.e. for 1 hour it is multiplied by 3.6 since there are 3600 seconds
! per hour and is divided by 1000 to conserve units
!**********************************************************************

!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model to WATFLOOD
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
!     rev. 10.1.49 Nov.  08/16  - NK: Overhauled lake evaporation
!     rev. 10.2.36 Oct.  17/18  - NK: Initialized evap_rate(n) in lake_evap.f 

      USE area_watflood
      implicit none
      save

      real*8      :: gam,delt,ddg,declin,eccent,ths,
     *               heatflux,short_wave,long_wave,
     *               net_solar_radv,lkalpha,tau,Ket,Kso,evapE
      integer*4   :: n,ii,l,hour_count,ios,rbin
      real*4      :: dlta,sindlta,cosdlta,tandlta,dr,ws
      logical     :: exists,firstpass,warningflag
      REAL*4, DIMENSION(:), ALLOCATABLE :: ftempold,tempw,ilat
     *        ,radlat,Temp_max,Temp_min,Temp_tot,Temp_avg,evap_rate
     *        ,heatstore,tau_so,tdum,logD
      REAL*4, DIMENSION(12,10) :: albedo(12,10)= RESHAPE((/ 100.,100.,
     * 30.1,29.3,17.1,14.8,16.0,24.6,34.2,100.,100.,100.,100.,30.1,31.9,
     * 22.5,16.0,13.1,14.5,20.6,29.4,30.5,100.,100.,30.1,33.8,22.9,14.8,
     * 11.6,11.2,11.4,13.4,20.2,31.3,30.1,100.,33.9,24.0,15.5,10.5,8.8,
     * 8.4,8.6,9.8,13.6,21.6,32.1,35.5,22.0,16.1,10.8,8.4,7.5,7.3,7.4,
     * 8.0,9.9,14.4,21.0,24.1,14.5,11.1,8.5,7.3,6.8,6.7,6.8,7.1,8.0,
     * 10.3,13.8,16.1,10.3,8.6,7.3,6.7,6.5,6.4,6.4,6.6,7.1,8.2,10.0,
     * 11.1,8.3,7.4,6.7,6.4,6.3,6.3,6.3,6.4,6.6,7.2,8.1,8.7,7.2,6.7,6.4,
     * 6.3,6.4,6.4,6.4,6.3,6.3,6.6,7.1,7.4,6.6,6.4,6.3,6.4,6.6,6.8,6.7,
     * 6.4,6.3,6.4,6.6,6.8/),(/12,10/))
      logical, DIMENSION(:), ALLOCATABLE ::  choice

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE
      
!     SET FIXED PARAMETERS, INITIAL TEMPERATURE IN FIRST PASS
      data firstpass/.true./
      data warningflag/.true./
!      data eff_bare_area/1.0/

      if(firstpass)then
        allocate(ftempold(na),tempw(na),radlat(na),ilat(na)
     *          ,Temp_tot(na),Temp_avg(na),Temp_max(na),Temp_min(na),
     *          evap_rate(na),heatstore(na),tau_so(na),choice(na),
     *          tdum(na),logD(noresv),stat=iAll)
	  if(iAll.ne.0) STOP 'Error allocating arrays in evap_lake.for'  

!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
!       NOTE:                   NK   Feb. 25/16      
!       for the Mackenzie River Basin Hydraulic Model LKdepth is used as a flag
!       for evaporation from large lakes. 
!       for small lakes and river reaches the depth is set to 1 m
!       so to use the lake evaporation model the lake depth has to be > 1 m
  
           alpha=1.28 !PTC alphaPT for lakes (could vary)
           hour_count=0
           do n=1,naa
              ! calculate latitude of grid in radians ->fine but this only works for basins in LATLONG (which works for Hydro)
      if(lakeflg.eq.'y')then
               ilat(n)=(iymin+yyy(n)*grdn)/60.
              radlat(n)=ilat(n)/180.*3.14159
             ! patch so lake evap will run (not very well..) when not in latlong
              if(radlat(n)<0)radlat(n)=0.
              if(radlat(n)>1.16)radlat(n)=1.16
              
              tau_so(n)=0.75+0.00002*elev(n)

              ftempold(n)=tempv(n)
              Temp_tot(n)=0
              Temp_max(n)=-999
              Temp_min(n)=999
              heatstore(n)=0
      endif         
      
              rbin=ireach(n)
             choice(n)=.false.
              if(rbin.gt.0)then
                if(LKdepth(rbin).gt.1.0.and.lakeflg.eq.'y')then
                  choice(n)=.true.
                  logD(rbin)=log(LKdepth(rbin))
                  if(logD(rbin)<1.)logD(rbin)=1.
                elseif(lakeflg.eq.'y')then
!     rev. 10.2.29 Aug.  21/18  - NK: Added warning for lake depths less than 1 m
                    if(warningflag)then
                    print*
                    print*,'WARNING'
                    print*,'lakeflg = y but lake depth < or = to 1 m'
                    print*,'Lake evaporation model will not work for '
                    print*,'lake depths less than 1 m and is turned off'
                    print*,'for lake #' ,rbin
                    print*
                    print*,'To turn on lake evaporation, set the'
                    print*,'lake depth > 1 m in the ill file'
                    print*,'Depth found =',LKdepth(rbin)
                    print*,'In lake # ',rbin
                    print*
                    print*,'This message given just once'
                    print*,'but be sure to fix all lake depths > 1 m'
                    print*,'to apply lake evaporation model to any lake'
                    print*
                    call BEEPQQ(500,500)
!       take out for unix - non standard
                    pause 'Hit enter to continue  - in lake_evap @ 115'
                    warningflag=.false.
                    endif
                endif  
              endif
             
!     rev. 10.2.36 Oct.  17/18  - NK: Initialized evap_rate(n) in lake_evap.f 
              evap_rate(n)=0.0  ! must be initialized because it's not calculated until time=24
              tdum(n)= !grid_area(n)/3600000.0*frac(n)*aclass(n,ii_water)
     *        1000.*step2/3600.*frac(n)*aclass(n,ii_water)

            end do    
      end if !first pass
      
!	ju=jul_day_now           ! added NK

      ! This is just an internal hour counter, REPLACE (0=midnight)   
      hour_count=hour_count+1
      if(hour_count>23)hour_count=0

!     qstream(n) = net precip in mm converted to cms
         
      do n=1,naa
!       add precip to all water areas        
        qstream(n)=r(n,ii_water)*tdum(n)*(1.0-sca(n,ii_water))
       if(snwflg.eq.'y') qstream(n)=qstream(n)
     *        +fexcess(n,ii_water)*tdum(n)*sca(n,ii_water)
        
       if(choice(n))then
       
!         lake_evap is called in sub for ALL grids
!         so river evaporation is calculated using the pet method 
!         when choice is false - i.e. whn not on a lake OR
!               when lake evap is not used - lakeflf=n   
  
!         choice is true when we're in a coded lake where :
!                           ireach(n)>0
!                      .and.LKdepth(rbin)>1.0
!                      .and.lakeflg.eq.'y'
c         if(hour_count.ne.0)then
        if(mod(int(totaltime),24).ne.0)then             
           Temp_tot(n)=Temp_tot(n)+tempv(n)
           if(tempv(n)>Temp_max(n))Temp_max(n)=tempv(n)
           if(tempv(n)<Temp_min(n))Temp_min(n)=tempv(n)
            
         else ! day end, new lake evap calc:
           Temp_avg(n)=Temp_tot(n)/24.
           
           dlta=0.4093*sin(0.0172142*float(jul_day_now)-1.405)
           dr=1.0+0.033*cos(0.0172142*float(jul_day_now))
           sindlta=sin(dlta)
           cosdlta=cos(dlta)
           tandlta=sindlta/cosdlta
           ws=acos(-tan(radlat(n))*tandlta)

           ! Shortwave Radation:
           declin=0.4093*sin(2*3.14159*(jul_day_now-81)/365)
           Ket=1367./3.14159*dr*
     *           (ws*sinlat(n)*sindlta+coslat(n)*cosdlta*sin(ws))
           Kso=tau_so(n)*Ket
           if(dlyflg)then
            tau=0.16*(1+0.000027*elev(n))*dly_diff(n)**0.5
           else
            tau=0.16*(1+0.000027*elev(n))*(Temp_max(n)-Temp_min(n))**0.5
           endif
        
           ! Water temperature:
           lkalpha=1./(6.7-0.829*logD(ireach(n)))
           ftempold(n)=lkalpha*Temp_avg(n)+(1-lkalpha)*ftempold(n)
           tempw(n)=0.62*logD(ireach(n))+.979*ftempold(n)
     *          +(0.0126-0.0059*logD(ireach(n)))*Kso
           if(tempw(n)<-0.5)tempw(n)=-0.5
           
           ! Net radiation:
           if(tempw(n).eq.-0.5)then
             short_wave=tau*Ket*(1-0.8)
           else
             short_wave=tau*Ket*(1-(albedo(month_now,9-int(ilat(n)/10))
     *             *(ilat(n)/10-int(ilat(n)/10))+albedo(month_now,10-
     *            int(ilat(n)/10))*(1-ilat(n)/10+int(ilat(n)/10)))/100.)
           endif
           long_wave=(5.67*10**-8.)*((Temp_avg(n)+273.16)**4.
     *          *(1-0.261*exp(-0.00077*Temp_avg(n)**2.))
     *          -0.98*(tempw(n)+273.16)**4.)
           net_solar_radv=short_wave+long_wave
             
           ! Heat storage flux:
           heatflux=-23+0.232*Ket+28.2*(Temp_avg(n)-tempw(n))
     *          -2.1*Temp_avg(n) !*log(LKdepth(ireach(n)))-1.1*Temp_avg(n)
           if(heatstore(n)<-heatflux) heatflux=-heatstore(n)
           heatstore(n)=heatstore(n)+heatflux
      
           ! Evaporation rate (daily)
           evapE=-heatflux
           
           if(net_solar_radv>0)evapE=evapE+net_solar_radv
           if(evapE.lt.0.0)evapE=0.0

             gam=.001013/0.622*(100*((44331.514-elev(n))
     *              /11880.516)**(1/0.1902632)/1000)/alamb
             delt=2508.3/(tempw(n)+237.3)**2*exp
     *           (17.3*tempw(n)/(tempw(n)+237.3))
             ddg=delt/(delt+gam)
        
            ! calculate lake evaporation (in mm/hr)
             evap_rate(n)=alpha*ddg*(evapE)/alamb/den*86.4/24
             
          if(heatstore(n)<0) heatstore(n)=0.0   ! fixed mixed mode
         
         ! Start up the temperature tracking for the next day
         Temp_tot(n)=tempv(n)
         Temp_max(n)=tempv(n)
         Temp_min(n)=tempv(n)
         
        end if !end of day if
             
        strloss(n)=evap_rate(n)*tdum(n)*(1.0-sca(n,ii_water))
        
!       use ftall(water) for optimization     
        strloss(n)=strloss(n)*ftall(ii_water)  ! dds optimization  
        
        if(r(n,classcount-1).gt.0.0) strloss(n)=0.0   ! fixed mixed mode
        
        if(strloss(n).lt.0.0) strloss(n)=0.0  ! fixed mixed mode
       
        evt(n,classcount-1)=evt(n,classcount-1)+evap_rate(n) 
        eloss(n)=eloss(n)+evap_rate(n)*aclass(n,classcount-1)
      
       else              ! end lake if                   PET method:
!         SURFACE WATER:
!         RAIN FALLS ON WATER SURFACE AND IS DIRECTLY ADDED TO 
!         RIVER FLOW. 
!         REV 7.9 NEW CALCULATION BASED ON NEW EVAPORATION
!         REV. 8.85 - Oct. 12/98 - FIXED RAIN & SNOW ON WATER CLASS
!         qstream(n) is added to channel inflow in route
!         qstream(n) = net precip in mm converted to cms
!          qstream(n)=r(n,ii_water)*tdum(n)        
!         if the water evaporation is read in, skip this part

!         For this option:
!         EVAPORATION IN THE WATER CLASS IS EQUAL TO
!         POTENTIAL EVAPORATION AND HAS TO BE TAKEN FROM
!         RIVER STORAGE AS A FLOW
!     rev. 9.1.66  Oct.  17/04  - NK; pet*fpet for loss from water instead of pet
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!change  nov 25/09 nk
          if(r(n,ii_water).gt.0.0)then
            ev(n,ii_water)=0.0    
          else
            ev(n,ii_water)=pet(n,ii_water)*fpet(ii_water)   ! pet calculated 
          endif
                
c!     REV. 10.1.16 Jan.  11/16  - NK: Added subroutine ice_factor.f
          if(ice_fctr(n).lt.0.0)then
            strloss(n)=0.0
            ev(n,ii_water)=0.0
          else !  ice_fctr(n)>.0.0

!     rev. 9.5.61  Sep.  03/09  - NK: bug/eloss - added water class for wfo weighted et
!           for the wfo file:

!            if(r(n,ii_water).le.0.001)then
              strloss(n)=ev(n,ii_water)*tdum(n)*(1.0-sca(n,ii_water))
!              strloss(n)=strloss(n)*fpet(ii_water)  ! dds optimization TH: fpet already in ev
!            else
!             no evaporation when it's raining            
!              strloss(n)=0.0
!            endif
            
!           strloss is taken from channel storage in route
c            if(ireach(n).eq.0)then
              if(store2(n).gt.strloss(n)*3600.0)then
                evt(n,ii_water)=evt(n,ii_water)+ev(n,ii_water)         ! for rffxx.txt
!               for ensim: prorated for sca 
!               100% snow cover - no et
!                0% snow cover - no reduction
                sum_et(n,ii_water)=sum_et(n,ii_water)+
     *                      ev(n,ii_water)*(1.0-sca(n,ii_water))   
     
!               for the watflood.wfo and watbal2.csv files:     
                eloss(n)=eloss(n)+ev(n,ii_water)*aclass(n,ii_water)         ! *(1-sca(n,ii_water))

!     REV. 10.1.17 Jan.  11/16  - NK: Added fpetLakeOverride factor
!               Temporary fix for dealing with rogue lakes - i.e. where fpet
!               just doesn't fit all the lakes
!               Maybve the lake evaporation model will fix this - maybe not.
!               It has to be done here so the tracer & isotope stuff 
!               is not affected
              else
                strloss(n)=0.0
              endif     !  store2(n).gt.............
c            endif     !  ireach(n).eq.0
          endif     !  ice_fctr(n).lt.0.0

   
!          snow covered ice - eventually.........................

c            else                    !IF(akfs(ii).gt.0.0)THEN
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

      end if !  ireach(n)>0.and.LKdepth(rbin)>1.0

      end do !end grid loop  

      firstpass=.false.
      RETURN
49900     format((i5,','),20(f10.5,','))     
 9992 format(7x,f9.2)
      END SUBROUTINE lake_evap
