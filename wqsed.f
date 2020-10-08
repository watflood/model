!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      subroutine wqsed(jan,aintvl) 
!***********************************************************************
!    Copyright (C) 1987 by Nicholas Kouwen and Luis Leon
        
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
     
! ************************************************************************ 
! modified: May/97 - Luis Leon
!          This subroutine calculates the sediment concentration 
!          for each element of the watershed yrot(n)
    
      use area_watflood
	implicit none

      integer :: jan,ii,n,i,j,i3,ii1,ii2
      real ::    aintvl,amin1,log10,psi  !<< replaced phi



! ****** Equations References And Units From:
!       Hartley, D.M.  1987.  "Simplified Process Model For Water
!            Sediment Yield From Single Storms.  Part I - Model
!            Formulation"  In Transactions Of The Asae 30(3):710-717.
!
!       naa = number of "n" elements
!       classcount = total number of "ii" land classes
!       gamma = specific weight of water [kg/m2/s2]
!       ro = density of the fluid [kg/m3]
!       viskin = kinemative viscosity [m2/s]
!       grav = standard acceleration [m/s2]
!       a,b = empirical constants used in transport capacity formula
!       kf(ii) = overland flow resistance [-]
!       gc(ii) = ground cover in percent
!       cf(ii) = cover factor in percent - determined by fig. 5.(hartley)
!       erod(n) = soil erodibility d [g/j]
!       d50(n) = average d50 of the detached sediment [mm] (lookup table)
!       spg(n) = specific gravity of detached sediment [-] (lookup table)
!       rey = shield's criterion reynold's number [-]
!       psi = shield's parmeter [-]
!       slope(n) = sqrt of channel slope
!       sl1(n) = sqrt of overland slope
!       al = length of element [m]
!       frac(n) = fraction of element within the basin
!       aclass(n,ii) = fraction of class ii in n
!       d1 = depth of water [mm]
!      ds = depression storage [mm]
! *note: hsed(n,ii) = d1-ds, average end depth [mm] {calc in runof5}
! *note: qs(n,ii) = overland flow [m3/s] = q1(ii) {calc in runof5}
!      ql(n,ii) = unit overland flow [m2/hr]
!      flow_shear = flow shear stress [units of gamma * mm]
!      crit_shear = critical shear [units of gamma * mm]
!      c_val = volumetric sediment concentration [-]
! *note: rf = runoff [mm] : rf(n,ii) = q1(ii)/tdum {calc in runof5}
!      y_crit = sediment transport capacity [kg/m2]
!      coef_cover = coefficient group for eq.44
!      ernfl = rate of rainfall energy [J/m2/hr]
!      ero = rate of runoff energy [J/m2/hr]
!      grf = rainfall soil detachment [kg/m2/hr]
!      gro = runoff soil detachment [kg/m2/hr]
!      y_pot = potential sediment supply/yield [kg/m2]
!      y(n,ii) = mininum of y_pot and y_crit [kg/m2]
!      yrot(n) = total yield from all classes for routing [kg/m3]
! 
!  ****** initialize values  ********************

       if(jan.eq.1)then

          do 40 ii=1,classcount
!            equation 16 (for each landclass)
             kf(ii)=60.0+3140.*gc(ii)**1.65
   40     continue
       endif

      do 50 n=1,naa
          yrot(n)=0.0
          do 50 ii=1,classcount
              y(n,ii) =0.0
              if(hsed(n,ii).le.0.0)then
                 hsed(n,ii)=0.0
              endif
!      convert to unit flow discharge (aintvl in hrs)
              ql(n,ii)=(qs(n,ii)/al)*(aintvl*3600.)
   50  continue
              
! ****** calculate rainffall energy *************

       do 60 n=1,naa
         i=yyy(n)
         j=xxx(n)
          if(p(i,j).le.0.0)then
             ernfl(n)=0.0
           else
!            eq43 [J/m2/hr] (for time step=1hr, p=rain intensity
!                 if not then divide p by time step)
             ernfl(n)=p(i,j)*(11.9+8.7*log10(p(i,j)))
           endif
   60  continue    

! ****** for each element ***********************

      do 120 n=1,naa
      if(slope(n).gt.0.0) then
          i=yyy(n)
          j=xxx(n)
          i3=ibn(n) 
          if(classcount.le.0) then
!           if classcount =0 then program runs in the lumped mode
            ii1=ibn(n)
            ii2=ibn(n)
          else
!           parameters grouped by land cover/use
            ii1=1
            ii2=classcount
          endif

! ****** for each land class type ***************
    
         do 110 ii=1,classcount
! *note: already change of qs [m3/s] to ql [m2/hr]
           if ((qs(n,ii).gt.0.0).and.(hsed(n,ii).gt.0.0))then
         
! * * * *  calculate the sediment tranport capacity
!             eq27 [gamma*mm] gamma cancels in eq21 : 5/8 is beta/beta+1
                flow_shear=(5.*60./8.*kf(ii))*hsed(n,ii)*sl1(n)
!             eq29 [-] including grav because shear has no gamma
                rey=((flow_shear*grav/1000.)**0.5)*d50(n)/viskin
!             from equation of shields diagram
                psi=(0.11/rey)+0.021*log10(rey)
!             eq20 [gamma*mm] gamma cancels in eq21
                crit_shear=(spg(n)-1.0)*psi*d50(n)
!             eq21 [-] for volumetric concentration
                c_val=a_wq*(flow_shear/crit_shear)**b_wq
!             eq42 [kg/m2] to compare with potential yield
                y_crit=2.65*ro*c_val*rf(n,ii)/1000.
                
!     changed a & b to a_wq and b_wq

! * * * *  calculate the sediment supply
!             eq44 [kg/m2/hr] divide by 1000 to convert g to kg
                coef_cover=(1-gc(ii))*cf(ii)
                grf=(ernfl(n)*coef_cover*erod(n))/1000.
!             eq45 [j/m2/hr] warning! the constant 60 has units
                ero=(60/kf(ii))*gamma*(ql(n,ii)/2.0)*sl1(n)
!             eq46 [kg/m2/hr] divide by 1000 to convert g to kg
                gro=(ero*erod(n))/1000.
!             eq47 [kg/m2]
                y_pot=(grf+gro)
                
! * * * *  calculate the sediment yield
!             compare capacity with supply and choose minimum
                y(n,ii)=amin1(y_crit,y_pot)
!             class weighted and converted to [kg/m3] /hsed in m (div/1000)
                y(n,ii)=y(n,ii)*frac(n)*aclass(n,ii)/(hsed(n,ii)/1000.) 
!             add for all classes to sediment supply per cell                
                yrot(n)=yrot(n)+y(n,ii)           

         endif

!      if(yrot(n).ne.0.0)then
!      print*,n,ii,yrot(n)
!      endif
  


110      continue
       endif
120    continue

!                   * * * end of main loop * * *
!  yrot(n)=total sediment inflow from element (n) in [kg/m3]
!  and ready to be routed into the downstream element ==> route

      return

      end
