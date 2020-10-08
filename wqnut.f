!************************************* 
      subroutine wqnut
! ************************************************************************ 

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
     
!     called every hour

! created Jul/98 - Luis Leon
!          this subroutine calculates the nutrient concentrations 
!          for each element of the watershed: cron(n,ii), crop(n,ii)

      use area_watflood
	implicit none

      integer :: n,ii,i,j
      real :: fpor,soln,pavr,peff,pavs,term1n,termexp1n,
     *        termexp2n,term2n,term3n,term1p,termexp1p,
     *        termexp2p,term2p,term3p,solp
       
! ***** Equations, References And Units From The Creams Model:
!      Frere Et.Al. (1980) The Nutrient Submodel. In: Knisel, W.G.
!           Creams: A Field Scale Model For Chemicals, Runoff, And Erosion
!           From Agricultural Management Systems, Usda, Rep No. 26
!      Young E.Al. (1986) Agricultural Nonpoint Source Pollution Model:
!           A Watershed Analysis Tool, Model Documentation, Usda.
!
! [general]
!      naa = number of "n" elements
!      classcount = total number of "ii" land classes
!      ieff = effective infiltration for storm [mm] - watflood value fake(ii)??
!      roff = total runoff for the storm [mm] - watflood value hsed(n,ii) = (d1-ds)
!      fpor = porosity factor [-]
!      peff = effective precipitation [mm]
!      por =  soil porosity [-]
!      spg(n) =  soil specific weight [-]
!      p(i,j) = storm precipitation for the time step [mm] - watflood value
!      er = nuetrient enrichment ratio
!      ysed = total sediment yield [kg/ha]
!      nofer = number of cells with fertilizer application
! [nitrogen]
!      cron(n,ii) = soluble nitrogen concentration in the runoff [kg/ha]
!      cronrot(n) = nitrogen concentration all classes for routing [kg/ha]
!      navs = available nitrogen in the surface [kg/ha]
!      navr = available nitrogen due to rainfall [kg/ha]
!      ndmv = rate for downward movement of nitrogen into the soil [1/mm]
!      nrmv = rate for nitrogen movement into the runoff [1/mm]
!      nrnc = nitrogen contribution due to rain [kg/ha]
!      soln = soluble nitrogen in the surface cm of the soil [kg/ha]
!      nfer(n) = nitrogen fertilizer application [kg/ha] - input data
!      nfa(n)  = fraction of nitrogen availability [%/100] - input data
!      ncpw = nitrogen concentration in pore water [ppm] - input data
!      ncrn = nitrogen concentration in rainfall [ppm] - input data
!      nlec = nitrogen leaching extraction coefficient - input data
!      nrec = nitrogen runoff extraction coefficient - input data
!      ndec = nitrogen decay fraction [%]
!      nsed = overland nitrogen transported by sediment [kg/ha]
!      nscn = soil nitrogen concentration [g n/g soil] - input data
! [phosphorus]
!      crop(n,ii) = soluble phosphorus concentration in the runoff [kg/ha]
!      croprot(n) = phosphorus concentration all classes for routing [kg/ha]
!      pavs = available phosphorus in the surface [kg/ha]
!      pavr = available phosphorus due to residual in soil [kg/ha]
!      pdmv = rate for downward movement of phosphorus into the soil [1/mm]
!      prmv = rate for phosphorus movement into the runoff [1/mm]
!      solp = soluble phosphorus in the surface cm of the soil [kg/ha]
!      pfer(n) = phosphorus fertilizer application [kg/ha] - input data
!      pfa(n)  = fraction of phosphorus availability [%/100] - input data
!      pcpw = phosphorus concentration in pore water [ppm] - input data
!      plec = phosphorus leaching extraction coefficient - input data
!      prec = phosphorus runoff extraction coefficient - input data
!      pdec = phosphorus decay fraction [%]      
!      psed = overland phosphorus transported by sediment [kg/ha]
!      pscn = soil phosphorus concentration [g p/g soil] - input data

!      ** initialize values      ****************

       do 10 n=1,naa
         cronrot(n)=0.0
         croprot(n)=0.0
         do 10 ii=1,classcount
              cron(n,ii)=0.0
              crop(n,ii)=0.0
   10  continue
 
   
!     ** calculate constants first **************

!      calculate the available nitrogen due to rainfall
      navr=ncrn*0.000001  
    
      do 20 n=1,naa
       i=yyy(n)
       j=xxx(n)
         if (p(i,j).gt.0.0) then
!         calculate soil porosity and porosity factor:
                por=1.0-(spg(n)/2.65)
                fpor=0.00001/por
!         calculate nutrients movement rates:
                ndmv=nlec/(10.*por)
                nrmv=nrec/(10.*por)
                pdmv=plec/(10.*por)
                prmv=prec/(10.*por)
!         calculate soluble nutrients in top cm of soil:
                soln=0.10*ncpw*por
                solp=0.10*pcpw*por
!         calculate available phosphorus due to soil residual:
                pavr=solp*fpor
!         calculate nitrogen contribution due to rainfall:
                nrnc=ncrn*p(i,j)
!         calculate the effective precipitation in the top cm:
                peff=p(i,j)-(10.*fpor)                
!         calculate the available nutrients in the surface:
                nfa(n)=nfa(n)/100.
                navs=(soln+(nfer(n)*nfa(n)))*fpor
                pfa(n)=pfa(n)/100.
                pavs=(solp+(pfer(n)*pfa(n)))*fpor

!     ** for each land class type ***************

      do 30 ii=1,classcount
           if ((hsed(n,ii).gt.0.0))then
!         calculate the nitrogen concentration in runoff
                term1n=(navs-navr)/fpor
                termexp1n=exp(-ndmv*fake(ii))
                termexp2n=exp(-ndmv*fake(ii)-nrmv*hsed(n,ii))
                term2n=termexp1n-termexp2n
                term3n=nrnc*hsed(n,ii)/peff
!         divide by 10 for effective precipitation in top centimeter                
                cron(n,ii)=((term1n*term2n)+term3n)/10.
!         calculate the phosphorus concentration in runoff
                term1p=(pavs-pavr)/fpor
                termexp1p=exp(-pdmv*fake(ii))
                termexp2p=exp(-pdmv*fake(ii)-prmv*hsed(n,ii))
                term2p=termexp1p-termexp2p
                term3p=pavr*prmv*hsed(n,ii)/fpor
!         multiply by 10 for porosity factor in top centimeter                    
                crop(n,ii)=((term1p*term2p)+term3p)*10.
!         class weighted by fractions of landclass (kg/ha)
                cron(n,ii)=cron(n,ii)*frac(n)*aclass(n,ii)
                crop(n,ii)=crop(n,ii)*frac(n)*aclass(n,ii)
!         conversion to kg/m3 (div/10,000) and then by hsed in m (/1000)
                cron(n,ii)=(cron(n,ii)/10000.)/(hsed(n,ii)/1000.)
                crop(n,ii)=(crop(n,ii)/10000.)/(hsed(n,ii)/1000.)
!         add for all classes to nutrient concentration per cell                
                cronrot(n)=cronrot(n)+cron(n,ii)           
                croprot(n)=croprot(n)+crop(n,ii)
           endif
30    continue


           endif
20    continue


! ****** for each element ***********************
!                  ***end of main loop** *
!  cronrot(n) and croprot(n)=total nutrients from element (n) in [kg/m3]
!           and ready to be routed into the downstream element.

      return
      
      end
