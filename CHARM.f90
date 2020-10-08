
        PROGRAM CHARM   ! Canadian Hydrological And Routing Model

!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen  
        
!    This file is part of WATFLOOD (R)      
        
!    WATFLOOD(R) is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.


!    WATFLOOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty ofdir runreport
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.

!    You should have received a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
     
!***********************************************************************
!    This program is the Canadian Hydrological And Routing Model (CHARM (R) )
!    It's purpose is to ingest precipitation and temperature time series data
!    and compute hydrographs at selected locations.
!***********************************************************************

        

!           the use of the main program is only to dimension the
!           variables and to decide a few options.
!
!           note that the dimensions of a variable must be the
!           same in all parts of the programs. eg: if qsyn is
!           dimensioned to (1,8,123) in the main calling program,
!           then is it is dimensioned qsyn(ls,js,ih), ls must be
!           set =1, js=8, & ih=123
!
!           the variables qi1 and qi2 have to be dimensioned to
!           the value of na while all the other variables can be
!           dimensioned to naa.  this is to allow the water to
!           run into something at the outlet.

    USE area_watflood
    use areacg
    use area_debug
!///////////////////////// 
!// Added by Dave
    USE EF_module
!// End Dave addition
!/////////////////////////

!     rev. 9.5.44  Oct.  27/08  - NK: removed code & obj modules for hasp & rainbow

      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

! TW  MAR. 5/98  - DATA STATEMENT MOVED FROM INPEVT.FOR

      CHARACTER*14    :: date
    character(3)    :: beepflg
      CHARACTER(1)    :: recflg
    character(4096)  :: line

      integer     :: n,ii,i,j,mon
      integer    ::  iallocate,attCount,ios,nhr,nhf,iendarg,ix
      real*4     ::  conv,scale,smc5(16),e1,recmax
      character(20) :: junk
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
      INTEGER(2) status1
      CHARACTER(1) buf,stopflg,ensimopenflg
      logical       :: exists
      logical       :: foundEndHeader

!      logical, parameter :: debug_output = .false. 
!!DIR$         IF DEFINED (DEBUG)
!                 debug_output = .true.
!!DIR$         ENDIF
 
                
!     rev. 7.2      sept. 19/94 - added ireach(n) for dwoper input 
!     rev. 7.3      dec.  20/94 - added uz & lz drainage in runof4
!     rev. 7.31     jan.  08/95 - set record length for 40 flow sta 
!     rev. 7.31.1   jan.  08/95 - set met data source for lapse rate
!     rev. 7.32     feb.  07/95 - added nopt to select opt flow sta 
!     rev. 7.33     feb.  20/95 - fixed flow initialization
!     rev. 7.4      feb.  24/95 - added 4 classes - max = 10
!     rev. not completed
!     rev. 7.41     apr.  15/95 - calc strmfl output /w inp fmt
!     rev. 7.42     may.  15/95 - check for div. by 0 in runof4
!     rev. 7.5      seperate snow covered and bare ground
!                   modified for separation of snowcovered ground and 
!                   bare ground by Frank Seglenieks  Feb/1995  new 
!                   runof5 debugged and intergrated by NK July/1995
!     rev. 7.51     oct.  08/95 - revise init channel flow in SUB
!     rev. 7.52     oct.  23/95 - check for opt constraints in main1
!     rev. 7.6      nov.  13/95 - added andrea's sediment routines
!     rev. 7.7      dec.  25/95 - added Allyson's Columbia routing
!     rev. 7.71     jan.  15/95 - fixed bug in uzs calculation 
!                                 uzs-retn =freely draining water
!     rev. 7.72     feb.  04/96 - took flowinit.for from sub.for
!     rev. 7.73     feb.  21/96 - fixed sca-continuity / runof5
!     rev. 7.74     may.  23/96 - include lapse rate & elv ref 
!                                 as part of .tmp file
!     rev. 7.75     may.  27/96 - added ak2fs in param & runof5
!     rev. 7.76     jun.  11/96 - # classes increased to 16 + urban
!     rev. 7.77     Jul.  02/96 - fixed snow redistribution
!     Rev. 7.78     Sept. 29/96 - fileio: modified for error checking 

!     rev. 7.80     Oct.  29/96 - CHARM7 added yymmdd.rin for res inflows
!                               - unit = 39   fln = 09
!     rev. 7.81     Nov.  07/96 - rdevt: added flags for stuff
!     rev. 7.83     Nov.  30/96 - fix div. by 0 - check - in lst.for
!     rev. 7.84     Dec.  16/96 - changed pmelt so that snowmelt only
!                                 occurs on snow covered area

!     rev. 8.0     Dec   18/96  - Added Todd Neff's evaporation
!     rev. 8.1     Feb.  15/97  - TBC & RSM (to be continued & resume) 
!     rev. 8.2     Feb.  15/97  - parameter selection for opt in main1
!     rev. 8.21    Mar.  15/97  - rain/snow choice tied to base temp  
!     rev. 8.22    Mar.  15/97  - glacier MF 2X when new snow=gone
!     rev. 8.23    Mar.  25/97  - fixed bug in route - keep qo2 for res
!     rev. 8.24    Apr.  07/97  - added glacier melt multiplier gladjust
!                               - used uzs-retn to determine freely 
!                                 draining water
!     rev. 8.25    May.  22/97  - fixed allocating the basin # in 
!                                 flowinit
!     rev. 8.3     May.  22/97  - added the simout/outfiles capability
!     rev. 8.31    June   3/97  - added initial uzs values in evap.par
!     rev. 8.32    June  13/97  - bypassed non-flagged parameters in OPT
!     rev. 8.4     July  16/97  - fixed melt routine and added init def
!     rev. 8.41    July  21/97  - added tipm to the optimization table
!     rev. 8.5     Oct.  09/97  - deleted the old interception stuff
!     rev. 8.51    Oct.  09/97  - fixed -ve qr() problem in runof5
!     rev. 8.52    Nov.  14/97  - replaced x4()= in runof
!     rev. 8.60    Nov.  14/97  - added sl2 to the interflow calculation
!     rev. 8.61    Dec.  12/97  - added contflg for statistics cont'n
!     rev. 8.62    Dec.  30/97  - fixed param s/r comb'd et & par flgs
!     rev. 8.70    Jan.  23/98  - added precip adjustment in rain.for
!     rev. 8.71    Feb.  24/98  - added evpflg2 to rdevt.for
!     rev. 8.72    Mar.   5/98  -tw: moved flgevp2 data statement to 
!                                CHARM.for
!     rev. 8.73    Mar.   1/98  - changed mhrd to mhtot in flowinit
!     rev. 8.74    Mar.  31/98  - reinvented fs stuff in opt
!     rev. 8.75    Apr.  27/98  - took da out of the resume file
!     rev. 8.76    May   26/98  - added precadj diagnostic to rain.for
!     rev. 8.77    June   1/98  - added sub-basin error calculation
!     rev. 8.78    July   7/98  - added scalesnw and scaletem to rdevt
!     rev. 8.79    July   7/98  - added 24 water survey format in strfw
!     rev. 8.80    July   9/98  - fixed precip shutdown after smearing
!     rev. 8.81    July  17/98  - precip adjust for T > 0 C only
!     rev. 8.82    July  10/98  - added runoff output option: routeflg
!     rev. 8.83    Sep.  23/98  - moved step args to area2.for
!     rev. 8.84    Sep.  28/98  - added runoff and evap fields to 
!                                 CHARM.txt
!     rev. 8.85    Oct.  12/98  - fixed rain & snow on water class
!     rev. 8.86    Nov.  02/98  - fixed opt problem found by ted.
!                               - fixed tto(n)=0 problem in etin
!     rev. 8.87    Nov.  17/98  - added watbal.for for water balance
!     rev. 8.88    Nov.  23/98  - fmadjust function of degree days
!     rev. 8.89    Nov.  30/98  - simplified uzs parameters
!     rev. 8.90    Dec.  04/98  - input to memory for opt runs
!     rev. 8.91    Dec.  07/98  - read rdevt in sub as well as CHARM!
!     rev. 8.92    Dec.  24/89  - check for 100% aclass coverage
!     rev. 8.93    Jan.  17/99  - sub modified for CHARM & watroute
!     rev. 8.94    Feb.  01/99  - crseflg to read resume & snow course
!     rev. 8.94a   Feb.  02/99  - reset heat deficit to 0.0 on Sept.01
!     rev. 8.94b   Feb.  06/99  - temperature correction and stop cmd
!     rev. 8.94c&d Feb.  20/99  - made paf.txt/error.txt default order
!     rev. 8.94e   Feb.  24/99  - added surfer output for error in lst
!     rev. 8.95    Mar.  15/99  - computed mean flows for time increment
!                               - involved getting rid of /kt throughout
!     rev. 8.96    Apr.  26/99  - lower zone function related to nbsn
!     rev. 8.96.1  May   12/99  - added ireport for reporting interval
!     rev. 8.97    July  12/99  - demonstration copy addition
!     rev. 8.98    July  15/99  - met grid shifting for weather models 
!     rev. 8.99    Aug.  18/99  - replaced err= with iostat= for f90
!     rev. 8.99a   Jul.     99  - lat-long watershed data
!     rev. 8.99b   Sept. 27/99  - divvy up interflow & drainae
!     rev. 8.99c   Oct.   5/99  - irough -> sl2 input in shed
!     rec. 8.99e   Nov.  29/99  - heat deficit initatialization
!     rev. 8.99f   Jan.   7/00  - changed uzs calcs re: shari's data
!     rev. 8.99g   Feb.   7/00  - added ttoinit to init evaporation
!     rev. 8.99k  feb. 15/2001  - fixex deficit calc in melt.for see9.06k
!     rev. 8.99l  Oct.    2001  - fixed reservoir release timing in CHARM8
!     rev. 8.99mm Dec. 13/2001-     added check for <= 0 init res flow
!     rev. 8.99n  Dec. 31/2001-     fixed nat. res initial flow (JW)
!     rev. 9.0     Mar.  21/00  - ts: converted to Fortran 90 
!                               - added dynamic memory allocation
!                               - added wfo file for ensim
!                  Fall   2000  - added wetland routing model         
!     rev. 9.01    Aug.   1/00  - added look up for minimum temperature
!                                 and function to calculate RH
!     rev. 9.02    Oct.   5/00  - added option to debug on one grid
!     rev. 9.03    Jan.   7/01  - set min precip rate for smearing
!     rev. 9.04    Jan    16/01 - fixed grid diagnosis in flowinit
!     rev. 9.05    Feb.   6/01  - chngd unit 61 to snw1.csv for surfer
!     rev. 9.06k   Feb.  15/01  - fixed deficit calc in melt (rem. qlz.txt) =8.99k
!     rev. 9.07    Mar.  14/01  - fixed use of opt par's  for numa=0  
!     rev. 9.08    Mar.  26/01  - checked limits on heat def.
!     rev. 9.08.01 Apr.   3/01  - check wetland designation in param
!     rev. 9.1     May    7/01  - updated Luis's sed & nutrient stuff
!     rev  9.1.02  July  12/01  - put in dacheck in flowinit for wetland flag
!     rev  9.1.03  July  24/01  - added polinomial to reservoir routing
!     rev. 9.1.04  Oct.   4/01  - added A7 for weighting old/new sca in melt
!                               - fixed Jan. 17/02 - didn't work before
!     rev. 8.99n   Dec.31/2001  -     fixed nat. res initial flow (JW)
!     rev. 9.1.05  Oct.   4/01  - new format parameter file
!     rev. 9.1.06  Oct.  16/01  - nrvr added to area3 to set # river types
!     rev. 9.1.07  Jan.   3/02  - check that outlet is in a lake
!     rev. 9.1.08  Jan.  17/02  - fixed rev. 9.1.04
!     rev. 9.1.09  Jan.  21/02  - fixed reservoir release timing in CHARM9 see8.99l
!     rev. 9.1.10  Jan.  29/02  - flow nudging added for nopt(l)=2
!     rev. 9.1.11  Feb.  07/02  - fixed bug in reservoir routing 
!     rev. 9.1.12  Mar.  15/02  - added xdelta and ydelta for ensim
!     rev. 9.1.13  Mar.  23/02  - fixed resv. timing, moved to beginning of dt
!     rev. 9.1.14  Mar.  24/02  - fixed wetland min time step & outflow
!     rev. 9.1.15  Apr.  02/02  - Luis' sediment stuff runs. Not checked with old version.
!     rev. 9.1.16  Apr.  03/02  - Added wetland conditional to select river w/wo wetland
!     rev. 9.1.17  May   05/02  - Some tidying up
!     rev. 9.1.18  Jun.  03/02  - Added sub-watershed modelling capability
!     rev. 9.1.19  Jun.  22/02  - Added A9 as the max heat deficit/swe ratio
!     rev  9.1.20  Jun.  25/02  - Added A10 as the power on the UZ discharge function
!     rev. 9.1.21  Jun.  28/02  - Added wetland storage & outflow to the wfo file
!     rev. 9.1.22  Jul.  22/02  - Added simout\error.r2s file for ENSIM_Hydrologic
!     rev. 9.1.23  Jul.  23/02  - Added control for nudging in event #1
!     rev  9.1.24  Sep.  11/02  - Added scaleallsnw to set snw scale in event 1
!     rev  9.1.25  Sep.  11/02  - Added A11 as bare ground equiv. vegn height  
!     rev. 9.1.26  Sep.  11/02  - fixed wetland evaporation re: uzsi
!     rev. 9.1.27  Sept. 19/02  - Added isbaflg
!     rev. 9.1.28  Sept. 19/02  - Added shedlfg to replace the bsnm.shd file
!     rev. 9.1.29  Nov.  07/02  - Changed the threshold flow values for error calculations
!     rev. 9.1.30  Nov.  08/02  - added q1, qint, drng & qlz to the wfo file
!     rev. 9.1.31  Nov.  13/02  - Fixed the wetland Q to account for wetland area
!     rev. 9.1.32  Nov.  20/02  - Fixed fpetmon() wrt. h()
!     rev. 9.1.33  Dec.  05/02  - Fixed instability in wetland flow    
!     rev. 9.1.34  Dec.  23/02  - Added ensim1flg - if ensimflg='a' for 1st id then 'y' for all events
!     rev. 9.1.35  Dec.  26/02  - Added wetland & channel heights to the wfo file
!     rev. 9.1.36  Jan.  28/03  - Fixed wetland init condition in flowinit
!     rev. 9.1.37  Mar.  22/03  - Option to turn off leakage by setting LZF < 0.0
!     rev. 9.1.38  Mar.  31/03  - revised str header and routing dt selectable
!     rev. 9.1.39  Apr.  06/03  - Fixed wetland routing when channel is dry
!     rev. 9.1.40  Apr.  24/03  - Min time step A6 read in strfw over rides the A6 from the par file
!     rev. 9.1.41  May   15/03  - Event average flows output to unit=75
!     rev. 9.1.42  May   31/03  - Tracer module added - first try
!     rev. 9.1.43  Jun.  01/03  - Fixed the qdwpr.txt function - re: last grid in lake
!     rev. 9.1.44  Jun.  11/03  - Added Cumulative precip to the wfo file
!     rev. 9.1.45  Jun.  11/03  - WATROUTE: runoff, recharge and leakage files added 
!     rev. 9.1.46  Jul.  17/03  - WATFLOOD LITE incorporated 
!     rev. 9.1.47  July  24/03  - TS: Tracer s/r deallocations added 
!     rev. 9.1.48  Dec.  08/03  - NK: sumrechrge() added to get total recharge
!     REV. 9.1.49  Nov.  23/03 - TS: Added wetlands to GW Tracer + Wetland Tracer
!     rev. 9.1.50  Jan.  14/04  - NK: version number added to the wfo_spec.txt file
!     rev. 9.1.51  Jan.  28/04  - NK: added iz.ne.jz conditional to ENSIM output  
!     rev. 9.1.52  Mar.  11/04  - NK: continuous water quality modelling
!     rev. 9.1.53  Mar.  14/04  - NK: hasp key configured
!     rev. 9.1.54  Apr.  12/04  - NK: SEDFLG set for multiple events at event No. 1
!     rev. 9.1.55  Jun.  12/04  - NK: write new str files to strfw\newfmt folder.
!     rev. 9.1.56  Jun.  18/04  - NK: write new rel & rin files to resrl\newfmt folder.
!     rev. 9.1.57  Jul.  06/04  - NK: Fixed major bug in shed.for max instead of min
!     rev. 9.1.58  Jul.  12/04  - NK: New header for the .shd file
!     rev. 9.1.59  Jul.  15/04  - NK: CHARMit rerout into two parts: rdresv & rerout
!     rev. 9.1.60  Jul.  27/04  - NK: reversed definitions for sl1 & sl2 Int. Slope
!     rev. 9.1.61  Aug.  25/04  - NK: Check for repeated met data in RAIN
!     rev. 9.1.62  Sep.  08/04  - NK: Fixed the conversion factor in SNW.FOR  (cnv)
!     rev. 9.1.63  Sep.  29/04  - NK: Added iopt_start as an arg for quick filecheck
!     rev. 9.1.64  Oct.  03/04  - NK: Coded up new header in ragmet.for
!     rev. 9.1.65  Oct.  03/04  - NK: Coded up new header for snow course file
!     rev. 9.1.66  Oct.  17/04  - NK; pet*ftall for loss from water instead of pet
!     rev. 9.1.67  Oct.  21/04  - NK; added unit 80 for lake_stor & lake_flow
!     rev. 9.1.68  Dec.  19/04  - NK: rewrote rdflow c/w memory allocation 
!     rev. 9.1.69  Dec.  19/04  - NK: rewrote rdresv c/w memory allocation 
!     rev. 9.1.70  Dec.  21/04  - NK: rewrote rdrain c/w memory allocation 
!     rev. 9.1.71  Dec.  28/04  - NK: rewrote rdtemp c/w memory allocation 
!     rev. 9.1.72  Dec.  28/04  - NK: fix bug in rdresv setting reach # 
!     rev. 9.1.73  Jan.  25/05  - NK: rewrote rdcrse c/w memory allocation 
!     rev. 9.1.74  Feb.  08/05  - NK: trashed rscrse replaced with rdswe
!     rev. 9.1.75  Feb.  08/05  - NK: added rdgsm (gridded soil moisture)
!     rev. 9.1.76  Mar.  09/05  - NK: separated glacier parameters in par file
!     rev. 9.1.77  Mar.  07/05  - NK: added .psm .gsm & .glz  files
!     rev. 9.1.78  Mar.  15/05  - NK: added WQD file to event file
!     rev. 9.1.79  Mar.  30/05  - NK: ktri to area2 for reservoir inflow dt
!     rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
!     rev. 9.1.81  Apr.  04/05  - NK: added sublimation,et and etfs to wfo file
!     rev. 9.2     Jun.  02/05  - NK: Numerous changes to program organization 
!     rev. 9.2.01  Jun.  29/05  - NK: Added write_r2s 
!     rev. 9.2.02  Jun.  29/05  - NK: Added read_r2s 
!     rev. 9.2.03  Jul.  11/05  - NK: Added s/r precip_adjust 
!     rev. 9.2.04  Jul.  13/05  - NK: allocation check for resrl 
!     rev. 9.2.05  Jul.  15/05  - NK: reversed order of reading resume file 
!     rev. 9.2.05  Jul.  27/05  - NK: initialized delta in s/r compute_error
!     rev. 9.2.06  Jul.  28/05  - NK: normalized error with da for optimization
!     rev. 9.2.07  Jul.  29/05  - NK: soilinit moved from runoff to sub 
!     rev. 9.2.08  Jul.  29/05  - NK: opt work-around in options 
!     rev. 9.2.09  Sep.  11/05  - NK: removed write_par.for from rdpar.for
!     rev. 9.2.10  Sep.  11/05  - NK: unlimited comments on .shd & .map files
!     rev. 9.2.11  Sep.  15/05  - NK: added Manning's n  r1n & r2n
!     rev. 9.2.12  Sep.  15/05  - NK: added EXCEL eqn to flowinit
!     rev. 9.2.13  Sep.  28/05  - NK: added freeze and break up to route
!     rev. 9.2.14  Sep.  29/05  - NK: Added control for opt in event #1
!     rev. 9.2.15  Sep.  30/05  - NK: Fixed bug for opt in flowinit
!     rev. 9.2.16  Oct.  10/05  - NK: Fixed bug for widep in rdpar
!     rev. 9.2.17  Oct.  11/05  - NK: Fixed bug for .str bounds in route
!     rev. 9.2.18  Oct.  27/05  - NK: Fixed bug in flowinit (init spike)
!     rev. 9.2.19  Oct.  28/05  - NK: Compute daily & monthly flows
!     rev. 9.2.20  Oct.  28/05  - NK: WFO_SPEC - reporting start & finish times 
!     rev. 9.2.21  Nov.  11/05  - NK: Set nopt in first event .str file 
!     rev. 9.2.22  Nov.  15/05  - NK: Fixed hmax bug in rdpar 
!     rev. 9.2.23  Nov.  22/05  - NK: Fixed res(n)=0 bug in route 
!     rev. 9.2.24  Dec.  07/05  - BT: DDS optimization 
!     rev. 9.2.25  Dec.  13/05  - NK: ENSIM r2c gridded soil moisture 
!     rev. 9.2.26  Dec.  23/05  - NK: Fixed reservoir outlet location bug 
!     rev. 9.2.27  Jan.  20/06  - NK: Separated header read in rdtemp
!     rev. 9.2.28  Jan.  30/06  - NK: Added low slope a4 for grids with water
!     rev. 9.2.29  Feb.  07/06  - NK: Read resv coeff first event only
!     rev. 9.2.30  Feb.  07/06  - NK: Added class_distribution.txt to output
!     rev. 9.2.31  Feb.  09/06  - NK: Added area chaeck to rdresume
!     rev. 9.2.32  Feb.  10/06  - NK: Added area_check.csv to output
!     rev. 9.2.33  Feb.  14/06  - NK: str stations from first event ONLY!!
!     rev. 9.2.34  Mar.  21/06  - NK: Activated glacier tracer1
!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
!     rev. 9.2.36  Mar.  30/06  - NK: Scaleallsnow changed to scale precip snow
!     rev. 9.2.37  Mar.  31/06  - NK: Removed impervious area as special class
!     rev. 9.2.38  Apr.  28/06  - NK: Lower bound set on a12 for smearing
!     rev. 9.2.39  May.  09/06  - NK: t added to route & rerout arg list
!     rev. 9.2.40  Jun.  09/06  - NK: added tto(),ttomin(),ttomax() to resume
!     rev. 9.2.41  Jun.  15/06  - NK: changed the resin.txt file to resin.csv
!     rev. 9.2.42  Jun.  20/06  - NK: water class included in the water balance
!     rev. 9.2.43  Jun.  21/06  - NK: fixed spikes in route
!     rev. 9.3.02  Jul.  18/06  - NK: converted runof, rchrg & lkage to r2c
!     rev. 9.3.03  Sep.  09/06  - NK: read s(i,j) from table instead of grid
!     rev. 9.3.04  Oct.  24/06  - NK: routing parameters dim to na in rte
!     rev. 9.3.05  Nov.  13/06  - NK: adder write_flowinit.for to flowinit.for
!     rev. 9.3.06  Dec.  17/06  - NK: added precip adjustment for bias
!     rev. 9.3.07  Dec.  29/06  - NK: added sum_precip for whole domain
!     rev. 9.3.08  Jan.  15/07  - NK: added lzs_init_new.r2c output to sub.for
!     rev. 9.3.09  Jan.  17/07  - NK: all file name lenghts = 60 in area12
!     rev. 9.3.10  Jan.  29/07  - NK: routing pars changed to gridded values
!     rev. 9.3.11  Feb.  28/07  - NK: ch_par added / event file ver = 9.5
!     rev. 9.4.01  Apr.  17/07  - NK: added deltat_report for gridflow.r2c
!     rev. 9.4.02  Apr.  18/07  - NK: moved rf, rffs from areawq to area1
!     rev. 9.4.03  Apr.  18/07  - NK: For water ev(n,ii)=pet(n,ii)*fpet(ii)
!     rev. 9.4.04  Apr.  23/07  - NK: moved allocate for melt from melt > CHARM
!     rev. 9.4.05  May.  04/07  - NK: revised timer for julian day calc.
!     rev. 9.4.06  May.  09/07  - NK: replaced por with spore(n,ii) in runof6
!     rev. 9.4.07  May.  15/07  - NK: converted opt to gridded routing parameters
!     rev. 9.4.08  May.  29/07  - NK: changed baseflow argument list
!     rev. 9.4.09  Jun.  19/07  - NK: added lake_area as a variable for iso
!     rev. 9.4.10  Jun.  19/07  - NK: adjusted frac for channel water area
!     rev. 9.4.11  Jun.  22/07  - NK: reordered rerout for glake 
!     rev. 9.4.12  Jul.  06/07  - NK: put qr + qstream - strloss back in runof6 
!     rev. 9.4.13  Jul.  09/07  - NK: modified lzs to account for lake area (flowinit) 
!     rev. 9.4.14  Jul.  09/07  - NK: added lake loss file 
!     rev. 9.4.15  Jul.  31/07  - NK: moved stuff from resume -> soil & flow init 
!     rev. 9.5     Sep.  07/07  - NK: changed wetland/channel routing 
!     rev. 9.5.01  Oct.  15/07  - NK: added wetland continuity check
!     rev. 9.5.02  Oct.  21/07  - NK: set init qdwpr=0.0 in route
!     rev. 9.5.03  Dec.  09/07  - NK: added reads for precip isotopes
!     rev. 9.5.04  Dec.  27/07  - NK: fixed bug in wetland routing
!     rev. 9.5.05  Jan.  13/08  - NK: added check for rec() in CHARM
!     rev. 9.5.06  Feb.  05/08  - NK: added pool and pool_o in rdpar & route
!     rev. 9.5.07  Feb.  05/08  - NK: fixed double counting of strloss & qstream
!     rev. 9.5.08  Feb.  08/08  - NK: new event parser
!     rev. 9.5.09  Feb.  12/08  - NK: added evap.r2c to the output files
!     rev. 9.5.10  Feb.  12/08  - NK: added water_area in lake_evap
!     rev. 9.5.11  Feb.  12/08  - NK: added -ve storage check for reservoirs
!     rev. 9.5.12  Feb.  13/08  - NK: added evaporation input file with read_r2c
!     rev. 9.5.13  Feb.  25/08  - NK: changed tolerance for coordinate check to .gt.0.001
!     rev. 9.5.14  Feb.  26/08  - NK: padded rel file for missing data
!     rev. 9.5.15  Feb.  28/08  - NK: fixed tdum & xdum for proper grid area in lat-long
!     rev. 9.5.16  Feb.  28/08  - NK: moved precip_adjust to sub
!     rev. 9.5.17  Feb.  28/08  - NK: moved scale snow from sub to process rain
!     rev. 9.5.18  Mar.  03/08  - NK: added conv to options & sub argument list
!     rev. 9.5.19  Mar.  05/08  - NK: prevented use of tracer * iso models with nudging
!     rev. 9.5.20  Mar.  06/08  - NK: added resvstore for iso model
!     rev. 9.5.21  Mar.  06/08  - NK: fixed dtmin for first time step each event
!     rev. 9.5.22  Mar.  12/08  - NK: added grdflg to print gridded flow, swe & evap
!     rev. 9.5.23  Mar.  12/08  - NK: fixed allocation error in read_resv_ef
!     rev. 9.5.24  Mar.  18/08  - NK: fixed missing data in read_resl_ef.f
!     rev. 9.5.25  Mar.  20/08  - NK: fixed lake initiation - moved code route -> flowinit
!     rev. 9.5.26  Apr.  04/08  - NK: added Julian day calc. to read_evt
!     rev. 9.5.27  Apr.  15/08  - NK: fixed allocation for chnl in rdpar
!     rev. 9.5.28  Apr.  15/08  - NK: fixed allocation for inbsnflg in flowinit
!     rev. 9.5.29  May.  26/08  - NK: fixed initialization in read_resv_ef
!     rev. 9.5.30  May.  26/08  - NK: conv back in read_rain & process_rain arg. list
!     rev. 9.5.31  May.  27/08  - NK: moved totsnw(n) computation in sub
!     rev. 9.5.32  Jun.  04/08  - NK: compute reservoir levels
!     rev. 9.5.33  Sep.  12/08  - NK: added column labels for grapher in flow_station_location.xyz
!     rev. 9.5.34  Sep.  17/08  - NK: fixed lake area in flowinit
!     rev. 9.5.35  Sep.  22/08  - NK: moved flow_sta_location to flowinit
!     rev. 9.5.36  Oct.  01/08  - NK: fixed ires bug for unevent dx & dy in read_resv
!     rev. 9.5.37  Oct.  14/08  - NK: added deltat_report to lake_sd.csv file write
!     rev. 9.5.38  Oct.  14/08  - NK: added optional coef6 & 7 to rel file for lake levels
!     rev. 9.5.39  Oct.  15/08  - NK: fixed bug in reservoit routing
!     rev. 9.5.40  Oct.  21/08  - NK: added diversions to rerout
!     rev. 9.5.41  Oct.  22/08  - NK: read in reservoir coefficients each event
!     rev. 9.5.42  Oct.  22/08  - NK: added b7() as the initial lake surface elevation
!     rev. 9.5.43  Oct.  27/08  - NK: changed bottom part of par file to be free format
!     rev. 9.5.44  Oct.  27/08  - NK: removed code & obj modules for hasp & rainbow
!     rev. 9.5.45  Dec.  16/08  - NK: added various error calculations - user's choice with errflg
!     rev. 9.5.46  Dec.  23/08  - NK: trying to fix problem with -ve storage. Changed conditional to .lt.
!     rev. 9.5.47  Dec.  26/08  - NK: add flwinitflg to warn about initial flows
!     rev. 9.5.48  Dec.  26/08  - NK: added event_fln() to allow unlimited events
!     rev. 9.5.49  Dec.  31/08  - NK: changed conditional to read releases in rerout
!     rev. 9.5.50  Jan.  05/09  - NK: read evap data for reaches only
!     rev. 9.5.51  Jan.  13/09  - NK: added reading yyyymmdd_ill.pt2 for all lakes
!     rev. 9.5.52  Jan.  20/09  - NK: added reading yyyymmdd_div.pt2 for diversions
!     rev. 9.5.53  Jan.  20/09  - NK: undid rev. 9.5.40
!     rev. 9.5.54  Feb.  11/09  - NK: undid rev. 9.2.28
!     rev. 9.5.55  Feb.  11/09  - NK: Correct R2n for instream lakes
!     rev. 9.5.56  Mar.  26/09  - NK: Fix bug with month in yearly events
!     rev. 9.5.57  Apr.  13/09  - NK: added ntrlflg for natural lake flows
!     rev. 9.5.58  Apr.  16/09  - NK: added nudgeflg for forcing gauge flows
!     rev. 9.5.59  Jul.  26/09  - NK: added fpet_lake for each lake in ill file
!     rev. 9.5.60  Sep.  01/09  - NK: added deltat_report for lake_sd.csv file
!     rev. 9.5.61  Sep.  03/09  - NK: bug/eloss - added water class for wfo weighted et
!     rev. 9.5.62  Sep.  04/09  - NK: new tb0 file for DW routing
!     rev. 9.5.63  Sep.  04/09  - NK: moved lapse rate from melt.f to process_temp.f
!     rev. 9.5.64  Sep.  16/09  - NK: corrected nudging wrt first event
!     rev. 9.5.65  Sep.  26/09  - NK: lapse rate changed from dC per 100 m to dC per m
!     rev. 9.5.66  Oct.  06/09  - NK: fixed bug in flowinit for init flows < 1.0
!     rev. 9.5.67  Oct.  06/09  - NK: fixed bug in rerout
!     rev. 9.5.68  Oct.  07/09  - NK: debugged read_resvin_ef.f
!     rev. 9.5.69  Oct.  10/09  - NK: added xcount & ycount to error & paf files
!     rev. 9.5.70  Oct.  11/09  - NK: fixed timer for r2c frames (use year_now)
!     rev. 9.5.71  Oct.  12/09  - NK: fixed bug in lst for setting value for nhyd(,)
!     rev. 9.5.72  Oct.  12/09  - NK: fixed bug in rdpar setting init values for fpet & ftal
!     rev. 9.5.73  Oct.  12/09  - NK: bypass using lake levels when optimizing
!     rev. 9.5.74  Oct.  21/09  - NK: in opt - made optim abs(optim)
!     rev. 9.5.75  Oct.  26/09  - NK: commented "deallocate in sub for watroute reads
!     rev. 9.5.76  Oct.  26/09  - NK: fixed basin exclusion for opt if resin present
!     rev. 9.5.77  Oct.  26/09  - NK: fixed some inits for out of basin gauges
!     rev. 9.5.78  Nov.  04/09  - NK: matched resvin locations to reach numbers
!     rev. 9.5.79  Nov.  04/09  - NK: added resumflg='s' for read_soilinit ONLY
!     rev. 9.5.80  Dec.  20/09  - NK: added swe_locations.txt file for swe input
!     rev. 9.5.81  Jan.  16/10  - NK: allow reservoirs outside watershed in resv file
!     rev. 9.5.82  Jan.  26/10  - NK: replaced error check for inflow locations
!     rev. 9.5.83  Feb.  17/10  - NK: non_basin exclusion for dds_flag=1
!     rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
!     rev. 9.6.01  Mar.  01/10  - NK: rlake parameter added for Manning n correction
!     rev. 9.6.02  Mar.  15/10  - NK: add sublimation to optimization
!     rev. 9.6.02  Mar.  23/10  - NK: add cumm_domain_precip
!     rev. 9.6.03  Mar.  31/10  - NK: replaced leakage.dat by nbs.tb0 fln(79)
!     rev. 9.6.04  Apr.  05/10  - NK: fixed filename carry over in read_evt
!     rev. 9.6.05  Apr.  06/10  - NK: added store_error_flag for -ve storage grids
!     rev. 9.6.06  Apr.  18/10  - NK: added glacier adjust for optimization
!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
!     rev. 9.7.01  Jun.  09/10  - NK: fixed error.xyz & error.r2s
!     rev. 9.7.02  Jun.  24/10  - NK: fixed bug in rdpar for classcount for imp area
!     rev. 9.7.03  Jun.  24/10  - NK: normalized SSE with station Qmean**2
!     rev. 9.7.04  Aug.  30/10  - NK: added to error message in read_rain & read_temp
!     rev. 9.7.04  Aug.  31/10  - NK: changed # decimal points for r2c files header
!     rev. 9.7.05  Aug.  31/10  - NK: changed error.r2s to error.r2c
!     rev. 9.7.06  Sep.  01/10  - NK: fixed subscript out of range errors in flowinit
!     rev. 9.7.07  Sep.  05/10  - NK: increased allowed # flow stations from 128 to 512
!     rev. 9.7.08  Sep.  21/10  - NK: revised mean squared error weighting for DDS
!     rev. 9.7.09  Sep.  29/10  - NK: corrected error.r2c file for sub-basin errors
!     rev. 9.7.09  Oct.  02/10  - NK: ensure fpet_lake is not assigned unintended values
!     rev. 9.7.10  Oct.  11/10  - NK: update flowflag in lst.f for subsequent events
!     rev. 9.7.11  Nov.  22/10  - NK: added monthly_climate_deltas.txt file
!     rev. 9.7.12  Nov.  10/10  - NK: fix array bugs for reservoir inflows
!     rev. 9.7.13  Nov.  22/10  - NK: Changed the outfiles.txt for more 30 rff classes
!     rev. 9.7.14  Nov.  22/10  - NK: Allow 30 land cover classes
!     rev. 9.7.15  Dec.  14/10  - NK: Create reduced precip & temp files for sub-basins
!     rev. 9.7.16  Jan.  05/11  - NK: Fixed init flows outside sub-basin
!     rev. 9.7.17  Jan.  05/11  - NK: Fixed diversions outside sub-basin
!     rev. 9.7.18  Jan.  17/11  - NK: Changed tolerance on the grid check in read_rain & read_temp
!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!     rev. 9.7.20  Jan.  31/11  - NK: Moved open statement for rdpar to rdpar/f
!     rev. 9.7.21. Mar.  07/11  - NK: Fixed delta_reort for longer periods in lst
!     rev. 9.7.22. Mar.  07/11  - NK: Changed diversion code: give/route take/rerout
!     rev. 9.7.23  Mar.  18/11  - NK: Revamped auto hydrograph fitting with precip icase=-2
!     rev. 9.7.24  Apr.  20/11  - NK: Added diverflg to indicate if a diversion is in grid
!     rev. 9.7.25  Apr.  28/11  - NK: Fixed daily flows
!     rev. 9.7.27  May.  26/11  - NK: Add lake_ice_factor
!     rev. 9.7.28  Jun.  14/11  - NK: Add degree_day for lake_ice_factor  dd_ice
!     rev. 9.7.29  Jul.  07/11  - NK: Add sublim_rate to set sublimation rate/day to par file
!     rev. 9.7.30  Jul.  13/11  - NK: imax > ycount & jmax > xcount also imin > 1 jnim > 1
!     rev. 9.8.00  Jul.  14/11  - NK: ntype+1 replaced by classcount (plus all derivatives)
!     rev. 9.8.01  Jul.  21/11  - NK: added ragmet optimization to dds setup
!     rev. 9.8.02  Jul.  26/11  - NK: reactivated meander length
!     rev. 9.8.02  Aug.  02/11  - NK: added lake level tb0 file
!     rev. 9.8.03  Aug.  08/11  - NK: check no of mean observed flows in file are ok
!     rev. 9.8.04  Sep.  02/11  - NK: Fix bug in write)par_10 when reading old par file
!     rev. 9.8.05  Oct.  18/11  - NK: New read_par_parser subroutine
!     rev. 9.8.06  Nov.  08/11  - NK: Added check for `water` class name
!     rev. 9.8.07  Oct.  10/11  - NK: area_check - removed unused stations
!     rev. 9.8.08  Nov.  18/11  - NK: added fratio for interception hight optimization
!     rev. 9.8.09  Nov.  22/11  - NK: nopt(l)=0 for area_error(l) > 10%
!     rev. 9.8.10  Dec.  06/11  - NK: Added message for FP overflow in route
!     rev. 9.8.11  Dec.  06/11  - NK: removed 30 char limit on find filetype 
!     rev. 9.8.12  Dec.  07/11  - NK: removed 30 char limit on find filetype 
!     rev. 9.8.12  Dec.  08/11  - NK: recognize kenueflg in the event file 
!     rev. 9.8.13  Jan.  17/12  - NK: modifications to read_r2c for single frame data
!     rev. 9.8.14  Jan.  27/12  - NK: dds_penalty added for swe not to zero in summer
!     rev. 9.8.15  Mar.  12/12  - NK: write error.txt for every dds evaluation
!     rev. 9.8.16  Mar.  21/12  - NK: reinstate reservoir inflow error for dds
!     rev. 9.8.17  Apr.  24/12  - NK: Moved dds flags to top of par file
!     rev. 9.8.18  Apr.  26/12  - NK: Added in-basin check in tracer4
!     rev. 9.8.19  May.  10/12  - NK: Added check on mising init flow for lakes
!     rev. 9.8.20  MAy.  15/12  - NK: fixed lake area in flowinit9.5.34
!     rev. 9.8.21  Jun.  18/12  - NK: Added swe observed date & report
!     rev. 9.8.22  Jul.  17/12  - NK: Added resetflg to reset cumm. precip Sept.1
!     rev. 9.8.23  Aug.  03/12  - NK: Added resinid1flg to use resinflg for id=1
!     rev. 9.8.24  Aug.  07/12  - NK: Added reading yyyymmdd_lvl.tb0 for lake levels
!     rev. 9.8.25  Sep.  26/12  - NK: Added warning for resumflg=y and ID > 1
!     rev. 9.8.26  Sep.  26/12  - NK: Added error check on # chained files for id>1'
!     rev. 9.8.27  Sep.  27/12  - NK: changed action on resumflg='s' - keep tbcflg='y'
!     rev. 9.8.28  Oct.  12/12  - NK: fixed heat deficit reset for resume
!     rev. 9.8.29  Oct.  15/12  - NK: added wetland_flag to speed up route.f
!     rev. 9.8.30  Oct.  16/12  - NK: remove p(i,j)=0.0 from precip_adjust
!     rev. 9.8.31  Oct.  16/12  - NK: continue rff files for contflg = y
!     rev. 9.8.32  Oct.  19/12  - NK: Fixed format for resin.csv in lst.f
!     rev. 9.8.33  Oct.  23/12  - NK: Deleted header for rff files with resumflg = y
!     rev. 9.8.34  Oct.  23/12  - NK: Added sums to the resume.txt file
!     rev. 9.8.35  Oct.  23/12  - NK: Fixed bug in read_soilinit_ef
!     rev. 9.8.36  Oct.  23/12  - NK: added fields to rff files
!     rev. 9.8.37  Oct.  27/12  - NK: added section to read_flow_ef to check # columns = no
!     rev. 9.8.38  Nov.  13/12  - NK: changed name level_plotting.xyz > level_station_location.xyz
!     rev. 9.8.39  Nov.  26/12  - NK: added check for flow stations in lakes
!     rev. 9.8.40  Jan.  14/13  - NK: convert interception cap: h(,)*fratio()
!     rev. 9.8.41  Jan.  28/13  - NK: fixed bug in lst for level print statement
!     rev. 9.8.42  Jan.  31/13  - NK: fixed bug in read_resvin: nopti int conversion
!     rev. 9.8.43  Jan.  31/13  - NK: fixed bug in lst.f : undefined output for iopt=99
!     rev. 9.8.44  Jan.  31/13  - NK: fixed bug in sub.f : uninitialized course_calc(n,j)
!     rev. 9.8.45  Jan.  31/13  - NK: disabled some writes for iopt = 99
!     rev. 9.8.46  Feb.  04/13  - NK: Fixed some write formats in lst,stats,watbal
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for CHARM & resin csv files
!     rev. 9.8.48  Feb.  12/13  - NK: Replaced CHARM.plt with CHARM.tb0 file
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals
!     rev. 9.8.50  Feb.  27/13  - NK: Initialize store1&2() for zero lake outflow
!     rev. 9.8.51  Mar.  11/13  - NK: Link skiphours in s/r stats to value1 in the str file
!     rev. 9.8.52  Mar.  20/13  - NK: deleted a pause for dds runs in route
!     rev. 9.8.53  Mar.  20/13  - NK: Add Lake St. Joseph diversion algorithm to REROUT.f
!     rev. 9.8.54  Apr.  02/13  - NK: deltat conversion seconds to hours
!     rev. 9.8.55  Apr.  10/13  - NK: fixed pause for dds runs in route
!     rev. 9.8.56  Apr.  10/13  - NK: Added check in rerout for -ve storage due to evaporation
!     rev. 9.8.57  Apr.  12/13  - NK: Added lakeEflg to stop lake evaporation whan levels very low
!     rev. 9.8.58  Apr.  12/13  - NK: REvised Family Lake (WPEGR) O/R in rerout
!     rev. 9.8.59  May   14/13  - NK: REmoved psmear & punused from the program
!     rev. 9.8.60  May   14/13  - NK: fixed ice factor for whole x-section
!     rev. 9.8.61  May   22/13  - NK: Introduced flag1 to speed up runof6
!     rev. 9.8.62  May   22/13  - NK: Fixed bug in runof6: (classcount-3) to (classcount-2) 
!     rev. 9.8.63  May   22/13  - NK: Fixed bug in s/r SUB.f argument list: "jan" missing 
!     rev. 9.8.64  May   28/13  - NK: Undocumented debug file 
!     rev. 9.8.65  May   28/13  - NK: Dimensioned firstpass_local()in REROUT
!     rev. 9.8.66  Jun   03/13  - NK: Added error_Dv.txt output in stats.f
!     rev. 9.8.67  Jun   06/13  - NK: Added allocation for flag1
!     rev. 9.8.68  Jun   17/13  - NK: Added dds_override file
!     rev. 9.8.69  Jun   17/13  - NK: Fixed bug in allocating clumnunits in SUB.f
!     rev. 9.8.70  Jun   17/13  - NK: for PAF: change error & PAF files to use GK formats
!     rev. 9.8.77  Jul   08/13  - NK: Made universal the use of wetland_flag(n)
!     rev. 9.8.78  Jul   16/13  - NK: Fixed divertflg to have the first event file value
!     rev. 9.8.79  Jul   19/13  - NK: Fixed wetland conditional screwed up with rev 9.8.77 in runof6
!     rev. 9.8.80  Aug   09/13  - NK: Added withdraw.r2c output file in route.f
!     rev. 9.8.81  Sep.  03/13  - NK: Add pafflg and update precip adjustment factors PAF!***********************************************************************
!     rev. 9.8.82  Sep.  07/13  - NK: Bypass of hard-coded lake rules when coeff1=0
!     rev. 9.8.83  Sep.  10/13  - NK: Set classcount=0 for fli.exe program only
!     rev. 9.8.84  Sep.  15/13  - NK: Added fratio to list of equal values for bog & fen
!     rev. 9.8.85  Sep.  30/13  - NK: Fixed the water balance for Lake St. Jo so diversion is taken care of
!     rev. 9.8.86  Oct.  16/13  - NK: Added version no to stats.txt output
!     rev. 9.8.87  Oct.  25/13  - NK: Added error message for mismatched resume file
!     rev. 9.8.88  Oct.  26/13  - NK: Fixed header writing sequence for CHARM.tb0 
!     rev. 9.8.89  Oct.  27/13  - NK: Fixed undefined (NAN) problem in flowint
!     rev. 9.8.90  Oct.  30/13  - NK: Added fetch to the shd file 
!     rev. 9.8.91  Oct.  30/13  - NK: Got rid of lzs_init.r2c - data is in flow_init.r2c already
!     rev. 9.8.92  Nov.  06/13  - NK: Changed output file swe.txt to swe.csv
!     rev. 9.8.93  Nov.  12/13  - NK: Added the routing initialization with yyyymmdd_fli.r2c
!     rev. 9.8.94  Nov.  20/13  - NK: Added check on interception capacity for water
!     rev. 9.8.95  Nov.  20/13  - NK: Changed unit 58 to 955 for CHARM.tb0
!     rev. 9.9.00  Dec.  08/13  - NK: Added Lake Evaporation model
!     rev. 9.9.01  Dec.  12/13  - NK: Added `pintwarning' in route added
!     rev. 9.9.02  Dec.  12/13  - NK: Changed format for origin in wfo code
!     rev. 9.9.03  Dec.  15/13  - NK: Change to gridded latitude for etharg
!     rev. 9.9.04  Dec.  17/13  - NK: Change over to gridded climate normals to diff
!     rev. 9.9.05  Jan.  02/14  - NK: Add check if in-basin in flowinit
!     rev. 9.9.06  Jan.  08/14  - NK: Add daily differences to Harfreaves ETHarg.f
!     rev. 9.9.07  Jan.  10/14  - NK: Overhaul of the frame numbers to EnSim specs
!     rev. 9.9.08  Jan.  10/14  - NK: Add check on diversion locations in read_divert'
!     rev. 9.9.09  Feb.  24/14  - NK: Fixed reading the time stame in r2c frame headers
!     rev. 9.9.10  Mar.  20/14  - NK: Update swe anytime a file is found
!     rev. 9.9.11  Mar.  20/14  - NK: Added lake_level_init.pt2 file for a resume
!     rev. 9.9.12  Apr.  04/14  - NK: Added min & max lake_level output file
!     rev. 9.9.13  Apr.  04/14  - NK: Fix water balance
!     rev. 9.9.14  Jun.  02/14  - NK: Fix water balance for water class
!     rev. 9.9.15  Jun.  02/14  - NK: Add lz to the water balance - it was missing
!     rev. 9.9.16  Jun.  06/14  - NK: Added location file for Root R. diversion
!     rev. 9.9.17  Jun.  07/14  - NK: Added check for allocation of outarray in Sub.f
!     rev. 9.9.18  Jun.  08/14  - NK: Fixed glacier_class check for wetlands
!     rev. 9.9.19  Jun.  11/14  - NK: Added a file for lat-long diversion locations for L. St. Jo
!     rev. 9.9.20  Jul.  15/14  - NK: Fix -ve lake storage when release data used
!     rev. 9.9.20  Jul.  24/14  - NK: Added dead storage for lakes "store_dead"
!     rev. 9.9.21  Jul.  27/14  - NK: Added allocation for outarray in sub
!     rev. 9.9.22  Jul.  29/14  - NK: Fixed basin no assignment in flowinit.f
!     rev. 9.9.23  Aug.  10/14  - NK: Uniform arguments for write_r2c
!     rev. 9.9.24  Aug.  20/14  - NK: Added monthly mean flow csv file CHARM_mly_nn.csv
!     rev. 9.9.25  Sep.  02/14  - NK: Finally fixed the error when nbasin=0
!     rev. 9.9.26  Sep.  16/14  - NK: Added precip adjust for forecast & fcstflg
!     rev. 9.9.27  Sep.  18/14  - NK: Added zero class bypass in intcept.f
!     rev. 9.9.28  Sep.  18/14  - NK: Added 'a' as option for ntrlflg & smrflg
!     rev. 9.9.29  Sep.  30/14  - NK: Remove unnecessary writes for watroute
!     rev. 9.9.30  Sep.  30/14  - NK: fixed allocation for qhyd_mly & qsyn_mly
!     rev. 9.9.31  Oct.  13/14  - NK: Changed flow initialization RE: zero init flows
!     rev. 9.9.33  Oct.  16/14  - NK: Added checks for files existing for a resume'
!     rev. 9.9.34  Oct.  17/14  - NK: Added re-compute of lake storage re: new lake levels
!     rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
!     rev. 9.9.36  Nov.  03/14  - NK: Revised error message for daily diff choices
!     rev. 9.9.37  Nov.  05/14  - NK: Added `newDataFlag` check to WATROUTE 
!     rev. 9.9.38  Nov.  12/14  - NK: Added LKdepth to ill file
!     rev. 9.9.39  Nov.  14/14  - NK: Modifications for watroute
!     rev. 9.9.40  Nov.  19/14  - NK: Modified the 'a' option for ntrlflg
!     rev. 9.9.41  Nov.  20/14  - NK: Added check if diversion = in-basin 
!     rev. 9.9.42  Nov.  26/14  - NK: Added errer check if diversion does not exist 
!     rev. 9.9.43  Nov.  26/14  - NK: Allocation for divertflg = 'g' 
!     rev. 9.9.44  Nov.  28/14  - NK: Added dead storage to reservoirs 
!     rev. 9.9.45  Dec.  03/14  - NK: Revamped read_pt2 for general use
!     rev. 9.9.46  Dec.  10/14  - NK: Added check on initial lake outflow
!     rev. 9.9.47  Dec.  24/14  - NK: Added lakeflg for lake evaporation option
!     rev. 9.9.48  Jan.  06/15  - NK: Added wetland cond. function for o/b flow
!     rev. 9.9.49  Jan.  06/15  - NK: Added courantflg
!     rev. 9.9.50  Jan.  07/15  - NK: Added zero - initial flow warning
!     rev. 9.9.51  Jan.  13/15  - NK: Added min channel area in flowinit
!     rev. 9.9.52  Jan.  14/15  - NK: Fixed bug for channel store < 0 for withdrawals
!     rev. 9.9.53  Jan.  18/15  - NK: Prevent mode switch during iteration in wetland routing
!     rev. 9.9.54  Jan.  19/15  - NK: Put par & shd file names for 1st event in the headers
!     rev. 9.9.55  Jan.  22/15  - NK: Added diversion upstream drainage area in div file
!     rev. 9.9.56  Feb.  04/15  - NK: Fixed missing initial rel data in read_resv
!     rev. 9.9.57  Feb.  08/15  - NK: Fixed resv inflow output resin & lake_sd
!     rev. 9.9.58  Feb.  13/15  - NK: Added time column to levels.txt
!     rev. 9.9.59  Mar.  06/15  - NK: In route: strloss option frcflg y/n
!     rev. 9.9.60  Mar.  06/15  - NK: In sub: fixed call write_r2c for close condition 
!     rev. 9.9.61  Mar.  06/15  - NK: In route: restored hcha2(n)=store2(n)/chaarea(n)
!     rev. 9.9.62  Mar.  21/15  - NK: Change zone from character to integer
!     rev. 9.9.63  Apr.  06/15  - NK: Changed reas_resv to carry on with lask known release(s)
!     rev. 9.9.64  Apr.  08/15  - NK: DDS bypass in sub for single runs
!     rev. 9.9.65  Apr.  03/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
!     rev. 9.9.66  Apr.  03/15  - NK: Added options to write_tb0 str files
!     rev. 9.9.67  Apr.  29/15  - NK: Deleted mid_file headers in with tbcflg=y
!     rev. 9.9.68  Apr.  29/15  - NK: Fixed tto reset with resume
!     rev. 9.9.69  Jun.  10/15  - NK: prevent write ro rff if there is no class area
!     rev. 9.9.70  Jun.  12/15  - NK: Add del_rain, and dSTRconc2  to the wfo file
!     rev. 9.9.71  Jun.  13/15  - NK: DDS obf function taken out of sub > s/r obj_fn
!     rev. 9.9.72  Jul.  21/15  - NK: Dave Newson additions to sub & process_rain
!     rev. 9.9.73  Aug.  31/15  - NK: Finshed rules s/r - ready for beta testing
!     rev. 9.9.74  Sep.  11/15  - NK: Added output to unit 53 in flowinit
!     rev. 9.9.75  Sep.  11/15  - NK: Added basin_no.r2c output to flowinit.f
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations
!     rev. 9.9.77  Sep.  11/15  - NK: S/r read_ts5 created
!     rev. 9.9.78  Sep.  16/15  - NK: Fixed wcl in melt.f
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
!     rev. 10.1.02 Oct.  09/15  - NK: Fixed allocation for qhyd_mly & qsyn_mly
!     rev. 10.1.03 Oct.  09/15  - NK: Added units 81-83 for isotope output
!     rev. 10.1.04 Oct.  10/15  - NK: Added year_last variable for use in reading isotope data
!     rev. 10.1.05 Oct.  11/15  - NK: Iso RMS error 
!     rev. 10.1.06 Nov.  19/15  - NK: Added area_check with can_discharge_sites.xyz 
!     rev. 10.1.07 Dec.  02/15  - NK: Added ice_fctr(n) to route 
!     rev. 10.1.08 Dec.  04/15  - NK: Added msg re: replacing "mean_observed_flows.txt"' 
!     rev. 10.1.09 Dec.  07/15  - NK: Add blank line for missing data in the precip.txt file in lst.f 
!     rev. 10.1.10 Dec.  09/15  - NK: Add blank line for missing data in the precip.txt file in lst.f 
!     rev. 10.1.11 Dec.  11/15  - NK: Revised ice factor initialization and calculation   
!     rev. 10.1.12 Dec.  12/15  - NK: Added Nash Efficiency nasheff.r2c file unit-66
!     REV. 10.1.13 Dec.  28/15  - NK: Rearranged the par file blocks & contents
!     REV. 10.1.14 Jan.  05/16  - NK: Added ice rules for Lakes Athabaska & Great Slave.
!     REV. 10.1.15 Jan.  08/16  - NK: Custom coding for Mackenzie River Basin Hydraulic Model
!     REV. 10.1.16 Jan.  11/16  - NK: Added subroutine ice_factor.f
!     REV. 10.1.17 Jan.  11/16  - NK: Added fpetLakeOverride factor
!     REV. 10.1.18 Jan.  15/16  - NK: Made opening of the master_inflow file optional with routeflg=q
!     REV. 10.1.19 Jan.  15/16  - NK: Fixed initialization of ice_factr - moved from lake_ice > runof6
!     REV. 10.1.20 Jan.  15/16  - NK: Fixed initialization of ice_factr - moved from lake_ice > runof6
!     REV. 10.1.21 Jan.  22/16  - NK: isotope updates
!     REV. 10.1.21 Jan.  23/16  - NK: Fixed lake init flow bug in flowinit
!     REV. 10.1.22 Jan.  25/16  - NK: Fixed flowinit for partial basins
!     REV. 10.1.23 Jan.  28/16  - NK: Added abort when water class not specified
!     REV. 10.1.24 Jan.  30/16  - NK: Added qUS1 & qUS2 for watbal
!     REV. 10.1.25 Feb.  21/16  - NK: Added nudge_flags.txt 
!     REV. 10.1.26 Mar.  23/16  - NK: Fixed comment for spinup period 
!     REV. 10.1.27 Apr.  19/16  - NK: Moved outfiles code in CHARM9 (below) 
!     REV. 10.1.28 Apr.  26/16  - NK: Fixed first day of output for master_inflows file 
!     REV. 10.1.29 May   04/16  - NK: Added parfile comments
!     REV. 10.1.30 May   08/16  - NK: Added smoothdist warning in read_par_parser
!     REV. 10.1.31 May   15/16  - NK: Revised output to precip.txt : include all str stations
!     REV. 10.1.32 May   18/16  - NK: Separate radinfl for precip & temperature
!     REV. 10.1.33 Jun   20/16  - NK: Change the time stamp in the watflood.wfo file
!     REV. 10.1.34 Jul   05/16  - NK: Added Obs. & Model mean flows to wfo file
!     REV. 10.1.35 Jul   07/16  - NK: Added simulation start time to the wfo file
!     REV. 10.1.36 Jul   12/16  - NK: Added results\LakeName.tb0
!     REV. 10.1.37 Jul   28/16  - NK: Added "Ellipsoid to the WFO header
!     REV. 10.1.38 Jul   28/16  - NK: Added noDataValue to WFO & tb0 files
!     REV. 10.1.39 Sep   16/16  - NK: Fixed stations outside the watershed for tb0
!     REV. 10.1.40 Oct   11/16  - NK: Fixed bug in read_divert for missing u/s DA
!     REV. 10.1.41 Oct   11/16  - NK: Added tb0flg to write lake_*.tb0 files
!     REV. 10.1.42 Oct   20/16  - NK: Reinstated read_ice_factor.f as default if present
!     REV. 10.1.43 Oct   21/16  - NK: lake_ice_facter changed from : to :,:  
!     rev. 10.1.44 Oct.  22/16  - NK: Reworked icerivflg & icelakeflg
!     rev. 10.1.45 Oct.  26/16  - NK: Added allocation check for qdivert in rerout
!     rev. 10.1.46 Nov.  08/16  - TH: Changed B1 - 5 to real*8
!     rev. 10.1.47 Nov.  08/16  - TH: Major changes in the ISO part of AET.f
!     rev. 10.1.48 Nov.  08/16  - TH: addet fpet(ii_water) to the wetland evaporation
!     rev. 10.1.49 Nov.  08/16  - TH: Overhauled lake evaporation
!     rev. 10.1.50 Nov.  08/16  - TH: Overhauled lst for new isotope output
!     rev. 10.1.51 Nov.  08/16  - TH: removed unused isotope related calculations, merged two 
!     rev. 10.1.51                TH: isotope related calc sections to reduce if statements
!     rev. 10.1.52 Nov.  08/16  - NK:
!     rev. 10.1.53 Nov.  09/16  - NK: Changed levels.txt to levels.csv
!     rev. 10.1.54 Nov.  25/16  - NK: Moved tdum under call timer in sub
!     rev. 10.1.55 Nov.  30/16  - NK: Fixed sumf & sumffs in runof6
!     rev. 10.1.56 Dec.  05/16  - NK: Fixed evt in AET.f to account for sca
!     rev. 10.1.57 Dec.  06/16  - NK: Added snwNN.txt files for iopt > 0
!     rev. 10.1.58 Dec.  06/16  - NK: corrected tdum > tdum1 for modelflg=i
!     rev. 10.1.59 Dec.  18/16  - NK: Fixed missing # channel correction chnl(1-5)
!     rev. 10.1.60 Jan.  03/17  - NK: Fixed conditional in route
!     rev. 10.1.61 Jan.  03/17  - NK: Changed results\peaks.txt to write peak flows
!     rev. 10.1.62 Jan.  08/17  - NK: Checkup on strloss effect on low flows
!     rev. 10.1.63 Jan.  25/17  - NK: Intel® Parallel Studio XE 2017 Update 1 
!     rev. 10.1.64 Jan.  26/17  - NK: Added XML output file 
!     rev. 10.1.65 Jan.  28/17  - NK: Fixed allocate lake_elv from read_flow 
!     rev. 10.1.66 Jan.  28/17  - NK: Fixed leap year in timer 
!     rev. 10.1.67 Feb.  18/17  - NK: Ignore start year in subsequent event files 
!     rev. 10.1.68 Mar.  03/17  - NK: Made midnight 00 instead of 24 
!     rev. 10.1.69 Mar.  03/17  - NK: Changed allocation for lvl_reach in read_lvl.f 
!     rev. 10.1.70 Mar.  03/17  - NK: Added year_now2 etc. for converting Grib2 files 
!     rev. 10.1.71 Mar.  14/17  - NK: Revised reading mean_observed_flows in sub 
!     rev. 10.1.72 Mar.  20/17  - NK: Fixed bug in sub for error_flag = 4 
!     rev. 10.1.73 Mar.  27/17  - NK: Advisory message set in precip.txt for iopt=0 
!     rev. 10.1.74 Apr.  01/17  - NK: Changed timer to fix 1 day-off problem 
!     rev. 10.1.75 Apr.  03/17  - NK: Fixed time & thr in runof6 arg list
!     rev. 10.1.76 Apr.  05/17  - NK: Reorganized the outfiles.* file
!     rev. 10.1.77 Apr.  17/17  - NK: Moved DDS err calcs to new dds_code s/r's
!     rev. 10.1.78 Apr.  17/17  - NK: New s/r dds_UZS to calculate low flow penalty
!     rev. 10.1.79 Apr.  18/17  - NK: Set trcflg=0 for all dds except errflg=10
!     rev. 10.1.80 Apr.  26/17  - NK: Fixed tracer turnoff for -ve resv. storage
!     rev. 10.1.81 May   05/17  - NK: Added snowg\yyyymmdd_swe.tb0 obs. swe
!     rev. 10.1.82 May   09/17  - NK: Added reservoir_fudge_factors.csv
!     rev. 10.1.83 May   09/17  - NK: Fixed lake evap bug - moved it outside lake-only loop
!     rev. 10.1.84 May   09/17  - NK: Put drng(n,ii)=drng(n,ii)*fraction back into runof6
!     rev. 10.1.85 May   17/17  - NK: Level_station_location.xyx for iopt > 0 only
!     rev. 10.1.86 May   17/17  - NK: Diversion_location.xyx for iopt > 0 only
!     rev. 10.1.87 May   18/17  - NK: Added DA to reservoir_location.xyz
!     rev. 10.1.88 May   23/17  - NK: Fixed Juliean_day problems for iso R/W
!     rev. 10.1.89 May   25/17  - NK: Added errflg = 11 for isotope DDS
!     rev. 10.1.90 Jul.  27/17  - NK: Added date_now for i/o files
!     rev. 10.1.91 May   25/17  - NK: Added errflg = 12 for isotope DDS
!     rev. 10.1.92 May   25/17  - NK: Changed to max 200 dds variables
!     rev. 10.1.92 Aug   12/17  - NK: delete store_dead in iso s/r's
!     rev. 10.1.93 Aug   17/17  - NK: allow year1 etc. to be passed for each event
!     rev. 10.1.94 Aug   29/17  - NK: Fixed col check bug in read_lvl
!     rev. 10.1.95 Sep   11/17  - NK: Fixed LKdepth bug in sub
!     rev. 10.1.96 Sep   11/17  - NK: Added variable lake depth calculation lake_elv()-LKinvert()
!     rev. 10.1.97 Sep   11/17  - NK: Moved hdrflg action in runof6.f
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
!     rev. 10.1.99 Oct   08/17  - NK: Added error check for # sdc classes in Melt.f
!     rev. 10.2.01 Oct   08/17  - NK: Moved ruleflg from sub.f to CHARM.f
!     rev. 10.2.02 Oct   24/17  - NK: Fixed xml output file
!     rev. 10.2.03 Oct   28/17  - NK: Revert to old G format for lakeSD.csv
!     rev. 10.2.04 Oct   28/17  - NK: Change to one xml output file for computed flow
!     rev. 10.2.05 Oct   28/17  - NK: Killed off stats_info.txt for iopt.ge.1
!     rev. 10.2.06 Oct   28/17  - NK: wfo_spec.txt in working OR basin directory
!     rev. 10.2.07 Nov.  03/17  - NK: New rt_pond subroutine for channel pond routing      
!     rev. 10.2.08 Nov.  04/17  - NK: New rt_channel & rt-_wetland subroutines      
!     rev. 10.2.09 Nov.  04/17  - NK: Reinstated old Manning's n correction for legacy files      
!     rev. 10.2.10 Nov.  04/17  - NK: Fixed XML file    
!     rev. 10.2.11 Dec.  18/17  - NK: 4 files added for BLEND.exe    
!     rev. 10.2.12 Dec.  30/17  - NK: Added frame headers to static r2c files incl. shd file   
!     rev. 10.2.13 Jan.  31/18  - NK: Re-wrote rules.f to mimic stop log operations    -> rules_sl.f
!     rev. 10.2.14 Jan.  31/18  - NK: Renamed rules.f to rules_tl.f - for use with target levels
!     rev. 10.2.15 Feb.  05/18  - NK: Added 'results\monthly_peaks'
!     rev. 10.2.16 Feb.  14/18  - NK: Added bankfull flow calculation
!     rev. 10.2.17 Feb.  26/18  - NK: Added WFruntime to read_flow_ef
!     rev. 10.2.18 Mar.  12/18  - NK: Fixed array fault in read_resv_ef and sub
!     rev. 10.2.19 Mar.  13/18  - NK: Fixed array fault read_divert.f
!     rev. 10.2.20 Apr.  06/18  - NK: Added res_next for NBS
!     rev. 10.2.21 Apr.  14/18  - NK: Added Lake Level update 
!     rev. 10.2.22 May   10/18  - NK: Set wetland classes uncoupled and coupled to same parameters
!     rev. 10.2.23 May   18/18  - NK: Revamped target level lake rules
!     rev. 10.2.24 May   21/18  - NK: Added error message in Read_rain & read_tmp
!     rev. 10.2.25 May   27/18  - NK: Fixed nash e calculation for value1=nopt=0
!     rev. 10.2.26 Jul.  03/18  - NK: Fixed nash e calculation for denominator = 0
!     rev. 10.2.27 Jul.  08/18  - NK: replaced lake_inflow_sum with temp_flow_sum
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
!     rev. 10.2.29 Aug.  21/18  - NK: Added warning for lake depths less than 1 m
!     rev. 10.2.30 Aug.  22/18  - NK: Added error for ftall(water) 
!     rev. 10.2.31 Aug.  23/18  - NK: Echo recorded levels for iopt>=99 only 
!     rev. 10.2.32 Aug.  23/18  - NK: Changed non-convergence value for qo2 
!     rev. 10.2.33 Sep.  14/18  - NK: Changed unit=42  fln(12) from clutter to model\*.r2c file 
!     rev. 10.2.34 Sep.  22/18  - NK: Initialize modelLastHour in rear_rain_ef
!     rev. 10.2.35 Oct.  08/18  - NK: Moved logical def. to area_watflood
!     rev. 10.2.36 Oct.  17/18  - NK: Initialized evap_rate(n) in lake_evap.f 
!     rev. 10.2.37 Oct.  18/18  - NK: Changed target levels to real*8 
!     rev. 10.2.38 Oct.  24/18  - NK: Added runReport.txt for forecast mode 
!     rev. 10.2.39 Nov.  15/18  - NK: changed snowc(n,ii) to snowc(n,ii)*sca(n,ii) in runof6
!     rev. 10.2.40 Nov.  22/18  - NK: GNU Lesser General Public License
!     rev. 10.2.41 Dec.  10/18  - NK: Added winter monthly peaks
!     rev. 10.2.42 Jan.  16/19  - NK: Added resume\yyyy-mm-dd folder with resume files
!     rev. 10.2.43 Jan.  17/19  - NK: Fixed bug in leapyear extra day in rules_echo.txt
!     rev. 10.2.44 Jan.  21/19  - NK: Added qOld() to allow previous weighted outflow on the resume file
!     rev. 10.2.45 Jan.  21/19  - NK: Fixed bug in reservoir initialization in sub - n was undefined
!     rev. 10.2.46 Jan.  21/19  - NK: Added netCDF capability
!     rev. 10.2.47 Feb.  10/19  - NK: Revised write_resume
!     rev. 10.2.48 Feb.  25/19  - NK: Moved temp correction to read_temp from process_temp
!     rev. 10.2.49 Jan.  21/19  - NK: Cleaned up write_2d_nc (netCDF)
!     rev. 10.2.50 Mar.  22/19  - NK: New resume files written to \resume\*.*
!     rev. 10.2.51 Apr.  03/19  - NK: New section to read ts5 format file for swe
!     rev. 10.2.52 Apr.  15/19  - NK: Added total UZS for reporting in FEWS
!     rev. 10.2.53 May.  09/19  - AJ: Added ensemble to read_2D for FEWS
!     rev. 10.2.54 June  07/19  - NK Fixed time stamp values for the wfo file when running FEWS
!     rev. 10.2.55 June  09/19  - NK Added read_swe_update.f90 for read in swe adjustment factors
!     rev. 10.2.56 June  13/19  - NK Added .nc output for grid_runoff & cumm ET
!     rev. 10.2.57 Jul.  07/19  - NK Tweaked coefficients in rules_tl
!     rev. 10.2.58 Jul.  17/19  - NK Remap nopt for changed flow stationlocations in FEWS
!     rev. 10.2.59 Aug.  17/19  - NK Convert read_flow_ef,read_lvl & read_divert to F90
!     rev. 10.2.60 Aug.  27/19  - NK Renamed read_lvl.f90 to read_level.f90 to get it added to the Repo
!     rev. 10.2.61 Aug.  28/19  - NK removed inbsnflg from read_resume
!     rev. 10.2.62 Sep.  09/19  - NK Added read_sm_update.f90 for read in sm adjustment factors
!     rev. 10.2.63 Sep.  09/19  - NK Fixed check for file exists for fln(54) - swe time series
!     rev. 10.2.64 Sep.  09/19  - NK added min_flow_cutoff for error calculations
!     rev. 10.2.65 Sep.  30/19  - NK Changed format of net_basin_supply outputs
!     rev. 10.2.66 Sep.  30/19  - NK Changed gage_temp from 7 to 8 char variable and used as flag ID
!     rev. 10.2.67 Nov.  03/19  - NK Fixed flow averaging in lst
!     rev. 10.2.68 Nov.  11/19  - NK New strfw\flowstation.xyz to set locations of the flow gages
!     rev. 10.2.69 Nov.  11/19  - NK New level\levelstation.xyz to set locations of the level gages
!     rev. 10.2.70 Nov.  11/19  - NK Lowered acceptable rvalue for Manning n correction 1.0 > 0.1
!     rev. 10.2.71 Nov.  18/19  - NK Bug fixes in wetland & reservoir routing
!     rev. 10.2.72 Nov.  21/19  - NK flow_sta_name: changed from 12 to 8 char  to accomodate US sta names
!     rev. 10.2.73 Dec.  14/19  - NK Convert rules_tl.f90 to Fortran 90 & fix drawdown comps.
!     rev. 10.3.00 Dec.  **/19  = NK Conversion f77 to f90      
!     rev. 10.3.01 Jan.  05/20  = Use nudge_flag.xyz file for station area check      
!     rev. 10.3.02 Jan.  17/20  = changed event\*.evt to event\*.cfg      
!     rev. 10.3.03 Jan.  21/20  = NK added dd_ice & dd_thaw to the resume.txt file      
!     rev. 10.3.04 Jan.  29/20  = NK added smrflg = "t" for FEWS only      
!     rev. 10.3.05 Feb.  06/20  = NK added climate data
!     rev. 10.3.06 Feb.  27/20  = NK Fixed temperature scaling factors for single events
!     rev. 10.3.07 Mar.  04/20  = NK Fixed weighted swe in wfo file for grids with water
!     rev. 10.3.08 Mar.  06/20  = NK For FEWS, add snow1\swe.nc for swe updating
!     rev. 10.3.09 Mar.  07/20  = NK Revise swe updating to maintain relative swe in classes
!     rev. 10.4.21 Apr.  21/20  = NK Add UZS deficit to wfo file = UZS(class=classcount)
!     rev. 10.4.22 Apr.  22/20  = NK fixed event precip scale factor  "readscale"
!     rev. 10.4.23 Apr.  30/20  = NK repurposed a2 for swe threshold"
!     rev. 10.4.24 Jun.  05/20  = NK ensemble revision
!     rev. 10.4.25 Jul.  25/20  = NK Changed rules_sl to use hourly rule intervales
!     rev. 10.4.26 Aug.  26/20  = NK Added Kling–Gupta efficiency (KGE) score
!     rev. 10.4.27 Sep.  28/20  = NK Write parfile.scv during dds run
!           
      program_name='CHARM     '
      program_version=' 10.4.27  '
      program_date='2020/10/05'

      
!     rev. 9.5.44  Oct.  27/08  - NK: removed code & obj modules for hasp & rainbow
!      NOT VALID IN UNIX
!      include 'watfile.fi'
 
!      DATA uzsinit/16*0.0/
!      DATA flgevp2/-1.0/
      
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      do i=1,100
          warning(i)=.false.
      end do
      
!     rev. 10.2.38 Oct.  24/18  - NK: Added runReport.txt for forecast mode 
      inquire(file='runReport.txt',EXIST=exists)
      if(exists)then
        open(unit=99,file='runReport.txt',status='unknown') 
        If(exists)then
          do while(.not.eof(99))
              read(99,*,iostat=ios)line
          end do
        endif
        write(99,*)'CHARM.exe - started'
        close(unit=99,status='keep')
      endif
      
      
!     stopflg is to keep the dos window open in WATFLOOD (VB)      
!     stopflg is to keep the dos window open in WATFLOOD (VB)      
!     If the program is run as:  CHARM 1  then the user has to hit return
!     when CHARM is finished.

      CALL GETARG(1, buf, status1)
      if(status1.ne.1)buf=' '
      if(buf.eq.'1')then
        stopflg='y'
      elseif(buf.eq.'9')then
!       rev. 9.1.63  Sep.  29/04  - NK: Added iopt_start as an arg for quick filecheck
      iopt_start=99
      elseif(buf.eq.'2')then
!       rev. 9.1.63  Sep.  29/04  - NK: Added iopt_start as an arg for quick filecheck
      iopt=2
      else
        stopflg='n'
      endif

!     increase stack size
!      editbin /stack:4000000 foo.exe

!     THIS SECTION GETS THE COMPUTER'S DATE AND TIME
!     THIS INFO USED TO COME FROM HEADER.FI, HOWEVER F90 INTRINSICS
!     AND CALLS ARE DIFFERENT AND THEREFORE IT NEEDED TO BE MODFIED

      call date_time(cday,time)

!     other unit numbers in use:
!     unit 29 - for domain_precip.txt in sub

!     INPUT FILES ARE UNITS 30-50  and  251- .....
!     THESE ARE OPENED MOSTLY IN SUB.F BECAUSE INPUT FILES
!     HAVE TO BE OPENED AND CLOSED FOR EACH MONTH IN THE ID LOOP

! unit=31  fln(1) - basin file (                    bsnm_shd.r2c)
! unit=32  fln(2) - parameter file (                bsnm.par)
! unit=33  fln(3) - point data location file (      bsnm.pdl)
! unit=34  fln(4) - not used in CHARM
! unit=35  fln(5) - point precipitation file        yyyymmdd_rag.tb0
! unit=36  fln(6) - streamflow data                 yyyymmdd_str.tb0  FEWS strfw\flow.nc
! unit=37  fln(7) - reservoir release data (        yyyymmdd_rel.tb0
! unit=38  fln(8) - reservoir inflow data           yyyymmdd_rin.tb0
! unit=39  fln(9) - unadjusted radar file           yyyymmdd_rad.r2c
! unit=39  fln(9) - reduced temperature grid        yyyymmdd_tem.r2c
! unit=40  fln(10)- precipitation data (            yyyymmdd_met.r2c  FEWS radcl\precip.nc
! unit=41  fln(11)- radar scan (cappi) file
!     rev. 10.2.33 Sep.  14/18  - NK: Changed unit=42  fln(12) from clutter to model\*.r2c file 
! unit=42  fln(12)- model r2c file                  yyyymmdd_capa.r2c
!                   note:  can be capa, regl_apcp or glb_apcp      
! unit=42  fln(12)- reduced precip grid             yyyymmdd_met.r2c
! unit=43  fln(13)- snow cover depletion curve      bsnm.sdc
! unit=44  fln(14)- point temperatures              yyyymmdd_tag.tb0
! unit=45  fln(15)- gridded temperatures            yyyymmdd_tem.r2c  FEWS tempr\tempt.nc
! unit=46  fln(16)- max temperatures                yyyymmdd.tmx
! unit=47  fln(17)- min temperatures                yyyymmdd.tmn
! unit=48  fln(18)- point daily snow files          yyyymmdd_dsn.tb0
! unit=49  fln(19)- radiation data gridded          yyyymmdd_flx.r2c
! unit=50  fln(20)- radiation data point            yyyymmdd_prn.tb0
! unit=251  fln(21)- humidity                       yyyymmdd_grh.r2c
! unit=252  fln(22)- wind speed                     yyyymmdd_gws.r2c
! unit=253  fln(23)- longwave radiation             yyyymmdd_glw.r2c
! unit=254  fln(24)- shortwave radiatin             yyyymmdd_gsw.r2c
! unit=255  fln(25)- atmospheric pressure           yyyymmdd_gpr.r2c
! unit=256  fln(26)- point relative humidity        yyyymmdd_prh.tb0       OR
! unit=256  fln(26)- point specific humidity        yyyymmdd_psh.tb0
! unit=257  fln(27)- point wind speed               yyyymmdd_pws.tb0
! unit=258  fln(28)- point longwave radiation       yyyymmdd_plw.tb0
! unit=259  fln(29)- point shortwave radiation      yyyymmdd_psw.tb0
! unit=260  fln(30)- point atmospheric pressure     yyyymmdd_ppr.tb0
! unit=261  fln(31)- gridded runoff files     runof\yyyymmdd_rff.r2c
! unit=262  fln(32)- gridded recharge         rchrg\yyyymmdd_rch.r2c
! unit=263  fln(33)- gridded leakage          lkage\yyyymmdd_lkg.r2c
! unit=264  fln(34)-  gridded lakeloss        eloss\yyyymmdd_lss.r2c
! unit=265  fln(35)- snow course data file          yyyymmdd_crs.pt2 
! unit=266  fln(36)- gridded snow water equivalant  yyyymmdd_swe.r2c
! unit=267  fln(37)- gridded soil moisture          yyyymmdd_gsm.r2c
! unit=268  fln(38)- gridded lower zone storage     yyyymmdd_lzs.r2c
! unit=269  fln(39)- point soil moisture            yyyymmdd_psm.pt2 
! unit=270  fln(40)- water quality data file  .wqd
! unit=271  fln(41)- routing parameter file         bsnm_ch_par.r2c
! unit=272  fln(42)- gridded wetland surface flux    - r2c file
! unit=273  fln(43)- gridded water surface flux file - 2rc file
! unit=274  fln(44)- point snow precip              yyyymmdd_snw.tbO
! unit=275  fln(45)- point drain (O18)              yyyymmdd_drn.tbO
! unit=276  fln(46)- point dsnow (O18)              yyyymmdd_dsn.tbO
! unit=277  fln(47)- gridded snow precip            yyyymmdd_snw.r2c
! unit=278  fln(48)- gridded drain (O18)            yyyymmdd_drn.r2c
! unit=279  fln(49)- gridded dsnow (O18)            yyyymmdd_dsn.r2c
! unit=280  fln(50)- point initial lake conditions  yyyymmdd_ill.pt2
! unit=281  fln(51)- gridded evaporation            yyyymmdd_evp.r2c
! unit=282  fln(52)- point diversion flow file      yyyymmdd_div.tb0
! unit=283  fln(53)- recorded lake level file       yyyymmdd_lvl.tb0  FEWS level\level.nc
! unit=284  fln(54)- swe time series pillows & crs  yyyymmdd_swe.ts5
! unit=285  fln(55)- init routing state variables   yyyymmdd_fli.r2c
! unit=286  fln(56)- point wind speed               yyyymmdd_spd.tb0
! unit=287  fln(57)- point wind direction           yyyymmdd_dir.tb0
! unit=288  fln(58)- gridded wind speed             yyyymmdd_spd.r2c
! unit=289  fln(59)- gridded wind direction         yyyymmdd_dir.r2c
! unit=290  fln(60)- mean daily temp diff    basin\mean_dly_diff.r2c
!           fln(61)- in usefor swe
! unit=292  fln(62)- daily temp differences         yyyymmdd_dif.r2c
!     rev. 10.2.11 Dec.  18/17  - NK: 4 files added for BLEND.exe    
! unit=293  fln(63)- Point Hourly Precip           yyyymmdd_pcp.tb0
! unit=294  fln(64)- Point Daily Precip            yyyymmdd_pcp.tb0
! unit=295  fln(65)- Gridded Hourly Precip         yyyymmdd_pcp.r2c
!           fln(66)- in used for results\nash_eff.r2c
!           fln(67)- in used for results\error.r2c
! unit=298  fln(68)- Gridded Daily precip          yyyymmdd_pcp.r2c

! unit=282  fln(52)- radar _1hr file    gzip1hr\yyyymm\yyyymmddhhmm_1hr
! unit=???  fln(201)-Gridded output for FEWS       CHARM_output.nc

! units 51-99 reserved for output files
! filenames 100-199 reserved for event file names

!     THE RESERVOIR INPUT FILE WAS ORIGINALLY 49 BUT THIS
!     WAS CHANGED WHEN TODD NEFF's EVAPORATION WAS ADDED

!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ALSO USED IN SHED: UNIT 9 FOR FLN(6)

    translateflg='n'  ! used on translate.for if 'y'
      ittoflg=0
      ssmc_firstpass='y'
      flgevp2=-1.0
      ensimopenflg='n'

!     TS - ALLOCATION FOR AREA12A ARRAYS
      allocate(fln(3000),filename(3000),outfln(3000),stat=iAllocate)
      if (iAllocate.ne.0)then
          write(*,*)'Error: error with allocation of fln arrays in CHARM'
          write(98,*)'Error: error with allocation of fln arrays in CHARM'
          STOP  
      endif

! SET DEFALUT FILENAMES FOR OUTPUT FILES:
!     these names may be replaced with the outfiles.txt file in the working
!     directory to send the files to a designated place for output.
!     A default outfiles.new file is created in the working directory each 
!     time this program is run.

!     unit number = file number

    filename(21)='results\monthly_peaks_winter.txt'          
    filename(22)='results\monthly_peaks.txt'          
      filename(23)='results\withdraw.r2c'
           fln(23)=filename(23)
    filename(25)='basin\new_par.csv'      
    filename(26)='results\parfile.csv'      
    filename(27)='results\precip.txt'             !dv vs classes
    filename(28)='results\stats.txt'              !statistics list
    filename(29)='results\domain_precip.txt'      !for BCHydro/CRA
    filename(30)='dds\dds_log.txt'                !for DDS pre-emption
!     units 31-50 reserved for input files
      filename(51)='debug\charm_info.txt'  !program information
      filename(52)='results\opt.txt'      !optimization dataF
      filename(53)='debug\res.txt'      !reservoir data
      filename(54)='debug\lake_error.txt'      
      filename(55)='debug\rte.txt'      !routing data
      filename(56)='56-notinuse'      !animation data for mapper
      filename(57)='results\snw.txt'      !snow pack info
      filename(58)='not_in_use'      
      filename(59)='results\stg.plt'      !output for stage plots
      filename(60)='results\spl.csv'      !paired observed/computed = default
!     for WATROUTE this file is named results\wrt.csv      
           fln(60)='basin\mean_dly_diff.r2c'
      filename(61)='results\swe.r2c'      !now in ensim format
           fln(61)=filename(61)             !used in write_r2c
      filename(62)='results\snw.csv'      !for snow plots
      filename(63)='debug\debug.txt'           
      filename(64)='debug\snwdebug.txt' !snow info
      filename(65)='results\watflood.wfo' !opened in wfocode
      filename(66)='results\nash_eff.r2c' !gridded nash efficiency
           fln(66)=filename(66)           !used in write_r2s
      filename(67)='results\error.r2c'    !gridded flow errors
           fln(67)=filename(67)           !used in write_r2s
      filename(68)='results\wetland.csv'  !wetland info
      filename(69)='results\sed.csv'      !sediment conc. obs/computed
      filename(70)='mrbhm\mrb_master_inflows.new'  !inflows for MRBHM
           fln(70)=filename(70)             !used in write_tb0
      filename(71)='results\CHARM_dly.csv'  !daily streamflow obs & comptd
      filename(72)='results\gridflow.r2c' !grid outflow fn(time) 
           fln(72)=filename(72)             !used in write_r2c
      filename(73)='results\resin.csv'    !lake inflow obs/computed
      filename(74)='results\evap.txt'     !evaporation debug data
      filename(75)='results\evt_means.csv'   ! used to be 'runoff.txt'
      filename(76)='results\peaks.txt'    !obs/computed Qp / event
      filename(77)='results\volumes.txt'  !obs/computed vol. / event
      filename(78)='results\CHARM_mly.csv'  !monthly streamflow          
      filename(79)='not_in_use'    ! for sum of all lzoutflow
      filename(80)='results\lake_sd.csv'  !  reservoir storage output
!     rev. 9.9.76  Sep.  11/15  - NK: Added recorded isotope concentrations


!     these 3 files are not included in the outfiles.txt file (for now)
      filename(81)='results\iso_info.txt'   !isotope messages, data echo, etc
      filename(82)='results\iso_18O_conc.csv'    ! 18O concentrations
      filename(83)='results\iso_2H_conc.csv'     ! 2H concentration
!     84-89 reserved for isotope output


      filename(90)='results\tracer.csv'     !tracer data
      filename(91)='results\tracerMB.csv'   ! added Oct.30/03 TS
      filename(92)='results\tracer_debug.csv'                       
      filename(93)='results\tracerWET.csv'    ! added Dec.01/03 TS
      filename(94)='results\tracerWETMB.csv'  ! added Dec.01/03 TS
      filename(95)='results\evapsep.txt'  ! added Apr.05/06 TS
      filename(96)='results\watbal1.csv'
      filename(97)='results\watbal2.csv'  ! runof6
      filename(98)='debug\warnings.txt'
      filename(99)='scratch5'             ! reserved as scratch file
      filename(100)='results\evap.r2c'    !weighted evaporation
           fln(100)=filename(100)             !used in write_r2c
      filename(201)='results\CHARM_flow_2D.nc'             !results\CHARM_output.nc
           fln(201)=filename(201)             !used in write_r2c
      filename(202)='results\CHARM_swe_2D.nc'             !results\CHARM_output.nc
           fln(202)=filename(202)             !used in write_r2c
      filename(203)='results\CHARM_uzs_2D.nc'             !results\CHARM_output.nc
           fln(203)=filename(203)             !used in write_r2c
      filename(204)='results\CHARM_grid_runoff_2D.nc'     !results\CHARM_output.nc
           fln(204)=filename(204)             !used in write_r2c
      filename(205)='results\CHARM_cumm_ET_2D.nc'         !results\CHARM_output.nc
           fln(205)=filename(205)             !used in write_r2c


!     reserve filename(2001-2999) for tb0 files with lake data: inflow/outflow/level

      
      
!     rev. 10.1.57 Dec.  06/16  - NK: Added snwNN.txt files for iopt > 0
      filename(801)='results\snw01.txt'     !snow plots class 1
      filename(802)='results\snw02.txt'     !snow plots class 2
      filename(803)='results\snw03.txt'     !snow plots class 3
      filename(804)='results\snw04.txt'     !snow plots class 4
      filename(805)='results\snw05.txt'     !snow plots class 5
      filename(806)='results\snw06.txt'     !snow plots class 6
      filename(807)='results\snw07.txt'     !snow plots class 7
      filename(808)='results\snw08.txt'     !snow plots class 8
      filename(809)='results\snw09.txt'     !snow plots class 9
      filename(810)='results\snw10.txt'     !snow plots class 10
      filename(811)='results\snw11.txt'     !snow plots class 11
      filename(812)='results\snw12.txt'     !snow plots class 12
      filename(813)='results\snw13.txt'     !snow plots class 13
      filename(814)='results\snw14.txt'     !snow plots class 14
      filename(815)='results\snw15.txt'     !snow plots class 15
      filename(816)='results\snw16.txt'     !snow plots class 16
      filename(817)='results\snw17.txt'     !snow plots class 17
      filename(818)='results\snw18.txt'     !snow plots class 18
      filename(819)='results\snw19.txt'     !snow plots class 19
      filename(820)='results\snw20.txt'     !snow plots class 20
      filename(821)='results\snw21.txt'     !snow plots class 21
      filename(822)='results\snw22.txt'     !snow plots class 22
      filename(823)='results\snw23.txt'     !snow plots class 23
      filename(824)='results\snw24.txt'     !snow plots class 24
      filename(825)='results\snw25.txt'     !snow plots class 25
      filename(826)='results\snw26.txt'     !snow plots class 26
      filename(827)='results\snw27.txt'     !snow plots class 27
      filename(828)='results\snw28.txt'     !snow plots class 28
      filename(829)='results\snw29.txt'     !snow plots class 29
      filename(830)='results\snw30.txt'     !snow plots class 30
      filename(831)='results\snw31.txt'     !snow plots class 30
      filename(832)='results\snw32.txt'     !snow plots class 30
      filename(833)='results\snw33.txt'     !snow plots class 30
      filename(834)='results\snw34.txt'     !snow plots class 30
      filename(835)='results\snw35.txt'     !snow plots class 30
      filename(836)='results\snw36.txt'     !snow plots class 30
      filename(837)='results\snw37.txt'     !snow plots class 30
      filename(838)='results\snw38.txt'     !snow plots class 30
      filename(839)='results\snw39.txt'     !snow plots class 30
      filename(840)='results\snw40.txt'     !snow plots class 30
      filename(841)='results\snw41.txt'     !snow plots class 30
      filename(842)='results\snw42.txt'     !snow plots class 30
      filename(843)='results\snw43.txt'     !snow plots class 30
      filename(844)='results\snw44.txt'     !snow plots class 30
      filename(845)='results\snw45.txt'     !snow plots class 30
      filename(846)='results\snw46.txt'     !snow plots class 30
      filename(847)='results\snw47.txt'     !snow plots class 30
      filename(848)='results\snw48.txt'     !snow plots class 30
      filename(849)='results\snw49.txt'     !snow plots class 30
      filename(850)='results\snw50.txt'     !snow plots class 30

!     rev. 9.7.13  Nov.  22/10  - NK: Changed the outfiles.txt for more 30 rff classes
      filename(901)='results\rff01.txt'     !runoff plots class 1
      filename(902)='results\rff02.txt'     !runoff plots class 2
      filename(903)='results\rff03.txt'     !runoff plots class 3
      filename(904)='results\rff04.txt'     !runoff plots class 4
      filename(905)='results\rff05.txt'     !runoff plots class 5
      filename(906)='results\rff06.txt'     !runoff plots class 6
      filename(907)='results\rff07.txt'     !runoff plots class 7
      filename(908)='results\rff08.txt'     !runoff plots class 8
      filename(909)='results\rff09.txt'     !runoff plots class 9
      filename(910)='results\rff10.txt'     !runoff plots class 10
      filename(911)='results\rff11.txt'     !runoff plots class 11
      filename(912)='results\rff12.txt'     !runoff plots class 12
      filename(913)='results\rff13.txt'     !runoff plots class 13
      filename(914)='results\rff14.txt'     !runoff plots class 14
      filename(915)='results\rff15.txt'     !runoff plots class 15
      filename(916)='results\rff16.txt'     !runoff plots class 16
      filename(917)='results\rff17.txt'     !runoff plots class 17
      filename(918)='results\rff18.txt'     !runoff plots class 18
      filename(919)='results\rff19.txt'     !runoff plots class 19
      filename(920)='results\rff20.txt'     !runoff plots class 20
      filename(921)='results\rff21.txt'     !runoff plots class 21
      filename(922)='results\rff22.txt'     !runoff plots class 22
      filename(923)='results\rff23.txt'     !runoff plots class 23
      filename(924)='results\rff24.txt'     !runoff plots class 24
      filename(925)='results\rff25.txt'     !runoff plots class 25
      filename(926)='results\rff26.txt'     !runoff plots class 26
      filename(927)='results\rff27.txt'     !runoff plots class 27
      filename(928)='results\rff28.txt'     !runoff plots class 28
      filename(929)='results\rff29.txt'     !runoff plots class 29
      filename(930)='results\rff30.txt'     !runoff plots class 30
      filename(931)='results\rff31.txt'     !runoff plots class 30
      filename(932)='results\rff32.txt'     !runoff plots class 30
      filename(933)='results\rff33.txt'     !runoff plots class 30
      filename(934)='results\rff34.txt'     !runoff plots class 30
      filename(935)='results\rff35.txt'     !runoff plots class 30
      filename(936)='results\rff36.txt'     !runoff plots class 30
      filename(937)='results\rff37.txt'     !runoff plots class 30
      filename(938)='results\rff38.txt'     !runoff plots class 30
      filename(939)='results\rff39.txt'     !runoff plots class 30
      filename(940)='results\rff40.txt'     !runoff plots class 30
      filename(941)='results\rff41.txt'     !runoff plots class 30
      filename(942)='results\rff42.txt'     !runoff plots class 30
      filename(943)='results\rff43.txt'     !runoff plots class 30
      filename(944)='results\rff44.txt'     !runoff plots class 30
      filename(945)='results\rff45.txt'     !runoff plots class 30
      filename(946)='results\rff46.txt'     !runoff plots class 30
      filename(947)='results\rff47.txt'     !runoff plots class 30
      filename(948)='results\rff48.txt'     !runoff plots class 30
      filename(949)='results\rff49.txt'     !runoff plots class 30
      filename(950)='results\rff50.txt'     !runoff plots class 30
!     rev. 9.8.92  Nov.  06/13  - NK: Changed output file swe.txt to swe.csv
      filename(951)='results\swe.csv'       !snow course comparison
      filename(953)='results\levels.csv'    !lake level comparison
      filename(954)='results\diversion.txt' !lake St. Jo. debug file
      filename(955)='results\CHARM.tb0'       !tb0 format output for CHARM
           fln(955)=filename(955)           !used in write_tb0
      filename(956)='results\lake_evap.txt' !lake evaporationdebug file
           fln(956)=filename(956)           !used in write_tb0
!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
      filename(957)='results\NetBasinSupply.txt'   !NBS
!     rev. 9.9.24  Aug.  20/14  - NK: Added monthly mean flow csv file CHARM_mly_nn.csv
      filename(960)='results\CHARM_mly_dlt.csv'      !daily comp-obs Q
      filename(961)='results\CHARM_mly_01.csv'       !monthly mean flow
      filename(962)='results\CHARM_mly_02.csv'       !monthly mean flow
      filename(963)='results\CHARM_mly_03.csv'       !monthly mean flow
      filename(964)='results\CHARM_mly_04.csv'       !monthly mean flow
      filename(965)='results\CHARM_mly_05.csv'       !monthly mean flow
      filename(966)='results\CHARM_mly_06.csv'       !monthly mean flow
      filename(967)='results\CHARM_mly_07.csv'       !monthly mean flow
      filename(968)='results\CHARM_mly_08.csv'       !monthly mean flow
      filename(969)='results\CHARM_mly_09.csv'       !monthly mean flow
      filename(970)='results\CHARM_mly_10.csv'       !monthly mean flow
      filename(971)='results\CHARM_mly_11.csv'       !monthly mean flow
      filename(972)='results\CHARM_mly_12.csv'       !monthly mean flow

      open(unit=51,file=filename(51),status='unknown',iostat=ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(ios.ne.0)call io_warning(51,filename(51),ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      open(unit=63,file=filename(63),status='unknown',iostat=ios)
        if(ios.ne.0)call io_warning(63,filename(63),ios)
      iendarg=0
!     id is used as a flag. when id=0 openfiles:
      id=1              !not used like this now  nk Apr. 8/03
      ni=1

! NOTE: UNIT 99 IS A SCRATCH FILE ONLY
!       CLOSE MUST FOLLOW OPEN IN SAME SEQUENCE
!           FILENAMES USING 99:
!             new.par
!             resume.par
!             tot.par
!             stop.txtc
!             outfiles.txt
!             losCHARMot.txt
!             newerror.txt
!             basin/evap.txt
!             event/event.evt
!             status='scratch'
!             fln(99)
!             basin/correct.tmp'
!             tempname.txt
!

! DURING EXECUTION, THE 'STOP' COMMAND CAN BE MADE FROM ANOTHER DOS
! WINDOW.  IT WILL CHECK AT THE END OF EACH EVENT.  THIS WILL CLOSE
! ALL FILES PROPERLY AS OPPOSED TO THE PROBLEMS WITH A CTRL/BRK CRASH
      open(unit=99,file='stop.txt',form='formatted',status='unknown',iostat=ios)
    if(ios.eq.0)then
        write(99,99001)
        close(unit=99,status='keep')
    else
      write(98,*)'Info: error opening stop.txt - new file not written'
    endif
      
! ioflg IS THE NUMBER OF OUTPUT FILES LISTED IN OUTFILES.TXT
! so it's value will be changed then

! OPEN CHARM_info.txt    MUST OPEN FIRST
! OPEN FILE FOR ALL CHARM ERROR MESSAGES:
      if(ioflg.gt.1)then
        filename(98)=outfln(98)
      endif
     
!     open warnings.txt     
      open(unit=98,file=filename(98),status='unknown',iostat=ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(ios.ne.0)call io_warning(98,filename(98),ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call date_and_time(cday,time)
      write(98,6016)'Info: ',time(1:2),time(3:4),time(5:6)
      write(98,6017)'Info: ',cday(1:4),cday(5:6),cday(7:8)

!     rev. 10.3.02 Jan.  17/20  = changed event\*.evt to event\*.cfg      
      fln(99)='event/event.cfg'
      INQUIRE(FILE=fln(99),EXIST=exists)
      if(exists)then
          fln(99)='event/event.cfg'
      else
          fln(99)='event/event.evt'
      endif
      write(98,*)'Warning: If program dies here, possibly there are forward '
      write(98,*)'Warning: slashes in the event file'
      write(98,*)'Warning: please convert / to \ in event files'
      write(*,*)'Warning: If program dies here, possibly there are forward '
      write(*,*)'Warning: slashes in the event file'
      write(*,*)'Warning: please convert / to \ in event files'
      
      iopt=1  
!     not read until par file is read which happens after
!     reading the event file to get the par file name
     
!**********************************************************************
!      call rdevt(date,conv,scale,smc5,nhr,nhf)
      call read_evt(date,conv,scale,smc5,nhr,nhf)
!**********************************************************************
      year_now=year1

       if(debug_output)write(63,*) ' In CHARM - before call rdpar'
  
! ORIGINAL FILEIO.FI SECTION - STARTS HERE: 

!   LEAVE FILEIO HERE BECAUSE rdevtA READS FILE NAMES THAT HAVE TO BE 
!                            OPENED

!   Rev. 7.78 modified for error checking - Sept.29/96  AC flight 
!   Rev  7.9  modified to open files for evaporation output
!   Rev. 8.3  - May.  22/97 -     added the simout/outfiles capability

       if(debug_output)write(63,*) ' In CHARM - 1190'

!     SHED READS IN ALL THE WATERSHED DATA, FROM SEPERATE PROGRAM:
 
       if(debug_output)write(63,*) ' In CHARM - before call shed'
      
!C//////////////////////////////////////////////
!C///////////////////////// 
!C// Added by Dave

!     rev. 9.9.35  Oct.  20/14  - NK: Added keyword & file checks
        if(IsFileTypeR2C(fln(1)))then
            call read_shed_ef(31,1)	
        else
!cc			call rdshed()
                 print*,'Old format shd files not accepted'
                 print*,'Please create EF ????_shd.r2c files & rerun'
                 write(98,*)'Error: Old format shd files not accepted-Program aborted in CHARM @ 994'
               stop 'Program aborted in CHARM9 @ 994'
        endif

!// End Dave addition
!C/////////////////////////
!C//////////////////////////////////////////////
       if(debug_output)write(63,*) ' In CHARM - before allocate'

!       TS - ALLOCATIONS OF AREA16A ARRAYS
        allocate(qrgrid(ycount+10,xcount+10),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: with allocation of area16a arrays in CHARM9'
            write(*,*)'Error: with allocation of area16a arrays in CHARM9'
            STOP  'Error with allocation of area16a arrays in CHARM9'
        endif
!       Initialize these values because many are outside grid
!       and would otherwise be undifined.  nk June 11/03
        do i=1,ycount+10
          do j=1,xcount+10
            qrgrid(i,j)=0.0
          end do
        end do

!       TS - ALLOCATIONS OF AREA6A ARRAYS
!       NK - nxtbasin is done in flowinit
        allocate(sn1(ycount,xcount),nhyd(ycount,xcount),&
       nbasin(ycount,xcount),basinerr(ycount,xcount),&
       nasheff(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*) 'Error: with allocation of area6a arrays in CHARM9'
            STOP 'Error with allocation of area6a arrays in CHARM9'
        endif
!       TS - ALLOCATION OF ARRAY RAD FROM AREAETA (REMAINDER)
        allocate(rad(ycount,xcount),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: with allocation of rad array in CHARM9'
            STOP 'Error with allocation of rad array in CHARM9'
        endif
!       TS - ALLOCATION OF AREAMELTA ARRAYS (REMAINDER)
!       ALLOCATION FOR SDCD,SDCSCA OCCUR IN RDSDCA.FOR
        allocate(snw(ycount,xcount),dsn(it,xcount-1),&
            tmx(it,xcount-1),tmn(it,xcount-1),el(it,xcount-1),&
            stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: with allocation of areamelta arrays in CHARM9 @ 1138'
            STOP 'Error with allocation of areamelta arrays in CHARM9 @ 1138'
        endif
        
    if(IsFileTypeR2C(fln(1))) then
!       TS - ALLOCATIONS OF AREAMELTA ARRAYS (PARTIAL)
!       SNW,DSN,TTEMP,TMX,TMN,EL ALLOCATED IN SHEDA.FOR
!       SDCD,SDCSCA ALLOCATED IN RDSDCA.FOR

!       done in rdshed for the old file types
!       now here so read_shed_ef can be the same for rte & CHARM
        allocate(snowc(na,classcount),dsnow(na),tmax(na),tmin(na),&
        tmin1(na),tmin2(na),sca(na,classcount),oldsca(na,classcount),&
        fexcess(na,classcount),snowcmin(na,classcount),&
        wcl(na,classcount),nsdc(classcount),snocap(classcount),&
        ati(na,classcount),def(na,classcount),&
        qtot(classcount),robg(classcount),rosn(classcount),qnet(classcount),&
        smelt(classcount),excess(classcount),extra(classcount),idump(classcount),&
        qrain(classcount),qsnow(classcount),qrn(na),qsn(na),glmelt(na),&
        qe(classcount),qh(classcount),qn(classcount),qp(classcount),&
        refrz(classcount),fmadj(na),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: Error with allocation of areamelta arrays in CHARM9 @1160'
        endif
!       TS - ALLOCATIONS OF AREAETA ARRAYS (PARTIAL)
!       RAD ALLOCATED IN SHEDA.FOR
!       TS - ADDED ALLOCATIONS FOR EVAP-SEPARATION PARAMS (22/03/06)
!       TS: CHANGED ALLOCATIONS OF alb,pet, evap TO classcount (27/03/06)
!       rev. 9.1.80  Mar.  31/05  - NK: added sublimation   (sublim)
!       parameter allocation moved to rdpar  27/07/06 nk

!       but oh my god. Not needed in rte so not in read_shed_ef
!       So moved here from rdshed. 

!       but strloss taken out!!!!!!

        allocate(vo(na,classcount),intev(na,classcount),&
        sublim(na,classcount),sum_sublim(na,classcount),&
        sum_et(na,classcount),sum_pet(na,classcount),pint(na,classcount),&
        ev(na,classcount),pet(na,classcount),fpet2(na),&
        x2(na,classcount),x3(na,classcount),intevt(na,classcount),&
        evt(na,classcount),ssumr(na,classcount),totint(na),eloss(na),&
        uzsinit(classcount),deficit(classcount),v1(na,classcount),&
        radv(na),sinlat(na),coslat(na),tanlat(na),&
        flgtemp(na),tto(na),ttomin(na),ttomax(na),&
        dd_ice(na),dd_thaw(na),stat=iAllocate)
 !      *rh(na),stat=iAllocate)
        if(iAllocate.ne.0)then
            write(98,*)'Error: Error with allocation of evt arrrays in CHARM9 @ 1185'
            STOP'Error with allocation of evt arrrays in CHARM9 @ 1185'
        endif

!     rev. 9.4.04  Apr.  23/07  - NK: moved allocate for melt from melt > CHARM
        if(snwflg.eq.'y')then
! TS - ADDED ARRAY ALLOCATIONS FOR MODULE MELTAONLY.FOR
            allocate(excess1(naa,classcount),raint(naa,classcount),&
            snowf(naa,classcount),snowt(naa,classcount),ta(naa,classcount),&
            deld(naa,classcount),dsno(naa,classcount),top(naa,classcount),&
            bot(naa,classcount),water(naa,classcount),wlmax(naa,classcount),&
            stat=iAllocate)
            if(iAllocate.ne.0)then
                write(98,*)'Error: Error allocating in CHARM @ 1195'
                STOP 'Error allocating in CHARM @ 1195'
            endif
        endif

        if(frcflg.eq.'y')then        
            allocate(evcg(na,classcount),acg(classcount),bcg(classcount),&
            relh(na),spech(na),delr(na),devcg(na,classcount),&
            icgflg(na,classcount),delsrf(na),delsm(na),&
            storeSW1(na,classcount),storeSW2(na,classcount),dels(na),&
            storeIF2(na,classcount),storeGW2(na),isodef(na,classcount),&
            isowcl(na,classcount),isowater(na,classcount),stat=iAllocate)
            if(iAllocate.ne.0)then
                write(98,*)'Error: Error with allocation of iso arrays in CHARM9 @ 1208'
                STOP 'Error with allocation of iso arrays in CHARM9 @ 1208'
            endif
        endif

!     rev. 9.5.03  Dec.  09/07  - NK: added reads for precip isotopes
        if(frcflg.eq.'y')then        
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
          allocate(dlt2H_snow(na),dlt2H_rain(na),dlt_rain(na),dlt_snow(na),stat=iAllocate)
          if(iAllocate.ne.0)then
            write(98,*)'Error: error with allocation of dlt_???? in CHARM9'
            STOP 'Error with allocation of dlt_???? in CHARM9'
          endif
        endif

!       AB Added to make sure that intevt(n,ii)=intevt(n,ii)+intev(n,ii) in intcept works May 10, 2002
        do i=1,ycount+10
          do j=1,xcount+10
            qrgrid(i,j)=0.0
          end do
        end do

      endif     ! if(IsFileTypeR2C(fln(1))) the
       if(debug_output)write(63,*) ' In CHARM - before call rdpar'

!     rev. 9.7.20  Jan.  31/11  - NK: Moved open statement for rdpar to rdpar/f
      call find_filetype(2)
      if(filetype.eq.'par'.or.filetype.eq.'PAR')then
      newparflg=.false.

!       **********************************************************************
        call rdpar(1,ix,e1)
!       **********************************************************************

!       rev. 9.5.05  Jan.  13/08  - NK: added check for rec() in CHARM
!       Check that rec is not too high in any grid & any class
!       Got rid of x(ii) in runof6
        recflg='n'
      recmax=1.0E+32
        do n=1,naa
        do ii=1,classcount-2
          if(sl1(n)*rec(ii).gt.1.0)then
            recflg='y'
            recmax=amin1(recmax,1.0/sl1(n))
            write(98,*)'Warning: Value of rec for class',ii,'is too high. ',sl1(n)*rec(ii)
          endif
        end do
      end do
      if(recflg.eq.'y')then
        write(98,*)'Warning: Max value for rec recommended for this watershed =',recmax
        write(98,*)'Warning: See charm.txt for infoon which grids'
        write(98,*)'Warning:You may try it at your own risk'
      endif

         if(debug_output)write(63,*) ' In CHARM: 3 - before call rdsdc'

!       **********************************************************************
        if(snwflg.eq.'y') call rdsdc()
!       **********************************************************************

!       temporary defaults
        radinfl=300
        smoothdist=35
        rainsnowtemp=0.0
!       replace the old per file format with a new format file
        if(iopt99)then
!       **********************************************************************
        call write_par_10(99,25)
!       **********************************************************************
        call write_par_10(99,26)
!       **********************************************************************
        call write_par_10(51,0)
!       **********************************************************************
        endif

        if(flgevp2.ne.4.0)then
!         not needed if daily differences are used        
          inquire(FILE='basin\monthly_climate_normals.txt',EXIST=exists)
          if(.not.exists)then
!           write a new monthly climate data file
!           previously part of the par file
          open(unit=99,file='basin\monthly_climate_normals.txt',status='unknown',iostat=ios)
          if(ios.ne.0)then
            write(98,*)'Error: Problem opening "basin\monthly_climate_normals.txt"'
            write(98,*)'Error: program aborted in CHARM @ 1509'
            print*,'program aborted in CHARM @ 1509'
          endif
          write(99,99011)
99011       format('month  jan  feb  mar  apr  may  jun  jul  aug  sep  oct  nov  dec')
            write(99,99012,iostat=ios)(diff(i),i=1,12)
99012       format('mxmn ',12f5.1)          
            write(99,99013,iostat=ios)(hu(i),i=1,12)
99013       format('humid',12f5.1)          
            write(99,99014,iostat=ios)(pres(i),i=1,12)
99014       format('pres ',12f5.1)          
          write(98,*)'Info: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          write(98,*)'Info: !the files:                          !'
          write(98,*)'Info: !basin\new_par.csv    &              !'
          write(98,*)'Info: !basin\monthly_climate_normals.txt   !'
          write(98,*)'Info: !have been written                   !'
          write(98,*)'Info: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          endif
        endif
   
      else

      newparflg=.true.

      debugflg=.true.  ! iopt from the par file will not be changed

!     rev. 9.8.05  Oct.  18/11  - NK: New read_par_parser subroutine
!       **********************************************************************
      call read_par_parser(32,2)
!	  call read_par(32,2)
!       **********************************************************************
!       write the par file in results\CHARM.txt
        if(iopt99)then
!         **********************************************************************
          call write_par_10(51,0)
!         **********************************************************************
        endif
!         write results\parfile.csv
!         **********************************************************************
          call write_par_10(99,26)
!         **********************************************************************

!       dds is a logical variable and can be used to suppress output
        if(dds_flag.eq.1)dds=.true.

!       open monthly peaks        
        if(iopt99)then
          open(unit=21,file=filename(21),status='unknown',iostat=ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(21,filename(21),ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          open(unit=22,file=filename(22),status='unknown',iostat=ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(22,filename(22),ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        endif
        
        

        if(flgevp2.ne.4.0)then
!         this file should be changed to an GK format file
!         and should get an entry in the event file
          inquire(FILE='basin\monthly_climate_normals.txt',EXIST=exists)
          if(exists)then
           write(98,*)'Info: found basin\monthly_climate_normals.txt file'
        else
             write(98,*)'Error: Could not find:"basin\monthly_climate_normals.txt"'
             write(98,*)'Error: This file not needed if daily differences were'
             write(98,*)'Error: calculated tempr\yyyymmdd_dif.r2c and flgevp = 4'
             write(98,*)'Error: in the par file to use the updated'
             write(98,*)'Error: Hargreaves ET solution'
             write(98,*)'Error: Program aborted in CHARM9 @ 1336'
           stop 'Program aborted in CHARM9 @ 1336'
        endif   
        
!         read the old monthly_climate_normals.txt file:      
        open(unit=99,file='basin\monthly_climate_normals.txt',status='unknown',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            write(98,*)'Error: Unable to open   basin\montly_climate_normals.txt'
            write(98,*)'Error: Possible cause(s):'
            write(98,*)'Error: file in use by another application'
            write(98,*)'Error: or target directory does not exist'
            write(98,*)'Error: Program aborted in CHARM.f @ 1587'
            write(*,*)'Error: Unable to open   basin\montly_climate_normals.txt'
            write(*,*)'Error: Possible cause(s):'
            write(*,*)'Error: file in use by another application'
            write(*,*)'Error: or target directory does not exist'
            write(*,*)'Error: Program aborted in CHARM.f @ 1587'
            stop 'Program aborted in CHARM.f @ 1587'
          endif
        write(51,*)  
        write(51,*)'Climate normals:'
        write(51,99011)
        read(99,*)line
          if(ios.ne.0)write(98,*)'Error: problems reading title line'
          if(ios.ne.0)write(*,*)'Error: problems reading title line'
          read(99,*,iostat=ios)junk,(diff(i),i=1,12)
          write(51,99012)(diff(i),i=1,12)
!          mean_dly_diff(n,ii) is not used - replaced by dly_diff(n)
          if(ios.ne.0)write(98,*)'Error: problems reading diff'
          if(ios.ne.0)write(*,*)'Error: problems reading diff'
          read(99,*,iostat=ios)junk,(hu(i),i=1,12)
          write(51,99013)(hu(i),i=1,12)
          if(ios.ne.0)write(98,*)'Error: problems reading humid'
          if(ios.ne.0)write(*,*)'Error: problems reading humid'
          read(99,*,iostat=ios)junk,(pres(i),i=1,12)
          write(51,99014)(pres(i),i=1,12)
          if(ios.ne.0)write(98,*)'Error: problems reading pres'
          if(ios.ne.0)write(98,*)'Program aborted in CHARM9 @ 1355'
          if(ios.ne.0)write(*,*)'Error: problems reading pres'
          if(ios.ne.0)write(*,*)'Program aborted in CHARM9 @ 1355'
        close(unit=99,status='keep')
        if(ios.ne.0)then
            
            stop 'Program aborted in CHARM9 @ 1355'
        endif
          write(51,*)
        
!         write the par file in results\CHARM.txt
!         **********************************************************************
        call write_par_10(51,0)
!         **********************************************************************
!         write results\parfile.csv
!         **********************************************************************
          call write_par_10(99,26)
!         **********************************************************************
        endif
      endif
!   THIS IS TO OPEN AND READ THE OPTIONAL OUTFILES.TXT FILE THAT SETS
!   THE LOCATIONS OF THE OUTPUT FILES.  IF FILE DOES NOT EXIST, DEFAULT
!   NAMES WILL BE USED.

      do i=1,999
          writeflg(i)=.true.
      end do
      
      if(iopt.eq.0.and.netCDFflg)then       
!           suppress output for netCDF production runs when iopt = 0          
            writeflg(1)=.false.      ! used as a generic flag that can be used anywhere
            writeflg(27)=.false.     ! precip.txt
            writeflg(51)=.false.  
            writeflg(59)=.false.     ! stage.plt
            writeflg(60)=.false.
            writeflg(63)=.false.
            writeflg(64)=.false.
            writeflg(65)=.false.     ! watflood.wfo
            writeflg(68)=.false.     ! used in rt_wetland
            writeflg(71)=.false.
            writeflg(78)=.false.
            writeflg(80)=.false.
            writeflg(953)=.false.    ! levels.csv
            writeflg(955)=.false.    ! CHARM.tb0
            ensimflg='n'
            nbsflg='n'
            write(98,*)'INFO: SUPRESSING OUTPUT for FEWS'
      else            
            write(98,*)'INFO: NOT SUPRESSING OUTPUT for FEWS'
      endif
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Moved all thhis from above so we have iopt99
      if(iopt99)then
          open(unit=52,file=filename(52),status='unknown',iostat=ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(52,filename(52),ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(51,5003)52,filename(52)
!         for dds we want to append - see errfg =9 in sub
        do while((.NOT.EOF(52)))
              read(52,*)line
          end do
      endif

      do i=53,55,2   ! don't open 54
      if(.not.i.eq.61)then
          if(ioflg.gt.1)then
            filename(i)=outfln(i)
          endif
          open(unit=i,file=filename(i),status='unknown',iostat=ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(i,filename(i),ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(51,5003)i,filename(i)
        endif
      end do
!     files 56 - 59 for wind and not checked here  
      
      do i=60,64
      if(.not.i.eq.61)then
          if(ioflg.gt.1)then
            filename(i)=outfln(i)
          endif

      
          if(writeFlg(i))open(unit=i,file=filename(i),status='unknown',iostat=ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(i,filename(i),ios)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(51,5003)i,filename(i)
        endif
      end do

       if(debug_output)write(63,*) 'passed 682 in CHARM'
      
!     File 61 (i=61) 'swe.r2c' is opened in write_r2c
!     File 65 (i=65) 'watflood.wfo' is opened in wfocode.for
!     File 66 (i=65) 'error.xyz' is opened in lst.for at end of run
!     File 67 (i=65) 'error.r2s' is opened in lst.for at end of run
      if(ioflg.gt.1)then
        do i=65,67
          filename(i)=outfln(i)
        end do
      endif

!     TS: CHANGED MAX LIMIT TO INCLUDE TRACER FILES + EVAPSEP FILE (APR 5/06)
      do i=68,80
        if(ioflg.gt.1)then
          filename(i)=outfln(i)
        endif
      if(i.ne.72.and.i.ne.70)then         ! 72 opened in write_r2c
!                                           ! 70 opened below if needed	  
          if(writeFlg(i))open(unit=i,file=filename(i),status='unknown',iostat=ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(i,filename(i),ios)
!          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(51,5003)i,filename(i)
      endif
!     REV. 10.1.18 Jan.  15/16  - NK: Made opening of the master_inflow file optional with routeflg=q
      if(i.eq.70.and.routeflg.eq.'q')then         ! open the MRB_master_inflow.tb0
          open(unit=i,file=filename(i),status='unknown',iostat=ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(i,filename(i),ios)
!          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(51,5003)i,filename(i)
          write(*,5003)i,filename(i)
      endif
      end do
      do i=90,95
        if(ioflg.gt.1)then
          filename(i)=outfln(i)
        endif
      if(i.ne.72)then         ! 72 opened in write_r2c
          open(unit=i,file=filename(i),status='unknown',iostat=ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(i,filename(i),ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(51,5003)i,filename(i)
      endif
      end do
      
!     This add the new run to the previous for CHARM.csv & lake_sd.csv
!     Added Feb. 23/19  NK      
      if(contflg.eq.'y')then
          do while(.not.eof(60))
              read(60,*)junk
          end do
          do while(.not.eof(80))
              read(80,*)junk
          end do
      endif
      
!     rev. 9.9.24  Aug.  20/14  - NK: Added monthly mean flow csv file CHARM_mly_nn.csv
      if(iopt99)then
          do i=960,972
              open(unit=i,file=filename(i),status='unknown',iostat=ios)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(ios.ne.0)call io_warning(i,filename(i),ios)
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              write(98,*)'Info: opened CHARM_mly files'
              write(51,5003)i,filename(i)
          end do
      endif

       if(debug_output)write(63,*) ' In CHARM -1001'

!     rev. 9.1.52  Mar.  11/04  - NK: continuous water quality modelling
!     wqual/yymmdd.wqd
      if(sedflg.eq.'y')then
        open(unit=256,file=fln(40),status='unknown',iostat=ios)
        if(ios.ne.0)then
          write(*,99119)fln(40)
          write(98,99119)fln(40)
99119     format(' Warning: Error opening or reading fln:',a30)
            write(98,*)'Error: Probable cause: read only yymmdd.wqd - wat qual file'
            write(98,*)'Error: This file is optional - used to enter nutrient loading'
            write(98,*)'Error: needed if sedflg=y in yymmdd.evt'
            write(98,*)'Error: OR: in config.sys have you set files=100 & buffers=50?'
            write(98,*)'Error: iostat code =',ios
            write(98,*)'Error: program aborted in CHARM.for @ 883'
            write(*,*)'Error: Probable cause: read only yymmdd.wqd - wat qual file'
            write(*,*)'Error: This file is optional - used to enter nutrient loading'
            write(*,*)'Error: needed if sedflg=y in yymmdd.evt'
            write(*,*)'Error: OR: in config.sys have you set files=100 & buffers=50?'
            write(*,*)'Error: iostat code =',ios
            write(*,*)'Error: program aborted in CHARM.for @ 883'
          stop 'program aborted in CHARM.for @ 883'
        endif
      endif

      write(51,6015)time(1:2),time(3:4),time(5:6),cday(1:4),cday(5:6),cday(7:8)

!      write(51,6015) hrs,mins,secs,day,month,year

99905 write(51,1000)fln(10),fln(1)

       if(debug_output)write(63,*) ' In CHARM - 1183'

!     don't know what this is used for <<<<<<nl  18/04/02
!      type=float(itype)

      if(iopt.ge.1)then
        write(51,6001) iopt,itype
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     REV. 10.1.27 Apr.  19/16  - NK: Moved outfiles code in CHARM9 (below) 
!     rev. 10.1.76 Apr.  05/17  - NK: Reorganized the outfiles.* file
! WRITE A NEW OUTFILES.TXT FILE THAT CAN BE MODIFIED BY THE USER:
      if(writeFlg(1))then
        open(unit=99,file='outfiles.new',form='formatted',status='unknown',iostat=ios)
        if(ios.eq.0)then
          write(99,*)'unit#  Filename'
          write(99,99002)(i,filename(i),i=26,30)
          write(99,99002)(i,filename(i),i=51,80)
          write(99,99002)(i,filename(i),i=90,100)
          write(99,99002)(i,filename(i),i=901,900+classcount)
          close(unit=99,status='keep')
        else
          write(98,*)'Error: Unable to open file  outfiles.new'
          write(98,*)'Error: Possible cause(s):'
          write(98,*)'Error: file in use by another application'
          write(98,*)'Error: Program aborted in CHARM.f @ 865'
          write(*,*)'Error: Unable to open file  outfiles.new'
          write(*,*)'Error: Possible cause(s):'
          write(*,*)'Error: file in use by another application'
          write(*,*)'Error: Program aborted in CHARM.f @ 865'
          stop 'Program aborted in CHARM.f @ 865'
        endif
      endif
      
!     rev. 10.1.76 Apr.  05/17  - NK: Reorganized the outfiles.* file
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call read_outfiles
!      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! THIS INCLUDE OPENS THE SIMOUT FILES AND WILL BE DIFFERENT FOR UNIX 
! OR DOS\WINDOWS
     
      call date_and_time(cday,time)
      if(numa.eq.0.and.dds_flag.eq.0)then
        write(51,6011)program_version,program_date
        write(51,6016)time(1:2),time(3:4),time(5:6)
        write(51,6017)cday(1:4),cday(5:6),cday(7:8)
        write(51,*)
        write(51,5003)i,filename(i)

        write(51,5002)
        write(51,1030)(i,i,filename(i),i=51,80)
        write(51,1030)(i,i,filename(i),i=901,930)
        write(51,1030)(i,i,filename(i),i=90,100)
    endif

      allocate (dly_diff(na),stat=iAllocate)
 
!     rev. 9.8.77  Jul   08/13  - NK: Made universal the use of wetland_flag(n)
!     added mar 28/06  nk
!     rev. 9.2.35  Mar.  22/06  - NK: Glacier flow bypasses wetlands
      if(glacier_class_number.ne.0)then
!       glacier_class_number is assigned in the rdpar file
        do n=1,naa
          if(aclass(n,glacier_class_number).gt.0.0)then
!         there is a glacier in this grid
          glacier_flag(n)='y'
          else
!         there is no glacier in this grid
          glacier_flag(n)='n'
          endif
        end do
      else  
!     rev. 9.9.18  Jun.  08/14  - NK: Fixed glacier_class check for wetlands
!       there are no glaciers at all      
        do n=1,naa
          glacier_flag(n)='n'
        end do
      endif

!     rev. 9.8.29  Oct.  15/12  - NK: added wetland_flag to speed up route.f
!     wetland_flag(n) replaced multiple logical checks by just 1 in each
!     time & grid loop - so will save a lot of time
      allocate(wetland_flag(naa),stat=iAllocate)
      if(iAllocate.ne.0)then
          write(98,*)'Error:Error with allocation of wetland_flag in CHARM9 @ 1496'
          STOP 'Error with allocation of wetland_flag in CHARM9 @ 1496'
      endif

      write(51,*)'Wetland flag settings:'
      do n=1,naa
        if(wetflg.eq.'y'&
                .and.grid_area(n).gt.10.0&   
                .and.aclass(n,classcount-2).gt.0.0&         
                .and.theta(n).gt.0.00001&        
                .and.glacier_flag(n).ne.'y')then     
          wetland_flag(n)=.true.
        else
          wetland_flag(n)=.false.
        endif
        write(51,*)n,grid_area(n),aclass(n,classcount-2),ibn(n),theta(n),glacier_flag(n),wetland_flag(n)
    end do  
       
!     rev. 9.8.02  Jul.  26/11  - NK: reactivated meander length
!     this can only be placed here after reqding both the shd & par files
!     as with the PS it could be recalculated with each iteration      
      do n=1,naa
        rl(n)=rl(n)*mndr(n)
      end do
      
!     rev. 10.1.57 Dec.  06/16  - NK: Added snwNN.txt files for iopt > 0
!     open the snw files      
      if(iopt.ge.2)then
        do i=801,800+classcount
          open(unit=i,file=filename(i),status='unknown',iostat=ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if(ios.ne.0)call io_warning(i,filename(i),ios)
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          write(51,5003)i,filename(i)
        end do
      endif

!     open the rff files      
      if(iopt.ge.1)then
          write(98,*)'Info: results\rff**.txt files - if any:'
          write(98,*)'Info:      class #         file name              class fraction'
        do ii=1,classcount
          if(ioflg.gt.1)then
            filename(ii)=outfln(ii)
          endif
!         Alwasy open a rff file for the water class as some area is given 
!         in flowinit          
          if(aclass(nnprint,ii).gt.0.00)then
              write(98,*)'Info: ',ii,filename(900+ii)(1:30),aclass(nnprint,ii)
            open(unit=900+ii,file=filename(900+ii),status='unknown',iostat=ios)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if(ios.ne.0)call io_warning(ii,filename(900+ii),ios)
!           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          endif
        end do
        
!     rev. 9.8.31  Oct.  16/12  - NK: continue rff files for contflg = y
        if(contflg.eq.'y')then
!         want rff files to continue from the spinup run
!         so read to end of files        
          do ii=1,classcount
            if(aclass(nnprint,ii).gt.0.00)then
              do while(.not.eof(900+ii))
                  read(900+ii,*)
              end do      
            endif
          end do      
        endif  
      else
!       delete existing files      
        if(.not.contflg.eq.'y')then
        do ii=1,classcount
          INQUIRE(FILE=filename(900+ii),EXIST=exists)
          if(exists)then
            open(unit=99,file=filename(900+ii),status='unknown',iostat=ios)
            close(unit=99,status='delete')
            write(98,*)'Info: rff??.txt file ',filename(900+ii)(1:25),' deleted'
          endif
        end do  
        endif
      endif

!     rev. 10.2.01 Oct   08/17  - NK: Moved ruleflg from sub.f to CHARM.f
!     rev. 9.9.65  Apr.  `3/15  - NK: Added rule s/r; resrl\rules.txt & ruleflg
      INQUIRE(FILE='resrl\rules.ts5',EXIST=exists)
      if(exists)then
!         this means there are rules for some or all of the lakes & reservoirs
!         Note:  not necessarily all        
          ruleflg=.true.
          ruletype='aaaaaaaaaaaa'
          open(unit=99,file='resrl\rules.ts5',status='old')
          foundEndHeader=.false.
          do while(.not.foundEndHeader)
            read(99,99000)line
99000       format(a4096)       
            if(line(1:10).eq.':endHeader')foundEndHeader=.true.   
            if(line(1:10).eq.':EndHeader')foundEndHeader=.true.   
            if(line(1:8).eq.':StopLog')     ruletype='StopLog     '  
            if(line(1:12).eq.':TargetLevel')ruletype='TargetLevel '  
          end do
          write(98,*)'Info: ruletype = ',ruletype
      else
!         this means there are no rules at all.        
          ruleflg=.false.
      endif      


       if(debug_output)write(63,*) ' In CHARM: 3 - before call options'

!     **********************************************************************
    call options(ix,e1,conv,scale,smc5,nhr,nhf)
!     **********************************************************************

       if(debug_output)write(63,*)'in CHARM9 @1'



!     rev. 10.2.28 Jul.  08/18  - NK: Revised for TSW = NBS without lake precip/evap'n
      if(nbsflg.eq.'y')then    ! for the NBS file
    write(957,*)
      write(957,*)'DISCLAIMER'
    write(957,*)'The WATFLOOD software and other material supplied' 
    write(957,*)'in connection herewith is furnished by N. Kouwen and the' 
    write(957,*)'University of Waterloo and is accepted by the' 
    write(957,*)'user upon the express understanding that N. Kouwen' 
    write(957,*)'or the University of Waterloo make no warranties, either' 
    write(957,*)'express or implied, ' 
    write(957,*)'concerning the accuracy, completeness,' 
    write(957,*)'reliability, usability, performance, or fitness for any' 
    write(957,*)'particular purpose.' 
    write(957,*)
    write(957,*)'The material is provided "as is". The entire risk as to' 
    write(957,*)'its quality and performance is with the user.'
      write(957,*)
    write(957,*)'The forecasts produced by the WATFLOOD software are for' 
    write(957,*)'information and discussion purposes only and are not to' 
    write(957,*)'be relied upon in any particular situation without the' 
    write(957,*)'express written consent of N. Kouwen or the' 
    write(957,*)'University of Waterloo.'
      write(957,*)
      endif
              
      

      do i=901,930
        close(unit=i,status='keep',iostat=ios)
        if(ios.ne.0)then
          write(98,*)'Warning:  Problems on unit ',i,filename(i)
          write(98,*)'Warning:  error in closing file name',filename(i)
          write(98,*)'Warning:  iostat code =',ios
        endif
      end do

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      CALL ERRSNS(io_err,SYS_ERR,STAT,UNIT,cond)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      call errormsg(io_err,SYS_ERR,STAT,UNIT,cond)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      beepflg='on'
      INQUIRE(FILE='beep.txt',EXIST=exists)
      IF(exists)THEN
      open(unit=99,file='beep.txt',status='unknown')
      read(99,99005)beepflg
    endif

      if(beepflg.ne.'off')then
!       take out for unix - non standard
!       CALL BEEPQQ (frequency, duration)
!       CALL SLEEPQQ(delay)
        call BEEPQQ(2000,100)
        call SLEEPQQ(100) 
        call BEEPQQ(2000,100)
!       take out for unix - non standard
      endif

      if(fratioflg.and.dds_flag.eq.0.and.iopt99)then
        write(98,*)'fratio(s) <> 1.0 so maybe you would like to correct'
        write(98,*)'the interception capacity values h(month,class)'
        write(98,*)'in the par if you are happy with the values'
        write(98,*)'Find the corrected h(,) values and fratio=1.0'
        write(98,*)'in results\parfile.csv  Copy this as you par file'
      endif
      
!     rev. 9.9.39  Nov.  14/14  - NK: Modifications for watroute
      if(.not.dlyflg.and.flgevp2.eq.4.and.modelflg.eq.'n')then
        write(98,*)'Warning: Samani & Hargreaves (1985) eqn. used without use'
        write(98,*)'Warning: of daily temeprature differences as intended'
        write(98,*)'Warning: If this is not intended, please run TMP to create'
        write(98,*)'Warning: the yyyymmdd_dif.r2c files.'
      endif  

      if(.not.courantflg)then
        write(98,*)'Warning: Courant criterion violated'
        write(98,*)'Warning: Likely problems with instability & FP overflows'
        write(98,*)'Warning: Please reduce the min time step A6 in the par file'
      endif

      if(iopt.ge.1)close (unit=51,status='keep')
      
!     rev. 10.1.01 Oct.  05/15  - NK: Isotope update: added 2H
      !TH: stuffing the output framework stuff here for 2H isotopes, since I only want output at the end, move if you know a better spot
      if(frcflg.eq.'y'.and.flg2H.eq.2)then
      
       do n=1,nisoframe
        do i=1,int(totaltime)
         sumRO(n)=sumRO(n)+bsnRO(i,n)*bsnR(i,n)
         sumR2H(n)=sumR2H(n)+bsnR2H(i,n)*bsnR(i,n)
         sumRP(n)=sumRP(n)+bsnR(i,n)
        enddo
        do i=1,int(totaltime)
         sumRNum(n)=sumRNum(n)+(bsnRO(i,n)-sumRO(n)/sumRP(n))*(bsnR2H(i,n)-sumR2H(n)/sumRP(n))*bsnR(i,n)
         sumRDen(n)=sumRDen(n)+(bsnRO(i,n)-sumRO(n)/sumRP(n))**2*bsnR(i,n)
        enddo
       enddo
       
       !mX=sumRO(n)/sumRP(n)
       !mY=sumR2H(n)/sumRP(n)
       !Bx=dA-5= (1000*dAtotal(n)/Etotal(n)-5)
       !Ex=d*+1= (dStotal(n)/Etotal(n)+2)
            
        open(2228,file='results\frame.csv',status='unknown',iostat=ios)
          if(ios.ne.0)write(98,*)'Error: Error opening frame.csv output file'
          if(ios.ne.0) STOP 'Error opening frame output file'
          write(2228,1997)
          write(2228,1998)((1000*dAtotal(n)/Etotal(n)-5),&
                8*(1000*dAtotal(n)/Etotal(n)-5)+10,&
                (1000*dAtotal(n)/Etotal(n)-5),&
                sumRNum(n)/sumRDen(n)*(1000*dAtotal(n)/Etotal(n)-5)&
                +sumR2H(n)/sumRP(n)-sumRNum(n)/sumRDen(n)*sumRO(n)/sumRP(n),&
                (1000*dAtotal(n)/Etotal(n)-5),&
                ((d2HEtotal(n)-d2HStotal(n))/(dEtotal(n)-dStotal(n)))&
                *(1000*dAtotal(n)/Etotal(n)-5)+d2HStotal(n)/Etotal(n)&
                -(((d2HEtotal(n)-d2HStotal(n))/(dEtotal(n)-dStotal(n)))&
                *dStotal(n)/Etotal(n)),&
                1000*dAtotal(n)/Etotal(n),1000*d2HAtotal(n)/Etotal(n),&
                dStotal(n)/Etotal(n),d2HStotal(n)/Etotal(n),&
                  n=1,nisoframe)
          write(2228,1999)((dStotal(n)/Etotal(n)+1),&
                8*(dStotal(n)/Etotal(n)+1)+10,&
                (dStotal(n)/Etotal(n)+1),&
                sumRNum(n)/sumRDen(n)*(dStotal(n)/Etotal(n)+1)+sumR2H(n)/sumRP(n)&
                -sumRNum(n)/sumRDen(n)*sumRO(n)/sumRP(n),&
                (dStotal(n)/Etotal(n)+1),((d2HEtotal(n)-d2HStotal(n))/&
                (dEtotal(n)-dStotal(n)))*(dStotal(n)/Etotal(n)+1)&
                +d2HStotal(n)/Etotal(n)-(((d2HEtotal(n)-d2HStotal(n))/&
                (dEtotal(n)-dStotal(n)))*dStotal(n))/Etotal(n),&
                n=1,nisoframe)
 1997 FORMAT(<nisoframe>('GMWL',',',',','LMWL',',',',','LEL',',',',','del a',',',',','del *',',',','))
 1998 FORMAT(<nisoframe>(10(f20.6,',')))
 1999 FORMAT(<nisoframe>(6(f20.6,','),',',',',',',','))
      end if



    print*,'DISCLAIMER'
    print*,'The WATFLOOD software and other material supplied' 
    print*,'in connection herewith is furnished by N. Kouwen and the' 
    print*,'University of Waterloo and is accepted by the' 
    print*,'user upon the express understanding that N. Kouwen' 
    print*,'or the University of Waterloo make no warranties, either' 
    print*,'express or implied, concerning the accuracy, completeness,' 
    print*,'reliability, usability, performance, or fitness for any' 
    print*,'particular purpose.' 
    print*
    print*,'The material is provided "as is". The entire risk as to' 
    print*,'its quality and performance is with the user.'
    print*
    print*,'The forecasts produced by the WATFLOOD software are for' 
    print*,'information and discussion purposes only and are not to' 
    print*,'be relied upon in any particular situation without the' 
    print*,'express written consent of N. Kouwen or the' 
    print*,'University of Waterloo.'
    print*

      
      if(error_msg.eq.1)then
        write(98,*)'Warning:'
        write(98,*)'Warning: # of events to follow is > 0 for 1 or more '
        write(98,*)'Warning: chained events other than the first event file.'  
        write(98,*)'Warning: Only the first event file is used to set the'
        write(98,*)'Warning: sequence of events - othere are ignored'
      endif
      
!     rev. 10.1.98 Oct   04/17  - NK: Deal with -ve flows in route
      if(warning(1).and.iopt99)then
          write(98,*)'Warning:'
          write(98,*)'Warning: -ve flows encountered in some grids'
          write(98,*)'Warning: Please see warnings.txt for grid #(s) and time'
          write(98,*)'Warning: Lowest value encountered = ',qlow
          write(98,*)'Warning: id ',lowid,' time ',lowtime,' grid# ',lowgrid
          write(98,*)'Warning: row ',yyy(lowgrid),' col ',xxx(lowgrid)
          write(98,*)'Warning: Probably cause: evaporation from streams & rivers'
          write(98,*)'Warning: too large after prolonged dry periods'
          write(98,*)'Warning: Please see qlow.xyz for locations with -ve flows'
      endif
      if(warning(2))then
            write(98,*)'Warning: Very low flows encountered in some grids'
            write(98,*)'Warning: and evaporation from streams set = 0 This may not be a problem. '
            write(98,*)'Warning: Set iopt=1 to see grid #(s) and time'
      endif

      write(98,*)
      write(98,*)'DISCLAIMER'
    write(98,*)'Info: The WATFLOOD software and other material supplied' 
    write(98,*)'Info: in connection herewith is furnished by N. Kouwen and the' 
    write(98,*)'Info: University of Waterloo and is accepted by the' 
    write(98,*)'Info: user upon the express understanding that N. Kouwen' 
    write(98,*)'Info: or the University of Waterloo make no warranties, either' 
    write(98,*)'Info: express or implied, ' 
    write(98,*)'Info: concerning the accuracy, completeness,' 
    write(98,*)'Info: reliability, usability, performance, or fitness for any' 
    write(98,*)'Info: particular purpose.' 
    write(98,*)
    write(98,*)'Info: The material is provided "as is". The entire risk as to' 
    write(98,*)'Info: its quality and performance is with the user.'
      write(98,*)
    write(98,*)'Info: The forecasts produced by the WATFLOOD software are for' 
    write(98,*)'Info: information and discussion purposes only and are not to' 
    write(98,*)'Info: be relied upon in any particular situation without the' 
    write(98,*)'Info: express written consent of N. Kouwen or the' 
    write(98,*)'Info: University of Waterloo.'
      write(98,*)
          
!     rev. 10.2.38 Oct.  24/18  - NK: Added runReport.txt for forecast mode 
      if(.not.netCDFflg)then
        inquire(file='runReport.txt',EXIST=exists)
        if(exists)Then
         open(unit=99,file='runReport.txt',status='unknown') 
           If(exists)then
             do while(.not.eof(99))
               read(99,*,iostat=ios)line
             end do
           endif
           write(99,*)'CHARM.exe - normal ending'
           close(unit=99,status='keep')
         endif
      endif
      
!     save the CHARM.csv file for forecast statistical analysis 
      if(ni.gt.1)then
!       event_fln not allocated for one event run          
        if(event_fln(ni)(1:13).eq.'event\glb.evt')then
          call date_and_time(cday,time)
          write(line,99006)cday(1:4),cday(5:6),cday(7:8)
99006     format('copy results\CHARM.csv results\',a4,a2,a2,'_CHARM.csv') 
          print*
          print*,line(1:45)
          CALL execute_command_line(line(1:45))
          print*
        endif
      endif
      
!     stopflg is to keep the dos window open in WATFLOOD (VB)      
      if(stopflg.eq.'y')then
        print*,' normal ending'
        PAUSE ' Hit enter to exit window'
        print*
        stop 
      else
        print*,'Please see debug\warnings.txt for messages'
        print*
        stop ' normal ending'
      endif



! FORMATS

 1000 format(' ',2(' ',a30))
 1001 format(39x,32a1)
 1002 format(' ','key = ',i10)
 1003 format(' ',33a1)
 1030 format(' ','Unit no. =',i3,' file no',i3,' = ',a30)
 1050 format(26x,a30)
 3000 format(a5,7i5,25x,f10.0)
 3001 format(a20,i12)
 3010 format(a1,a79)
 5001 format(a30)
 5002 format(/' output files')
 5003 format(' opened unit',i5,' file name ',a30)
 5004 format(i10,a30)
 5006 format(/' Output file names ')
 5007 format(' runtime & date ',2a20)
 5008 format(' closed unit',i5,' file name ',a30)
 6001 format(' CHARM7: iopt=',i5,' itype=',i5)
 6040 format(' .met file not found. run radmet or ragmet again')
 6010 format(12i5)
 6011 format(1x,'*      ver=',2a10,'          *',2x,a30)
 6015 format(' runtime  ',a2,':',a2,':',a2,2x,a4,'-',a2,'-',a2)
! 6015 format(' runtime ',2(i2,':'),i2,2x,2(i2,'/'),i4)
 6016 format(a10,'  runtime    ',2(a2,':'),a2)
 6017 format(a10,'  rundate  ',a4,'-',a2,'-',a2)
 9021 format(f10.0,f5.0,3i5)
 9804 format(a5,f10.0,a60)
 9805 format(a5,i10,a60)
99001 format('  0.0')
99002 format(i5,1x,a256)
99003 format(i1,5x,a50)
99004 format(i5)
99005 format(a3)

      END PROGRAM CHARM     
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine io_warning(unit_number,file_name,ios)


      use area_watflood
      implicit none

      integer  :: unit_number,ios
      character*256  :: file_name

      write(98,*)'Error: Problems on unit',unit_number
      write(98,*)'Error:  error in opening file name  '
      write(98,*)'Error: ',file_name(1:60)
      write(98,*)'Error: possible cause: existing file is read-only'
      write(98,*)'Error:     which happens if it is open by say excel'
      write(98,*)'Error:     or some other program'
      write(98,*)'Error: or folder does not exist    <<<'
      write(98,*)'Error: or file does not exist      <<<'
      write(98,*)'Error: or DISK does not exist      <<<'
      write(98,*)'Error: iostat code =',ios
      write(98,*)'Error: program aborted in io-warning.for (CHARM)  @1054'
      
      print*,'Problems on unit',unit_number
      print*,'Warning:'
      print*,' error in opening file name  '
      print*,file_name(1:60)
      print*,'possible cause: existing file is read-only'
      print*,'    which happens if it is open by say excel'
      print*,'    or some other program'
      print*,'or folder does not exist    <<<'
      print*,'or file does not exist      <<<'
      print*,'or DISK does not exist      <<<'
      Print*
      print*,'iostat code =',ios
      print*
 
      STOP 'program aborted in io-warning.for (CHARM)  @1054'


      end subroutine io_warning

