!****************************************************************************
!
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
     
!  PROGRAM: Stats
!
!  PURPOSE:  To calculate Nash, r^2, RMS, RMS/avg and %Dv using watflood output files
!
!  WRITEN BY: Angela MacLean, September 2005
!
!  Revised  NK Sept. 6/06
!
!****************************************************************************

!     rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
!     rev. 9.8.51  Mar.  11/13  - NK: Link skiplines in s/r stats to value1 in the str file
!     rev. 9.8.51  Mar.  11/13  - NK: Link skiplines in s/r stats to value1 in the str file
!     rev. 9.8.66  Jun   03/13  - NK: Added error_Dv.txt output in stats.f
!     rev. 9.9.54  Jan.  19/15  - NK: Put par & shd file names for 1st event in the headers
!     rev. 10.4.26 Aug.  26/20  = NK Added Kling–Gupta efficiency (KGE) score

c      implicit none

      subroutine Stats(unitNum,flnNum)
      
	use area_watflood
      use areacg
	implicit none


	save

      real*4, dimension(:),   allocatable :: 
     *    sumAct,sumSA,sumSAsqr,ActAvg,
     *    sumAAsqr,SimAvg,SumSSsqr,sumSim,
     *    sumSimsqr, sumActsqr, stdAct, stdSim,
     *    Correl,Dv,RMS,rNr,RMS2,
     *    APB,sumAbsSA,aMAE,bias,garAvg,
     *    sumGar,Garrick,sumAA,sumSS,sumAASS,
     *    SdevAct,SdevSim,KGE 
	integer nash_count,unitNum,flnNum
         
      real*4, dimension(:,:),allocatable :: Act,Sim     
	real*4  nash_mean,r2_mean,dv_mean,count,KGE_mean
      real*4  mean_18O,mean_2H,count_18O,count_2H

      real nsum,nash_sum,Dv_Sum,rms_sum,sumR2,sumRMS2,sumRMS

      integer, dimension(:), allocatable :: junk,numValues
      integer iflag, intValues, daily, intdays, intRemain, intTotal,
     *      intGar,anallines,hourly
	integer  iAllocate,nstations,ii,jj,kk,ios,i,j
      logical*1 exists,firstpass,printflg
      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
	character(1)   :: yn
	character(80) :: line

      data firstpass/.true./

      print*
      print*,'Calculating the statistics instats:'

      nstations = no + noresvi
!     no = # flow stations
!     noresvi = # reservoir inflows      
	anallines=int(totaltime)/kt

      if(firstpass)then
      allocate(junk(anallines),numValues(nstations),
     *  Act(nstations,anallines),Sim(nstations,anallines),
     *  sumAct(nstations),sumSA(nstations),
     *  sumSAsqr(nstations),ActAvg(nstations),
     *  sumAAsqr(nstations),SimAvg(nstations),
     *  SumSSsqr(nstations),sumSim(nstations),
     *  sumSimsqr(nstations), sumActsqr(nstations), 
     *  stdAct(nstations), stdSim(nstations),
     *  Correl(nstations),Dv(nstations),RMS(nstations),
     *  rNr(nstations),RMS2(nstations),
     *  APB(nstations),sumAbsSA(nstations),aMAE(nstations),
     *  bias(nstations),garAvg(anallines),
     *  sumGar(nstations),Garrick(nstations),
     *  sumAA(nstations),sumSS(nstations),
     *  sumAASS(nstations),
     *  SdevAct(nstations),SdevSim(nstations),
     *  KGE(nstations), 
     *  stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *         'Error with allocations  in stats @ 76'
          
      endif

      if(.NOT.allocated(R2))then 
        allocate(r2(nstations),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *         'Error with allocation of R2 in stats @ 30'
      endif

      OPEN(UNIT=unitNum,FILE=filename(flnNum), status='unknown')
!     rev. 10.2.05 Oct   28/17  - NK: Killed off stats_info.txt for iopt.ge.1
      if(iopt99)OPEN(UNIT=99,FILE='debug\stats_info.txt',
     *                          status='unknown')

      write(unitNum,*)'***************************************'
      write(unitNum,*)'*                                     *'
      write(unitNum,*)'*          WATFLOOD (R)               *'
      write(unitNum,*)'*          CHARM    (TM)              *'
      write(unitNum,*)'*                                     *'
      write(unitNum,*)'*      Version ',program_version,'             *'
      write(unitNum,*)'*      Compile date ',program_date,'        *'
      write(unitNum,*)'*                                     *'
      write(unitNum,*)'*       User conditions stated        *'
      write(unitNum,*)'*    in the WATFLOOD manual apply     *'
      write(unitNum,*)'*                                     *'
      write(unitNum,*)'*      (c) N. Kouwen, 1972-2018       *'
      write(unitNum,*)'*                                     *'
      write(unitNum,*)'***************************************'
      write(unitNum,*)

      call date_and_time(cday,time)
      write(unitNum,5026)time(1:2),time(3:4),time(5:6),
     *              cday(7:8),cday(5:6),cday(1:4)
5026  format(' Created     :      ',
     *        2(a2,':'),a2,2x,a2,'-',a2,'-',a4/)

	if(kt.eq.1)daily=2
	if(kt.eq.24)daily=1
c	skiplines=0

      if(iopt99)then
      write(99,*)'kt=',kt
	write(99,*)'no=',no
      write(99,*)'noresvi=',noresvi
      write(99,*)'daily=',daily
      write(99,*)'skiplines=',skiplines
	write(99,*)'anallines=',anallines
      endif

d     print*,'In stats'

11111 continue

	rewind 60
	rewind 73

!     rev. 9.8.47  Feb.  04/13  - NK: Headers added for spl & resin csv files
      if(hdrflg0)read(60,*,iostat=ios)(line,ii=1,no*2+1)
  
!     skip over the specified no of  record
      if(skiplines.gt.0)then
        do jj=1,skiplines
          read(60,9002,iostat=ios)junk(jj)
        end do
      endif
!      pause 11

      ios=0
      jj=1
	do ii=1,no
	  numvalues(ii)=0.0
	end do

      if(iopt99)write(99,*)'obs & computed streamflows:'
      do jj=(skiplines+1),anallines
        read(60,*,iostat=ios)
     *      kk,(Act(ii,jj),Sim(ii,jj),ii=1,no)
97001   format(i10,999f10.3)
        if(ios.ne.0)then
          print*,'end of data in line(1)',skiplines+jj
        else
          if(iopt99)write(99,97001)
     *        kk,(Act(ii,jj),Sim(ii,jj),ii=1,no)
        endif 
      end do
       
!     skip over the specified no of  record
!     the same as for the str file
      if(skiplines.gt.0)then
        do jj=1,skiplines
          read(73,*,iostat=ios)junk(jj)
        end do
      endif
      print*,'in stats - # lines in resin.txt =',anallines
!     read the lake inflow values: actual,computed
      if(noresvi.gt.0)then
        do jj=(skiplines+1),anallines
          read(73,*,iostat=ios)
     *    junk(jj),(Act(ii,jj),Sim(ii,jj),ii=no+1,nstations)
97002     format(i10,999f10.3)
          if(ios.ne.0)then
            print*,'end of data in line(1)',skiplines+jj
          else
            if(iopt99)write(99,97001)
     *        jj,(Act(ii,jj),Sim(ii,jj),ii=no+1,nstations)
          endif 
c          if(ios.ne.0)then
c            print*,'end of data in line(2)',skiplines+jj
c            go to 15
c          endif 
        end do
      endif

15    continue

c      pause 2

!     Next step is to preform the calclations to do this we need to sum 
!     the values to preform the necessary calculations 
      
c      intValues = jj-1
      intValues=anallines

      do ii=1,nstations
	  numvalues(ii)=0
!       Initialize the data for the calculationsc
c        do jj=1,intValues
          sumAct(ii) = 0.0
          sumSim(ii) = 0.0
          sumActsqr(ii) = 0.0
          sumSimsqr(ii) = 0.0
          sumSA(ii) = 0.0
          sumSAsqr(ii) = 0.0
          sumAbsSA(ii) = 0.0
	    sumAASS(ii)=0.0
          sumAA(ii) = 0.0
          sumSS(ii) =0.0
          sumAAsqr(ii) = 0.0
          sumSSsqr(ii) =0.0
          Correl(ii) = 0.0
          sumGar(ii) = 0.0
c        end do
      end do

c      pause 3

      if(iopt99)then
          write(99,*)
          write(99,*)'Sums:'
      endif
      
      do ii=1,nstations
!     Use only values if both actual & computed are > 0.0
!     Sum the data for the calculations
        do jj=1,intValues
c          if(Act(ii,jj).gt.0.0000.or.ii.gt.no)then
c          if(Act(ii,jj).gt.0.0000)then
!     rev. 10.2.64 Sep.  09/26  - NK added min_flow_cutoff for error calculations
          if(Act(ii,jj).gt.min_flow_cutoff)then
            numValues(ii) = numValues(ii) + 1
            sumAct(ii) = sumAct(ii) + Act(ii,jj)
            sumSim(ii) = sumSim(ii) + Sim(ii,jj)
            sumActsqr(ii) = sumActsqr(ii) + (Act(ii,jj))**2
            sumSimsqr(ii) = sumSimsqr(ii) + (Sim(ii,jj))**2
            sumSA(ii) = sumSA(ii) + (Sim(ii,jj)-Act(ii,jj))
            sumSAsqr(ii) = sumSAsqr(ii) + (Sim(ii,jj)-Act(ii,jj))**2
            sumAbsSA(ii) = sumAbsSA(ii) + abs(Sim(ii,jj)-Act(ii,jj))
	    endif
        end do
        if(iopt99)write(99,9901)ii,numvalues(ii),sumAct(ii),sumSim(ii),
     *   sumActsqr(ii),sumSimsqr(ii),sumSA(ii),sumSAsqr(ii),sumAbsSA(ii)
9901   format(2i10,999g15.3)
      end do

!      pause 4
      if(iopt99)then
      write(99,*)
      write(99,*)'Averages:'
      write(99,*)'  location no observations  sum_actual  sum_computed'
     *          ,'actual mean  computed mean'
      endif     
      do ii=1,nstations
!     Calclate the averages
        if(numValues(ii).gt.0)then
          ActAvg(ii)= sumAct(ii)/numValues(ii)   ! mean observed
          SimAvg(ii)= sumSim(ii)/numValues(ii)   ! mean computed
        else
          ActAvg(ii)=0.0   
          SimAvg(ii)=0.0   
	  endif
	  if(iopt99)write(99,9902)ii,numValues(ii),sumAct(ii),sumSim(ii),
     *	  ActAvg(ii),SimAvg(ii)
9902   format(2i10,999f15.3)
      end do
      
!     Finish summing the values
      do ii=1,nstations
        do jj=1,intValues
          if(Act(ii,jj).gt.0.0000)then
            sumAA(ii) = sumAA(ii) + (Act(ii,jj)-ActAvg(ii))
            sumSS(ii) = sumSS(ii) + (Sim(ii,jj)-SimAvg(ii))
	      sumAASS(ii)  = sumAASS(ii)+
     *              (Act(ii,jj)-ActAvg(ii))*(Sim(ii,jj)-SimAvg(ii))
            sumAAsqr(ii) = sumAAsqr(ii) + (Act(ii,jj)-ActAvg(ii))**2
            sumSSsqr(ii) = sumSSsqr(ii) + (Sim(ii,jj)-SimAvg(ii))**2
            
	      if(sumAAsqr(ii).gt.0.0.and.sumSSsqr(ii).gt.0.0)then
              R2(ii)=
     *        (sumAASS(ii)/(sqrt(sumAAsqr(ii))*sqrt(sumSSsqr(ii))))**2
	      else
	        r2(ii)=-9.99
            endif
          endif
c          Correl(ii) = Correl(ii) + 
c     *             (Act(ii,jj)-ActAvg(ii))*(Sim(ii,jj)-SimAvg(ii))
        end do  
      end do



      
!      pause 15

!     Now that the sums for all the values have been calculated the stats can be calculated
!     Deviation of Runoff volumes(Dv%)
!     Root mean squared Error (RMS)
!     Nash-Sutcliffe Coefficient (rNr)
!     Correlation Coefficient squared (R2)
!     RMS/qbar(avg observed flows) (RMS2)
!     Absolute percent bias (APB)
!     Mean Absolute Error (aMAE)
!     Bias (bias)
      
      do ii=1,nstations
        if(sumact(ii).gt.0.0)then
          Dv(ii) = (sumSA(ii)/sumAct(ii))*100
        else
          Dv(ii)=-999.99
        endif
        if(numValues(ii).gt.0)then
          RMS(ii) = SQRT(sumSAsqr(ii)/numValues(ii))
	  else
          RMS(ii) = -999.0
        endif
        
!     rev. 10.2.26 Jul.  03/18  - NK: Fixed nash e calculation for denominator = 0
        if(sumAAsqr(ii).ne.0.0)then
          rNr(ii) = 1.0 - (sumSAsqr(ii)/sumAAsqr(ii))
        else
          rNr(ii) = 0.0  
        endif  
c          R2(ii)=((Correl(ii)/numValues(ii))/(stdAct(ii)*stdSim(ii)))**2
        RMS2(ii) = RMS(ii)/ActAvg(ii)
        APB(ii) = (sumAbsSA(ii)/sumAct(ii))*100
        aMAE(ii) = sumAbsSA(ii)/numValues(ii)
        bias(ii) = sumSA(ii)/numValues(ii)
        
        SdevSim(ii) =  sqrt(sumSSsqr(ii)/(numValues(ii)-1))
        SdevAct(ii) =  sqrt(sumAAsqr(ii)/(numValues(ii)-1))
        KGE(ii)=1.0-sqrt((R2(ii)-1.0)**2+
     *                   (SdevSim(ii)/SdevAct(ii)-1.0)**2+
     *                   (simAvg(ii)/actAvg(ii)-1.0)**2)
        
        if(iopt99)write(99,*)
     *        ii,numValues(ii),sumSAsqr(ii),sumAAsqr(ii),ActAvg(ii)

      end do
        
!     Print out on Screen  
      print*,'Nash-Sutcliffe Coefficient (rNr)'
      print*,'Correlation Coefficient squared (R2)'
      print*,'Root mean squared Error (RMS)'
      print*,'RMS/qbar(avg observed flows) (RMS2)'
      print*,'Deviation of Runoff volumes(Dv%)'
      print*,'Absolute percent bias (APB)'
      print*,'Bias (bias)'
      print*,'Mean Absolute Error (MAE)'
      print*,'Weight w = abs(Dv)'
      write(*,*)
      write(*,*)'No of streamflow data points ',anallines
      write(*,*)'Min flow cutoff  ',min_flow_cutoff
      write(*,*)'No of skipped data points ',skiplines
      write(*,*)'Parameter file: ',par_fln(1:60)
	write(*,*)'Precipitation scale factor =',scaleall
      write(*,*)
      write(*,6002)

 !    Print to File
      write(unitNum,*)'Deviation of Runoff volumes(Dv%)'
      write(unitNum,*)'Root mean squared Error (RMS)'
      write(unitNum,*)'Nash-Sutcliffe Coefficient (rNr)'
      write(unitNum,*)'Correlation Coefficient squared (R2)'
      write(unitNum,*)'RMS/qbar(avg observed flows) (RMS2)'
      write(unitNum,*)'Absolute percent bias (APB)'
      write(unitNum,*)'Bias (bias)'
      write(unitNum,*)'Mean Absolute Error (aMAE)'
      write(unitNum,*)'Weight w = abs(Dv)'
      write(unitNum,*)
      write(unitNum,*)'No of streamflow data points ',jj-1
      write(unitNum,*)'Min flow cutoff  ',min_flow_cutoff
      write(unitNum,*)'No of skipped data points ',skiplines
      write(unitNum,*)'Parameter file: ',par_fln(1:60)
	write(unitNum,*)'Precipitation scale factor =',scaleall
      write(unitNum,*)
      write(unitNum,6003)

      do ii=1,no
!       add error check iostat Feb. 6/13 NK          
	  if(numValues(ii).gt.0.and.SimAvg(ii).gt.0.001.and.nopt(ii).eq.1)then
          write(*,6004,iostat=ios)
     *       ii,rNr(ii),KGE(ii),R2(ii),RMS(ii),RMS2(ii),
     *                 Dv(ii),APB(ii),bias(ii),aMAE(ii)
        else
          write(*,6006,iostat=ios)ii,'no data or not selected'
	  endif
      end do

      printflg=.true.
      do ii=no+1,nstations
	  if(ii.gt.no.and.printflg)then
	    print*,'Lake/Reservoir inflows'
	    printflg=.false.
	  endif    
!       add error check iostat Feb. 6/13 NK          
	  if(numValues(ii).gt.0.and.SimAvg(ii).gt.0.001.
     *                            and.nopt(no+ii).eq.1)then
          write(*,6004,iostat=ios)
     *       ii,rNr(ii),KGE(ii),R2(ii),RMS(ii),RMS2(ii),
     *                 Dv(ii),APB(ii),bias(ii),aMAE(ii)
        else
          write(*,6006,iostat=ios)ii,'no data or not selected'
	  endif
      end do
	print*
	print*,'Please see stats.txt file for copy of above stats'
	print*

!      pause 17

	nash_count=0
	nash_sum=0.0
	rms_sum=0.0
	sumR2=0.0


!     Print to File
      do ii=1,nstations
c      do ii=1,no
	  if(numValues(ii).gt.0.and.SimAvg(ii).gt.0.001)then
	    if(ii.eq.no+1)write(unitNum,*)'Lakes/Reservoir inflows'
          write(unitNum,6004,iostat=ios),ii,rNr(ii),KGE(ii),R2(ii),
     *      RMS(ii),RMS2(ii),Dv(ii),APB(ii),bias(ii),aMAE(ii),ActAvg(ii)
     *      ,R2(ii)*abs(Dv(ii))/100.0,numValues(ii)
!         no point trying to fit hydrographs with N-S efficiency less than -1
!     rev. 9.6.01  Mar.  01/10  - NK: DDS capability added
!         calculate the modified Nash-Sutcliffe coefficient
c          if(nopt(ii).gt.0)then
            if(rNr(ii).gt.-1.0.and.rNr(ii).le.1.0)then
              nash_count=nash_count+1
              if(abs(Dv(ii)).gt.a4)then
                nash_sum=nash_sum+rNr(ii)-a3*(Dv(ii)-a4)*(Dv(ii)-a4)
	        else
!               no penalty for small Dv error
	          nash_sum=nash_sum+rNr(ii)
	        endif
              Dv_Sum=Dv_Sum+Dv(ii)*Dv(ii)
	        rms_sum=rms_sum+RMS2(ii)
	        sumR2=sumR2+R2(ii)
	        sumRMS=sumRMS+RMS(ii)
	        sumRMS2=sumRMS2+RMS2(ii)
c	        print*,nash_count,a3,a4,rNr(ii),Dv(ii)
c	        print*,nash_sum
	      endif
c	    endif
        else
          write(unitNum,6004,iostat=ios)ii
	  endif
      end do

!     This sectin added May 24, 2011  NK
      nash_mean=0.0
	r2_mean=0.0
	dv_mean=0.0
	count=0.0
      KGE_mean=0.0
      optim=0.0

	do ii=1,int(nstations)
!     rev. 10.2.25 May   27/18  - NK: Fixed nash e calculation for value1=nopt=0
	  if(numValues(ii).gt.0.and.SimAvg(ii).gt.0.001.and.nopt(ii).eq.1)then
	    count=count+1.0
	    nash_mean=nash_mean+rNr(ii)
	    r2_mean=r2_mean+R2(ii)
	    dv_mean=dv_mean+abs(dv(ii))
          optim=optim+(1.0-KGE(ii))
          KGE_mean=KGE_mean+KGE(ii)
	  endif
	end do

	nash_mean=nash_mean/count
	r2_mean=r2_mean/count
	dv_mean=dv_mean/count
      KGE_mean=KGE_mean/count
      optim=optim/count

	write(*,*)
	write(*,*)'Mean Nash   = ',nash_mean
      write(*,*)'Mean KGE    = ',KGE_mean
	write(*,*)'Mean R^2    = ',r2_mean
	write(*,*)'Mean abs(Dv)= ',dv_mean
	write(unitnum,*)
	write(unitnum,*)'Mean Nash   = ',nash_mean
      write(unitnum,*)'Mean KGE    = ',KGE_mean
	write(unitnum,*)'Mean R^2    = ',r2_mean
	write(unitnum,*)'Mean abs(Dv)= ',dv_mean
      

!     rev. 10.1.89 May   25/17  - NK: Added errflg = 11 for isotope DDS
      if(frcflg.eq.'y')then
        mean_18O=0.0
        count_18O=0.0
        do j=1,n18O
          if(iso_rms_18O(j).gt.0.0)then
              mean_18O=mean_18O+iso_rms_18O(j)
              count_18O=count_18O+1.0
          endif
        end do
        mean_18O=mean_18O/count_18O
        
        mean_2H=0.0
        count_2H=0.0
        do j=1,n2H
          if(iso_rms_2H(j).gt.0.0)then
              mean_2H=mean_2H+iso_rms_2H(j)
              count_2H=count_2H+1.0
          endif
        end do
        mean_2H=mean_2H/count_2H
        
        write(*,*)
        write(*,*)'Isotope RMS error       '
        write(*,*)
     * 'Location  # readings   error 18O    # readings   error 2H'
        do j=1,n18O
          write(*,90100)j,
     *     iso_n_18O(j),iso_rms_18O(j),iso_n_2H(j),iso_rms_2H(j)
90100     format(i8,i12,f12.2,i12,f12.2)          
        end do
        write(*,*)
      
        write(unitnum,*)  
        write(unitnum,*)'Isotope RMS error      '
        write(unitnum,*)
     * 'Location   # readings   error 18O   # readings   error 2H'
        do j=1,n18O
          write(unitnum,90100)j,
     *     iso_n_18O(j),iso_rms_18O(j),iso_n_2H(j),iso_rms_2H(j)
        end do
        write(unitnum,*)
      
	  write(*,*)'Mean_18O    = ',mean_18O
	  write(*,*)'Mean 2H     = ',mean_2H
	  write(unitnum,*)'Mean_18O    = ',mean_18O
	  write(unitnum,*)'Mean 2H     = ',mean_2H
      endif


	open(unit=99,file='debug\stats_means.csv',status='unknown')
	write(*,99999)
c	write(*,99999)nash_mean,r2_mean,dv_mean
c	write(99,99999)nash_mean,r2_mean,dv_mean
	if(iopt99)write(99,*)nash_mean,r2_mean,dv_mean
99999	format(f7.3,',',f7.3,',',f7.2)
      close(unit=99,status='keep')
      close(unit=unitNum,status='keep')
      firstpass=.false.
      
!     rev. 9.8.66  Jun   03/13  - NK: Added error_Dv.txt output in stats.f
      open(unit=99,file='error_Dv.txt',status='unknown',iostat=ios)
      if(ios.ne.0)then
        print*,'Error opening error_Dv.txt in working directory'
        print*,'File not written'
        print*
        return
      endif
      if(iopt99)then
      do ii=1,nstations
        write(99,99001)ii,dv(ii)
99001   format(i5,f10.1)        
      end do
      close(unit=99,status='keep')  
      endif
          
      
	return

6002  format(' locn  Nash   KGE       r^2     rms   rms/qbar  %Dv'
     +     '     APB     Bias     MAE')
6003  format(' locn  Nash   KGE       r^2     rms   rms/qbar  %Dv'
     +     '     APB     Bias     MAE   Qbar  w*r^2   # observed')
6004  format(i4,11f8.2,i10)
c6004  format(i4,2f7.2,2f10.2,8f7.2)
6005  format(a18,a60)
6006  format(i4,a60)
    
7001  format(i5,30f10.3)
9001  format(5X,100(f7.0,1X))
9002  format(i4,1X,100(f9.0,1X))
9003  format(10f10.0)
      
      
      end subroutine Stats

