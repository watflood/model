      subroutine dds_options(dds_error)
      
!***********************************************************************
!    Copyright (C) 2003 by Nicholas Kouwen  
        
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
      

!     rev. 10.1.77 Apr.  17/17  - NK: Moved DDS err calcs to new dds_code s/r's
      
!     DDS only <<<
!     DDS only <<<
!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!     error calculation for dds with pre-emption for each event
      
      use area_watflood
      use areacg
      USE EF_module
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE
      
      integer    ::  i,j,k,l,n,ii,ios
!     DDS error functions
     	real*4    ::   log_qhyd,log_qinfl,log_mean_obs           
      real*4    ::   sum_swe
      real*4    ::   score
      real*4    ::   dds_penalty
      real*4    ::   temp_value,sum_sq_error,
     *               temp_value_sum,dds_error,
     *               error_2H,error_18O 
      CHARACTER(10) :: ctime
      CHARACTER(8)  :: cday

      data dds_penalty/1.0/
      
!       use abs(dds_flag) so it works for DDS & sensitivity runs
!     DDS only <<<
!     DDS only <<<
!     rev. 9.7.00  May.  26/10  - NK: dds with pre-emption
!     rev. 9.7.19  Jan.  18/11  - NK: Added sensitivity analysis
!     error calculation for dds with pre-emption for each event
!       use abs(dds_flag) so it works for DDS & sensitivity runs
      
c        if(dds_flag.eq.1.and.id.gt.idskip)then
      if(errflg.eq.1)then     ! weighted MSE          
!           calculate the mean squared error for each chosen station
            do l=1,no
              if(nopt(l).eq.1)then
	          sum_sq_error=0.0
                i=0
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
!                 note: qsyn -1 when not in watershed
!     rev. 9.7.08  Sep.  21/10  - NK: revised mean squared error weighting for DDS
                    temp_value=(qhyd(l,k)-qsyn(l,k))
                    sum_sq_error=sum_sq_error+temp_value*temp_value
	              i=i+1
                  endif
                end do
	          if(i.gt.0)mse(l)=sum_sq_error*sta_weight(l)
	          dds_error=dds_error+mse(l)
              endif
            end do
            
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
		  ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
c              if(inbsnflg(l+no).eq.1)then
                if(nopti(l).eq.1)then
                  do k=ktri,mhtot,ktri
                    if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                      temp_value=(qinfl(l,k)-qdwpr(l,k))
                      sum_sq_error=sum_sq_error+temp_value*temp_value
                      i=i+1
                    endif
                  end do
	          if(i.gt.0)mse(l)=sum_sq_error*sta_weight(l+no)
	          dds_error=dds_error+mse(l)
                endif
c              endif
            end do
            
      elseif(errflg.eq.2)then  ! unweighted SSE
            i=0
            do l=1,no
              if(nopt(l).eq.1)then
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
                    temp_value=(qhyd(l,k)-qsyn(l,k))
                    dds_error=dds_error+temp_value*temp_value
                    i=i+1
                  endif
                 end do
              endif
            end do
            
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
c              if(inbsnflg(l+no).eq.1)then
                if(nopti(l).eq.1)then
                  do k=ktri,mhtot,ktri
                    if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                      temp_value=(qinfl(l,k)-qdwpr(l,k))
                      dds_error=dds_error+temp_value*temp_value
                      i=i+1
                    endif
                  end do
                endif
c              endif
            end do
            
      elseif(errflg.eq.3.or.errflg.eq.11)then  ! sse weighted with mean flow
!     rev. 10.1.89 May   25/17  - NK: Added errflg = 11 for isotope DDS
        
            do l=1,no
              if(nopt(l).eq.1.and.mean_observed(l).gt.0.0)then
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
!     rev. 9.7.03  Jun.  24/10  - NK: normalized SSE with station Qmean**2
c                  sum_sq_error=sum_sq_error+(qhyd(l,k)-qsyn(l,k))**2/
c     *          			  mean_observed(l)/mean_observed(l)
                    temp_value=(qhyd(l,k)-qsyn(l,k))/mean_observed(l)
                    dds_error=dds_error+temp_value*temp_value
                  endif
                 end do
              endif
            end do
       
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            if(noresvi.gt.0)then
            do l=1,noresvi
c              if(inbsnflg(l+no).eq.1)then
                if(nopti(l).eq.1.and.mean_observed(no+l).gt.0.0)then
                  do k=ktri,mhtot,ktri
                    if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                      temp_value=(qinfl(l,k)-qdwpr(l,k))
     *                 					/mean_observed(no+l)
                      dds_error=dds_error+temp_value*temp_value
                      i=i+1
                    endif
                  end do
                endif
c              endif
            end do
            endif
            
            
c!     rev. 9.8.21  Jun.  18/12  - NK: Added swe observed date & report
c            if(courseflg.and.id.eq.ni)then
c!             this is a big penalty!!            
c              write(949,*)'penalty',
c     *                   dds_error,swe_penalty,dds_error*swe_penalty
c              dds_error=dds_error*swe_penalty
c            endif
            
      elseif(errflg.eq.4)then  ! optimize on volume only

              do l=1,no
              if(nopt(l).eq.1)then
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
                    qhyd_sum(l)=qhyd_sum(l)+qhyd(l,k)
                    qsyn_sum(l)=qsyn_sum(l)+qsyn(l,k)
                    num_obs(l)=num_obs(l)+1
                  endif
                end do
              else
                qhyd_sum(l)=0.0
                qsyn_sum(l)=0.0
              endif
              end do

!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
              if(nopti(l).eq.1)then
                do k=ktri,mhtot,ktri
                  if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                    qhyd_sum(no+l)=qhyd_sum(no+l)+qinfl(l,k)
                    qsyn_sum(no+l)=qsyn_sum(no+l)+qdwpr(l,k)
                    num_obs(no+l)=num_obs(no+l)+1
                  endif
                end do
              else
                qhyd_sum(no+l)=0.0
                qsyn_sum(no+l)=0.0
              endif
            end do

            if(id.eq.ni)then
              i=0
              do l=1,no+noresvi
!                last event              
                if(num_obs(l).gt.0)then
c	            dds_error=dds_error+
c     *	                ((qhyd_sum(l)-qsyn_sum(l))/float(num_obs(l)))**2
	            dds_error=dds_error+
     *	                (qhyd_sum(l)-qsyn_sum(l))**2
                  i=i+1
                endif
              end do

!     rev. 10.1.72 Mar.  20/17  - NK: Fixed bug in sub for error_flag = 4 
              if(i.gt.0)then 
                dds_error=sqrt(dds_error)/float(i)
              else
                  print*,'Error: no data to calculate an error'
                  print*,'Program paused in sub @ 4460'
                  pause 'hit enter to continue'
              endif
            endif

      elseif(errflg.eq.5)then  ! optimize on weighted volume Dv
            do l=1,no
              if(nopt(l).eq.1)then
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
                    qhyd_sum(l)=qhyd_sum(l)+qhyd(l,k)
                    qsyn_sum(l)=qsyn_sum(l)+qsyn(l,k)
                    num_obs(l)=num_obs(l)+1
                  endif
                end do
              else
                qhyd_sum(l)=0.0
                qsyn_sum(l)=0.0
              endif
            end do

!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
              if(nopti(l).eq.1)then
                do k=ktri,mhtot,ktri
                  if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                    qhyd_sum(no+l)=(qhyd_sum(no+l)+qinfl(l,k))
                    qsyn_sum(no+l)=(qsyn_sum(no+l)+qdwpr(l,k))
                    num_obs(no+l)=num_obs(no+l)+1
                    i=i+1
                  endif
                end do
              else
                qhyd_sum(no+l)=0.0
                qsyn_sum(no+l)=0.0
              endif
            end do

            if(id.eq.ni)then
!             if in the last event:
              do l=1,no+noresvi
                if(num_obs(l).gt.0)then
                  qhyd_mean(l)=qhyd_sum(l)/num_obs(l)
                  qsyn_mean(l)=qsyn_sum(l)/num_obs(l)
                  temp_value=
     *               abs((qhyd_mean(l)-qsyn_mean(l)))/qhyd_mean(l)*100.0
	            dds_error=dds_error+temp_value
                endif
              end do
              dds_error=dds_error/float(no+noresvi)
            endif

      elseif(errflg.eq.6)then  ! Differences weighted with mean flow
            do l=1,no
              if(nopt(l).eq.1.and.mean_observed(l).gt.0.0)then
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
!     rev. 9.7.03  Jun.  24/10  
                    dds_error=dds_error
     *                     +abs(qhyd(l,k)-qsyn(l,k))/mean_observed(l)
                  endif
                end do
              endif
            end do

!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
              if(nopti(l).eq.1)then
                do k=ktri,mhtot,ktri
                  if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                    dds_error=dds_error
     *                   +abs(qhyd(l,k)-qdwpr(l,k))/mean_observed(no+l)
                  endif
                end do
              endif
            end do

      elseif(errflg.eq.7.or.errflg.eq.10.or.errflg.eq.12)then  ! Nash efficiency
            ii=0
            do l=1,no
	        if(inbsnflg(l).eq.1)then
                if(nopt(l).eq.1.and.mean_observed(l).gt.0.0)then
                  ii=ii+1
                  do k=kt,nl,kt
                    if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
                      sum_num(l)=sum_num(l)
     *        				  +(qhyd(l,k)-qsyn(l,k))**2
                      sum_den(l)=sum_den(l)
     *        				  +(qhyd(l,k)-mean_observed(l))**2
                      num_obs(l)=num_obs(l)+1
                    endif
                  end do
                endif
	        endif
            end do
           
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
              if(inbsnflg(l+no).eq.1)then
                if(nopti(l).eq.1.and.mean_observed(l+no).gt.0.0)then
                  ii=ii+1
                  do k=ktri,mhtot,ktri
                    if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                      sum_num(l+no)=sum_num(l+no)
     *         				  +(qinfl(l,k)-qdwpr(l,k))**2
                      sum_den(l+no)=sum_den(l+no)
     *        				  +(qinfl(l,k)-mean_observed(l+no))**2
                      num_obs(l+no)=num_obs(l+no)+1
                    endif
                  end do
                endif
              endif
            end do
            
            if(id.eq.ni)then
              ii=0
!             sum the errors            
              do l=1,no+noresvi
                if(num_obs(l).gt.nl/kt/2)then     ! need at least 1/2 of the record
	            ii=ii+1
	            dds_error=dds_error+sum_num(l)/sum_den(l)
                else
                  if(inbsnflg(l))then     ! need at least 1/2 of the record
                    print*,'WARNING'
                    print*,'Number of observed flows for location ',l
                    print*,'fewer than half the record length so no'
                    print*,'error is calculated for this station'
                    print*,'Check value1 flag in the str file is 1'
                  endif
                endif
              end do
              if(ii.gt.0)then
                  dds_error=dds_error/float(ii)
              endif
            endif
            
      elseif(errflg.eq.8)then  ! Nash for log(q)
            ii=0
            do l=1,no
              if(nopt(l).eq.1.and.mean_observed(l).gt.0.0)then
                ii=ii+1
	          log_mean_obs=log(mean_observed(l))
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.001.and.qsyn(l,k).gt.0.001)then
	              log_qhyd=log(qhyd(l,k))
                    sum_num(l)=sum_num(l)
     *        				  +(log_qhyd-log(qsyn(l,k)))**2
                    sum_den(l)=sum_den(l)
     *        				  +(log_qhyd-log_mean_obs)**2
                    num_obs(l)=num_obs(l)+1
                  endif
                end do
              endif
              if(num_obs(l).gt.nl/kt/2)then
	          dds_error=dds_error+sum_num(l)/sum_den(l)
              endif
            end do
           
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
              if(inbsnflg(l+no).eq.1)then
                if(nopti(l).eq.1.and.mean_observed(l+no).gt.0.0)then
                  ii=ii+1
	          log_mean_obs=log(mean_observed(l))
                  do k=ktri,mhtot,ktri
                    if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                      log_qinfl=log(qinfl(l,k))
                      sum_num(l+no)=sum_num(l+no)
     *         				  +(log_qinfl-log(qdwpr(l,k)))**2
                      sum_den(l+no)=sum_den(l+no)
     *        				  +(log_qinfl-mean_observed(l+no))**2
                      num_obs(l+no)=num_obs(l+no)+1
                    endif
                  end do
                endif
              endif
            end do
                            
            if(id.eq.ni)then
              ii=0
!             sum the errors            
              do l=1,no+noresvi
                if(num_obs(l).gt.nl/kt/2)then     ! need at least 1/2 of the record
	            ii=ii+1
	            dds_error=dds_error+sum_num(l)/sum_den(l)
                endif
              end do
              dds_error=dds_error/float(ii)
            endif
            
      elseif(errflg.eq.9)then  ! optimize on rms of event errors
!           first calculate event sums  
            do l=1,no+noresvi
              qhyd_sum(l)=0.0
              qsyn_sum(l)=0.0
              num_obs(l)=0
            end do
            do l=1,no
              if(nopt(l).eq.1)then
                do k=kt,nl,kt
                  if(qhyd(l,k).gt.0.001.and.qsyn(l,k).gt.0.001)then
                    qhyd_sum(l)=qhyd_sum(l)+qhyd(l,k)
                    qsyn_sum(l)=qsyn_sum(l)+qsyn(l,k)
                    num_obs(l)=num_obs(l)+1
                  endif
                end do
              endif
!              store the values of each event mean 
              if(num_obs(l).gt.0)then           
                qhyd_mean_evt(l,id)=qhyd_sum(l)/num_obs(l)
                qsyn_mean_evt(l,id)=qsyn_sum(l)/num_obs(l)
              endif 
c      write(777,77777)id,l,nopt(l),num_obs(l),
c     *               qhyd_mean_evt(l,id),qsyn_mean_evt(l,id)
c77777 format(4i5,2f12.3)      
            end do
	      write(52,52001)id,(l,num_obs(l),
     *               qhyd_mean_evt(l,id),qsyn_mean_evt(l,id),l=1,no)
52001       format(i5,<l>(2i10,2f10.3))
            write(52,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~'  
            if(id.eq.ni)then
!             we're done & calculate rms of event errors  
	        temp_value_sum=0.0
              do i=1,ni 
                temp_value=0.0         
                do l=1,no
                  if(nopt(l).eq.1)then
                    temp_value=
     *                 (qhyd_mean_evt(l,i)-qsyn_mean_evt(l,i))**2
     *                        /(mean_observed(l)*mean_observed(l))   
	              if(temp_value.gt.10)then
                      print*,'for event=',i,'sta=',l,
     *                   'probable bad data ignored!!!'
	              else
                      temp_value_sum=temp_value_sum+temp_value
	              endif
	              write(52,*)
	              write(52,*)
                    write(52,*)i,l,'temp_value=',
     *                               temp_value,temp_value_sum
                  endif
                end do
                dds_error=dds_error+temp_value_sum
              end do
            endif
            write(52,*)dds_error
            dds_error=sqrt(dds_error/float(no))/float(ni)
            write(52,*)'sqrt',dds_error
            write(52,*)'______________________________________________'

      elseif(errflg.eq.13)then  ! NK  Kling-Gupta Efficiency 
                
!     rev. 10.4.26 Aug.  26/20  = NK Added Kling–Gupta efficiency (KGE) score
            continue   ! error to be calculated in stats()
                
!     REV. 10.1.21 Jan.  22/16  - NK: isotope updates
      elseif(errflg.eq.100)then  ! TH: Trying isotopes + Nash efficiency
            ii=0
            do l=1,no
	        if(inbsnflg(l).eq.1)then
                if(nopt(l).eq.1.and.mean_observed(l).gt.0.0)then
                  ii=ii+1
                  do k=kt,nl,kt
                    if(qhyd(l,k).gt.0.0000.and.qsyn(l,k).gt.0.0000)then
                      sum_num(l)=sum_num(l)
     *        				  +(qhyd(l,k)-qsyn(l,k))**2
                      sum_den(l)=sum_den(l)
     *        				  +(qhyd(l,k)-mean_observed(l))**2
                      num_obs(l)=num_obs(l)+1
                    endif
                  end do
                endif
	        endif
            end do
           
!     rev. 9.8.16  Mar.  21/11  - NK: reinstate reservoir inflow error for dds
            ktri=min(ktri,mhtot)  !  ensure deltat .le. event length 
            do l=1,noresvi
              if(inbsnflg(l+no).eq.1)then
                if(nopti(l).eq.1.and.mean_observed(l+no).gt.0.0)then
                  ii=ii+1
                  do k=ktri,mhtot,ktri
                    if(qinfl(l,k).gt.0.0.and.qdwpr(l,k).gt.0.0)then
                      sum_num(l+no)=sum_num(l+no)
     *         				  +(qinfl(l,k)-qdwpr(l,k))**2
                      sum_den(l+no)=sum_den(l+no)
     *        				  +(qinfl(l,k)-mean_observed(l+no))**2
                      num_obs(l+no)=num_obs(l+no)+1
                    endif
                  end do
                endif
              endif
            end do

            if(id.eq.ni)then
              ii=0
!             sum the errors            
              do l=1,no+noresvi
                if(num_obs(l).gt.nl/kt/2)then     ! need at least 1/2 of the record
	            ii=ii+1
	            dds_error=dds_error+sum_num(l)/sum_den(l)
                endif
              end do
              dds_error=dds_error/float(ii)/5.
              ii=0
              do l=1,n18O
                if(iso_n_18O(l).gt.0)then
                  ii=ii+1
                  temp_value=temp_value+iso_rms_18O(l)/
     *                     (abs(iso_sumO_18O(l))/iso_n_18O(l))
                  if(flg2H.eq.2)temp_value=temp_value+iso_rms_2H(l)/
     *                     (abs(iso_sumO_2H(l))/iso_n_2H(l)) 
                endif
              end do
              if(flg2H.eq.2)then
                dds_error=dds_error+temp_value/2./float(ii)
              else
                dds_error=dds_error+temp_value/float(ii)
              endif
            endif

      else
            print*,'errfg =',errflg
            print*,'No error function specified'
            pause 'hit return & then "ctrl C"'
            stop 'Program aborted in sub @ 2824'
      endif
          
          
!     rev. 10.1.89 May   25/17  - NK: Added errflg = 11  for isotope DDS
!     rev. 10.1.91 May   25/17  - NK: Added errflg = 12 for isotope DDS
      if(errflg.eq.11.or.errflg.eq.12.and.id.eq.ni)then
!         add 2nd part of the objective function of isotope error
!         this is in conjunction with errflg = 3 & 7          
          temp_value_sum=0.0
          if(n18O.gt.0)then
              do j=1,n18O
                  temp_value_sum=temp_value_sum+iso_rms_18O(j)
              end do
              error_18O=temp_value_sum/float(n18O)
          endif  
              
          temp_value_sum=0.0
          if(n2H.gt.0)then
              do j=1,n2H
                  temp_value_sum=temp_value_sum+iso_rms_2H(j)
              end do
              error_2H=temp_value_sum/float(n2H)
          endif
          if(errflg.eq.11)then
              dds_error=dds_error*error_18O*error_2H
          elseif(errflg.eq.12)then
!             This puts a little more weight on Nash's e:
!             otherwise, not enough.              
              dds_error=dds_error*dds_error*error_18O*error_2H
          endif
      endif
          
!     rev. 9.8.14  Jan.  27/11  - NK: dds_penalty added for swe not to zero in summer
!         a penalty is assigned when swe does not go to zero each summer.
          dds_error=dds_error*dds_penalty

!         write progress to dds_log.txt
      if(dds_flag.eq.1)then
            call date_and_time(cday,ctime)
            write(30,30011)id,cday(1:4),cday(5:6),cday(7:8),
     *                 ctime(1:2),ctime(3:4),ctime(5:6),
     *                 dds_error,pre_emption_value
            write(*,30011)id,cday(1:4),cday(5:6),cday(7:8),
     *                 ctime(1:2),ctime(3:4),ctime(5:6),
     *                 dds_error,pre_emption_value
30011       format(i5,5x,a4,'-',a2,'-',a2,2x,2(a2,':'),a2,2e15.6)
      endif

!         check to see if program can be aborted it error > pre-emption value
!         no pre-emption for using Nash as the objective function
!         CAN ONLY BE USED FOR MONATOMICALLY INCREASING ERROR FUCTION!!!
      if(dds_flag.eq.1)then
          if(dds_error.gt.pre_emption_value.and.errflg.le.3)then
!           kill the run - write a large value to dds\function_out.txt
            open(unit=99,file='dds\function_out.txt',
     *          status='unknown',iostat=ios)
            if(ios.ne.0)then    ! added Nov. 10/14  nk
              print*
              print*,'Unable to open file  dds\function_out.txt'
              print*,'Possible cause(s):'
              print*,'file in use by another application'
              print*,'or target directory does not exist'
              pause 'Program will abort with enter'
              stop 'Program aborted in sub.f @ 4216'
            endif
	      if(id.lt.ni)then
!             write a large error so it can't be confused by DDS for a low value
              write(99,*)dds_error*100.0
	      else
              write(99,*)dds_error
	      endif
            close(unit=99,status='keep')
!           if we get to the end of the event list, there is no pre-emption
            if(id.lt.ni)then
              print*,'Program snuffed due to pre-emption for DDS',id
            endif
            write(30,*)  !blank line between trials in dds_log.txt
!           stop here for DDS run when pre_empted
!           For sensitivity, keep going
            if(dds_flag.eq.1)stop 
          endif
      endif
        
      RETURN
          
      end subroutine dds_options

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
          
      subroutine dds_uzs(score)
      
!     rev. 10.1.77 Apr.  17/17  - NK: Moved DDS err calcs to new dds_code s/r's
!     rev. 10.1.78 Apr.  17/17  - NK: New s/r dds_UZS to calculate low flow penalty
      
!     Subroutine to calculate the UZS index for dds penalty
!     to ensure flow = close to Qlz only for very low river flows
!     Written by NK  April 16/2017      
      
      use area_watflood
      USE EF_module
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      integer    ::  i,j,k,l,n,ii,count,count_total
      integer    ::  kk(99)
      real*4     ::  sum_sq_error,temp_value,score
      Real*4     ::  UZS_retn(9999),UZSindex(99)
      logical    ::  firstpass
      
      data firstpass/.true./
      
      if(firstpass)then
          do l=1,no
!     rev. 9.7.16  Jan.  05/11  - NK: Fixed init flows outside sub-basin
              if(inbsnflg(l).eq.1)then    ! added nk Jan. 05/11
!             STORES THE GRID NO. AS FXN OF THE GAUGE(BSN) NO.
                  nn(l)=s(iy(l),jx(l))
              else
!     rev. 9.8.18  Apr.  26/11  - NK: Added in-basin check in tracer4
                  nn(l)=-1
              endif
          end do
          count=0
          count_total=1
      endif
      
!     Calculate UZS/retn for each grid   
!     assume the grid areas small & large will average out to frac = 1
!     UZS_retn is a ratio of the average ratio UZS/retn       
      do n=1,naa
          UZS_retn(n)=0.0
          do ii=1,classcount-3
              UZS_retn(n)=UZS_retn(n)+
     *                    (uzs(n,ii)/retn(ii)*(1.0-sca(n,ii))+
     *                    uzsfs(n,ii)/retn(ii)*sca(n,ii))*
     *                    aclass(n,ii)
          end do
      end do
      
      do l=1,no
          kk(l)=0
          UZSindex(l)=0.0
      end do
      
      do n=1,naa
          j=xxx(n)
          i=yyy(n)
          l=nbasin(i,j)
          UZSindex(l)=UZSindex(l)+UZS_retn(n)
          kk(l)=kk(l)+1
      end do
      
      do l=1,no
          UZSindex(l)=UZSindex(l)/float(kk(l))
      end do
      
!     calculate the average GW concentration during the low flows      
      do l=1,no
          count_total=count_total+1
          if(nopt(l).eq.1.and.UZSindex(l).le.0.4)then   ! depleted UZS
c          if(UZSindex(l).le.0.3)then   ! depleted UZS
c              if(isoconcGW(nn(l),l).lt.0.80)count=count+1
              if(isoconcGW(nn(l),l).lt.0.70)count=count+1
          endif
      end do

      score=float(count)/float(count_total)*10.0
      
C     WRITE(799,*)score
C      write(799,79900)count,count_total,
C     *       (UZSindex(l),l=1,no),score,(isoconcGW(nn(l),l),l=1,no)
C79900 format(2i10,999f8.3)      
     
      firstpass=.false.
      
      return
      
      end subroutine dds_uzs
      
