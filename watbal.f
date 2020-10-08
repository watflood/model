      SUBROUTINE watbal(jan)

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

!     REV. 9.00  - Mar.   2000 -  TS: CONVERTED TO FORTRAN 90
!     rev. 9.9.13  Apr.  04/16  - NK: Fix water balance

!***********************************************************************

      use area_watflood
	implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      integer(4)  ::  io_err,sys_err,stat,unit,cond,n,jan,l,j,i,ii,ios
      integer(4)  ::  m
      real(4)     ::  freewat(17500),balance(17500)
      integer     ::  iallocate
      real*4, dimension(:),   allocatable :: lz_frac


!     rev. 9.8.45  Jan.  31/13  - NK: disabled some writes for iopt = 99
      if(iopt.eq.99)return
      if(.not.allocated(lz_frac))then
        allocate(lz_frac(na),stat=iAllocate)
        if(iAllocate.ne.0)STOP 
     *    'Warning: error with allocation of lz_area arrays in watbal' 
      endif
!     START OF WATER BALANCE CHECK
!     ALL AMOUNTS HERE ARE FOR a whole grid

      do n=1,naa
	   totint(n)=0.0
         totd1(n)=0.0
         totuzs(n)=0.0
         totsnw(n)=0.0
         totchnl(n)=0.0
         freewat(n)=0.0
      end do

d	if(iopt.eq.2)print*,' checkpoint 1 in watbal. JAN=',jan

      if(jan.eq.1)then   !firstpass
         do n=1,naa
            storinit(n)=0.0
            sump(n)=0.0
            sumrff(n)=0.0
            do ii=1,classcount
               ssumr(n,ii)=0.0
               sumf(n,ii)=0.0
               sumffs(n,ii)=0.0	   
            end do
            eloss(n)=0.0
	     netinflow(n)=0.0
	     netoutflow(n)=0.0
         end do
         do ii=1,classcount
            sr(ii)=0.0
            sqint(ii)=0.0
            sqintfs(ii)=0.0
            sdrng(ii)=0.0
            sdrngfs(ii)=0.0
            sq1(ii)=0.0
            sq1fs(ii)=0.0           
            sexcess(ii)=0.0
         end do
         sqlz=0.0
         slzinflw=0.0    ! not used

c!     rev. 9.9.15  Jun.  02/14  - NK: Add lz to the water balance - it was missing
c!        calculate the fraction of the grid occupied by LZS
c         allocate(lz_frac(na),stat=iAllocate)
c         if(iAllocate.ne.0)STOP 
c     *    'Warning: error with allocation of lz_area arrays in watbal' 
c         do n=1,naa
c           lz_frac(n)=0.0
c           if(wetflg.eq.'y')then
c             do ii=1,classcount-3
c               lz_frac(n)=lz_frac(n)+aclass(n,ii)
c               write(813,81302)n,ii,aclass(n,ii),lz_frac(n)
c             end do
c           else
c             do ii=1,classcount-2
c               lz_frac(n)=lz_frac(n)+aclass(n,ii)
c               write(813,81302)n,ii,aclass(n,ii),lz_frac(n)
c             end do
c           endif
c           lz_frac(n)=lz_frac(n)+aclass(n,classcount)
c           write(813,81302)n,classcount,aclass(n,classcount),lz_frac(n)
c81302      format(2i5,2f15.3)
c           write(813,*)
c         end do
      endif     !  jan = 1  =  firstpass

      do n=1,na
        lz_frac(n)=1.0    ! lzs is spread over the whole grid, including water & wetlands
      end do


d	if(iopt.eq.2)print*,' checkpoint 2 in watbal. JAN=',jan

      if(wetflg.eq.'y')then
        m=classcount-2
      else
        m=classcount-1
      endif

      do n=1,naa
        do ii=1,classcount
           if(aclass(n,ii).gt.0.0)then
!              for all classes:           
               totsnw(n)=totsnw(n)
     *                  +snowc(n,ii)*aclass(n,ii)*sca(n,ii)
               if(abs(wcl(n,ii)).lt.0.0001)wcl(n,ii)=0.0  
!                     prevent underflow - 22/04/02 nk
               freewat(n)=freewat(n)+wcl(n,ii)*aclass(n,ii)*sca(n,ii)
c              freewat(n)=freewat(n)+water(n,ii)*aclass(n,ii)*sca(n,ii)
             if(ii.eq.classcount-1)then   !water
!              water water water water water water water water water             
               if(ireach(n).eq.0)then
!                not in a lake or reservoir              
c                totchnl(n)=(store2(n))/grid_area(n)*1000.*aclass(n,ii)  
c                 totchnl(n)=store2(n)/grid_area(n)*1000.    !*aclass(n,ii) 
c                 totchnl(n)=store2(n)/(grid_area(n)/frac(n))*1000.    !*aclass(n,ii) 
                 totchnl(n)=store2(n)/al/al*1000.     
               else
                 if(res(n).gt.0)then   ! outlet grid
c                   totchnl(n)=(store2(n))/grid_area(n)*1000.    !*aclass(n,ii)   
c                   totchnl(n)=(store2(n))/(grid_area(n)/frac(n))*1000.    !*aclass(n,ii)   
                   totchnl(n)=(store2(n))/al/al*1000.    !*aclass(n,ii)   
                 else
                   totchnl(n)=0.0   
                 endif
               endif
             elseif(ii.eq.classcount-2.and.wetflg.eq.'y')then 
!              coupled wetland   coupled wetland   coupled wetland 
!              convert wetland storage into depth in mm on the nominal grid
               totwetl(n)=wstore2(n)/(grid_area(n)/frac(n))*1000.   
               totint(n)=totint(n)+v(n,ii)*aclass(n,ii)
             else
!              for all classes except water and coupled wetland:  
!     rev. 10.2.52 Apr.  15/19  - NK: Added total UZS for reporting in FEWS
!              this stuff added to sub for reporting in FEWS
               totint(n)=totint(n)+v(n,ii)*aclass(n,ii)
               totd1(n)=totd1(n)
     *                 +d1(n,ii)*aclass(n,ii)*(1.0-sca(n,ii))
     *                 +d1fs(n,ii)*aclass(n,ii)*sca(n,ii)
               totuzs(n)=totuzs(n)
     *                  +uzs(n,ii)*aclass(n,ii)*(1.0-sca(n,ii))
     *                  +uzsfs(n,ii)*aclass(n,ii)*sca(n,ii)
               totwetl(n)=0.0
             endif
           endif         !  class.gt.0.0
        end do          !  ii=1,classcount
         
!       add all the storages in mm
        totgrid(n)=totd1(n)
     *              +totint(n)
     *              +totuzs(n)
     *              +totsnw(n)
     *              +freewat(n)
     *              +lzs(n)*lz_frac(n)
     *              +totwetl(n)
     *              +totchnl(n)
!       need to scale the netoutflow to nominal grid frac = 1

        netoutflow(n)=netoutflow(n)/frac(n)

c        if(n.eq.nnprint)then
c           write(812,81102)totd1(n),totint(n),totuzs(n),totsnw(n)
c     *              ,freewat(n),lzs(n)*lz_frac(n),totwetl(n),totchnl(n),
c     *               lz_frac(n)  
c81102      format(9f15.3)
c        endif   
         
      end do

d	if(iopt.eq.2)print*,' checkpoint 3 in watbal. JAN=',jan

      if(jan.eq.1)then
         do n=1,naa
            storinit(n)=totgrid(n)
            balance(n)=0.0
         end do
      else
         do n=1,naa
!          need to adjust sump because water area is NOT included in 
!          the water balance. All amounts are for a whole grid
!          but water is added only to the non-water area
!          & spread over the whole grid for this accounting.
           if(ireach(n).eq.0)then
!     rev. 9.9.13  Apr.  04/16  - NK: Fix water balance
             balance(n)=storinit(n)
     *          	    +sump(n)               ! *(1.0-aclass(n,classcount))
     *                  -eloss(n)
     *                  -netoutflow(n)
     *                  -totgrid(n)
           else
             balance(n)=storinit(n)
     *                  -totgrid(n)
     *                  -netoutflow(n)
           endif
         end do
      endif

d	if(iopt.eq.2)print*,' checkpoint 4 in watbal. JAN=',jan

      if(jan.eq.1)then
        open(unit=97,file=filename(96),status='unknown',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file',filename(96)(1:40)
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
          stop 'Program aborted in watbal.f @ 211'
        endif
      else
        open(unit=97,file=filename(97),status='unknown',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file',filename(97)(1:40)
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
          stop 'Program aborted in watbal.f @ 221'
        endif
      endif

      write(97,6200)
 6200 format
     *(' n, basin, row, col, reach, fracn,  elev,   stinit, interc,
     *  surf,   swe,       uzs,   lzs,  chlstore,  wetstore,  totstore,
     *  precip,  eloss ,netoutflow, balance, class_1_to_classcount')

      do n=1,naa
        if(ireach(n).eq.0)then
          i=yyy(n)
          j=xxx(n)
            if(n.eq.nnprint)write(97,6201,iostat=ios)
            write(97,6201,iostat=ios)
     *        n,nhyd(i,j),i,j,ireach(n),frac(n),elev(n),
     *        storinit(n),totint(n),totd1(n),totsnw(n),totuzs(n),
     *        lzs(n)*lz_frac(n),totchnl(n),totwetl(n),totgrid(n),
     *        sump(n),eloss(n),netoutflow(n),balance(n),
     *        (aclass(n,ii),ii=1,classcount)
            if(n.eq.nnprint)write(97,6201,iostat=ios)
        endif    
      end do
      
      write(97,*)
      write(97,*)'Water balance for reservoirs & lakes'
      write(97,6200)
      do n=1,naa
        if(ireach(n).gt.0.and.res(n).eq.0)then
          i=yyy(n)
          j=xxx(n)
            if(n.eq.nnprint)write(97,6201,iostat=ios)
            write(97,6201,iostat=ios)
     *        n,nhyd(i,j),i,j,ireach(n),frac(n),elev(n),
     *        storinit(n),totint(n),totd1(n),totsnw(n),totuzs(n),
     *        lzs(n)*lz_frac(n),totchnl(n),totwetl(n),totgrid(n),
     *        sump(n),eloss(n),netoutflow(n),balance(n),
     *        (aclass(n,ii),ii=1,classcount)
            if(n.eq.nnprint)write(97,6201,iostat=ios)
        endif
      end do
      
      write(97,*)
      write(97,*)'Water balance for reservoir & lake outlets'
      write(97,6200)
      do n=1,naa
        if(ireach(n).gt.0.and.res(n).gt.0)then
          i=yyy(n)
          j=xxx(n)
            if(n.eq.nnprint)write(97,6201,iostat=ios)
            write(97,6201,iostat=ios)
     *        n,nhyd(i,j),i,j,ireach(n),frac(n),elev(n),
     *        storinit(n),totint(n),totd1(n),totsnw(n),totuzs(n),
     *        lzs(n)*lz_frac(n),totchnl(n),totwetl(n),totgrid(n),
     *        sump(n),eloss(n),netoutflow(n),balance(n),
     *        (aclass(n,ii),ii=1,classcount)
            if(n.eq.nnprint)write(97,6201,iostat=ios)
        endif
      end do
 6201 format(i4,4(',',i4),',',f7.3,',',f7.1,',',f10.0,3(',',f6.1),
     *           ',',f10.3,',',f6.1,
     *           3(',',f10.1),4(',',f8.0),'  ',99(',',f6.3))



d	if(iopt.eq.2)print*,' checkpoint 5 in watbal. JAN=',jan

      close (unit=97,status='keep')


d	if(iopt.eq.2)print*,' checkpoint 5 in watbal. JAN=',jan

!      CALL ERRSNS(io_err,SYS_ERR,STAT,UNIT,cond)

d	if(iopt.eq.2)print*,' checkpoint 6 in watbal - before return'


!     * * * *  END OF WATER BALANCE CHECK  * * * * 

! FORMATS:



      RETURN

      END SUBROUTINE watbal
