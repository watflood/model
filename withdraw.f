
      SUBROUTINE withdraw()

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
!
!   THIS S/R  takes water from storage for municipal or irrigation use
!
!***********************************************************************
!     rev. 9.8.49  Feb.  20/13  - NK: Added n=municipal & irrigation withdrawals

!     This s/r is called only once a year even if events are shorter
!     in the case of shorter events, 
!     the file is read in the first event: id=1      


      use area_watflood
	implicit none

      logical        :: exists,firstpass
      character(32)  :: tempname(3),use0
      character(256) :: line
      integer        :: no_frames
      integer        :: ios,i,j,k,l,m,n,ntake,iallocate
      real*4         :: lat0,long0,value0
      real*4         :: latS(6),longS(6),valueS(12)
	real*4         :: crop0(12),urban(12),other(12),
     *              	winter(12),summer(12)
	real*4         :: mohours(12),mosecs(12)


	data firstpass/.true./
	data crop0/0.0,0.0,0.0,0.05,0.05,0.14,0.39,0.29,0.04,0.04,0.0,0.0/
	data urban/0.06,0.06,0.06,0.06,0.06,0.13,
     *           0.13,0.13,0.13,0.06,0.06,0.06/
	data other/0.08333,0.08333,0.08333,0.08333,0.08333,0.08333,
     *           0.08333,0.08333,0.08333,0.08333,0.08333,0.08333/
	data summer/0.00,0.00,0.00,0.125,0.125,0.125,
     *           0.125,0.125,0.125,0.125,0.125,0.00/
	data winter/0.20,0.20,0.20,0.0,0.0,0.0,
     *           0.0,0.0,0.0,0.0,0.20,0.20/
      DATA mohours
     *     /744.,672.,744.,720.,744.,720.,744.,744.,720.,744.,720.,744./

      write(line,10001)year1
10001 format('diver\ab',i4,'.xyz')
      read(line,10002)tempname(1)
10002 format(a16)
      write(line,10003)year1
10003 format('diver\mb',i4,'.xyz')
      read(line,10002)tempname(2)
      tempname(3)='diver\SR12011.txt' ! SK withdrawal file name.
!     withdrawals can be made one year and not the next
      withdrawflg=.false.
      if(firstpass)then
        allocate(qwdr(na,12),stat=iAllocate)
        if(iAllocate.ne.0) STOP
     *  'Error with allocation of withdraw(,) in withdraw.f'  
      endif
!     number of seconds in each month        
      do l=1,12
        mosecs(l)=mohours(l)*3600.0
      end do
!     for leap years extra day        
      if(mod(year1,4).eq.0)then
        mosecs(2)=mosecs(2)+86400.0
      endif

!     initialize withdrawals every year      
      do n=1,na
	  do l=1,12
          qwdr(n,l)=0.0
	  end do
      end do
      
!     Alberta & Manitoba Alberta & Manitoba Alberta & Manitoba Alberta & Manitoba        
c      if(iopt.ge.1.and.id.eq.1)then
c        open(unit=51,file='withdrawal_debug.csv',
c     *                              status='unknown',iostat=ios)
c      endif
        
      do k=1,2  ! alberta & manitoba (so far)
!       Look for the data file and read if present
!       and accumulate withdrawals in each grid
        INQUIRE(FILE=tempname(k),EXIST=exists)
        IF(exists)THEN
	    if(iopt.ge.1)
     *	    print*,'withdrawal file found file name =',tempname(k)(1:32)
          if(iopt.ge.1)then
            write(51,*)'Province #',k,tempname(k)
            write(51,*)
     *          '    i    j    k    m    n  monthly wthdrwls m**3 ---->'
          endif
          withdrawflg=.true.
          open(unit=99,file=tempname(k),status='old',iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',tempname(k)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in withdraw.f @ 104'
          else 
            if(iopt.ge.1)print*,'opened ',tempname(k)(1:32)
          endif
          ios=0
          i=0
!         count the number of lines in the file & find # locations        
          do while(ios.eq.0)  
            read(99,*,iostat=ios)line
            i=i+1
          end do
          ntake=i-2
          if(iopt.ge.1)then
            print*,'No of withdrawal locations =',ntake
          endif
          rewind 99
          read(99,*)line
!         accumulate the withdrawals in m**3 in each grid        
          do m=1,ntake
!           read the withdrawal amount          
            read(99,*)long0,lat0,value0,use0,source
            if(iopt.ge.0)then
              write(51,*)long0,lat0,value0,use0,source
            endif
            i=int((lat0-yorigin)/ydelta)+1
            j=int((long0-xorigin)/xdelta)+1
            if(i.ge.1.and.i.le.ycount.and.j.ge.1.and.j.le.xcount)then
            n=s(i,j)     !  grid number
            if(n.ge.1.and.n.le.naa)then
!             populate the array with monthly withdrawals in cubic metres
!             accumulate the amounts in each grid
              if(use0(1:14).eq.'growing_season')then
	          do l=1,12
	            qwdr(n,l)=qwdr(n,l)+value0*crop0(l)
	          end do
	        elseif(use0(1:11).eq.'open_season')then
 	          do l=1,12
	            qwdr(n,l)=qwdr(n,l)+value0*summer(l)
	          end do
	        elseif(use0(1:5).eq.'urban')then
 	          do l=1,12
	            qwdr(n,l)=qwdr(n,l)+value0*urban(l)
	          end do
	        elseif(use0(1:5).eq.'winter')then
 	          do l=1,12
	            qwdr(n,l)=qwdr(n,l)+value0*urban(l)
	          end do
              else      ! default - equal monthly withdrawals
!               this includes the annual class in the files          
 	          do l=1,12
	            qwdr(n,l)=qwdr(n,l)+value0*other(l)
	          end do
	        endif
	        if(iopt.ge.1)then
	          write(51,69900)i,j,k,m,n,(qwdr(n,l),l=1,12)
69900           format(5i5,12f12.0)	 
              endif  
            endif 
	      endif  ! in the watershed bounds
          end do
	    close(unit=99,status='keep')
        
	    withdrawflg=.true.
        endif  ! if data files exists
      end do ! m=1,2 provinces 

!     convert to withdrawal flows in m**3/s
!     qwdr(n,l) is used in route where it is subtracted from 
!     channel or lake inflow along with other outflows (e.g. streamloss)
!This conversion applies only to Alberta and Manitoba, Not Sask.
!---------------------------------------------------------------
	do l=1,12
        do n=1,na
          qwdr(n,l)=qwdr(n,l)/mosecs(l)
        end do
      end do

!     Saskatchewan Saskatchewan Saskatchewan Saskatchewan Saskatchewan Saskatchewan 
      if(iopt.ge.1)then
c        write(*,*)'sask'
c        print*,tempname(3)
        write(51,*)'sask withdraw'
      endif
      INQUIRE(FILE=tempname(3),EXIST=exists)
      IF(exists)THEN
        if(iopt.ge.1)
     *    print*,'withdrawal file found file name =',tempname(3)(1:32)
        if(iopt.ge.1)then
          write(51,*)'Saskatchewan file=',tempname(3)
          write(51,*)
     *          '    i    j    k    m    n  monthly wthdrwls m**3 ---->'
        endif
        withdrawflg=.true.
        open(unit=99,file=tempname(3),status='old',iostat=ios)
        if(ios.ne.0)then    ! added Nov. 10/14  nk
          print*
          print*,'Unable to open file',tempname(3)
          print*,'Possible cause(s):'
          print*,'file in use by another application'
          print*,'or target directory does not exist'
          stop 'Program aborted in withdraw.f @ 205'
        else 
          if(iopt.ge.1)print*,'opened ',tempname(3)(1:32)
        endif
!       read the coordinates  
        if(iopt.ge.1)write(51,*)'Saskatchewan withdrawals:'      
        do i=1,5
          read(99,*)line,latS(i),longS(i)
          if(iopt.ge.0)then
            write(51,*)line,latS(i),longS(i)
          endif
        end do
!       scan the file for the year we want
        do while(.not.eof(99))
          read(99,*,iostat=ios)i
          if(i.eq.year1)then
            read(99,99000)line   !read the blank line
99000       format(a1)            
            do k=1,5
              read(99,99001,iostat=ios)line,(valueS(l),l=1,12)
	        if(ios.ne.0)then
	          print*,'problems reading sask withdrawals'
	          print*
	          stop 'program aborted in withdraw @ 207'
	        endif
99001         format(a25,12f9.3)
              write(51,*)line,(valueS(l),l=1,12)
              i=int((latS(k)-yorigin)/ydelta)+1
              j=int((longS(k)-xorigin)/xdelta)+1
              n=s(i,j)
	        do l=1,12
	          qwdr(n,l)=qwdr(n,l)+valueS(l)
	        end do
            end do
          endif
        end do
        close(unit=99)
      endif  
      
      
            
!     write the withdrawals to       
      if(withdrawflg)then
      if(iopt.ge.1)then
        do l=1,12
          write(line,10005)l
10005     format(i2,'.xyz')
          read(line,10006)tempname(1)
10006     format(a6)
          open(unit=700+l,file=tempname(1),iostat=ios)
          if(ios.ne.0)then    ! added Nov. 10/14  nk
            print*
            print*,'Unable to open file',tempname(1)
            print*,'Possible cause(s):'
            print*,'file in use by another application'
            print*,'or target directory does not exist'
            stop 'Program aborted in withdraw.f @ 259'
          endif
          write(700+l,*)'month=',l,' unit =',700+l,' ID=',id
	    do i=1,ycount
	      do j=1,xcount
              n=s(i,j)
	        if(n.gt.0)then
                if(qwdr(n,l).gt.0.0)then
	          write(700+l,70000)
     *	                    xorigin+xdelta*float(j)-xdelta/2.0,',',
     *                      yorigin+ydelta*float(i)-ydelta/2.0,',',
     *                      qwdr(n,l),n,i,j
70000           format(2(f10.3,a1),f10.3,3(',',i10))     
                endif
              endif
	      end do
	    end do
	    if(id.eq.ni)close(unit=700+l,status='keep')
	  end do
	endif
	endif

      
      
      
      
      firstpass=.false.

      return

      end subroutine withdraw