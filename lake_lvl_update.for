      SUBROUTINE lake_lvl_update(jz,date,time)

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
     
        
      use area_watflood
	implicit none
      save
      integer(4)     :: i,j,n,l,jz,ios,iAllocate
      integer(4)     :: fileDD,fileMM,fileYY
      real(4)        :: div,thr,at,dtmin,time
      CHARACTER(14)  :: date
      character(256) :: line,junk
      logical        :: firstpass,exists
      
      data firstpass/.true./
      
      if(firstpass)then
                allocate(lake_release(noresv,366),stat=iAllocate)
                if(iAllocate.ne.0) STOP
     *			'Allocation Error: lake_release in lake_lvl_update @ 38'
                do l=1,noresv
                    do j=1,366
                        lake_release(l,j)=-999.0
                    end do
                end do
      endif
        
            if(mod(jz,24)+6.eq.6)then
!             update 1 am local time for Ont                
              if(month_now.lt.10.and.day_now.lt.10)then
                write(line,77771)year_now,month_now,day_now
77771           format('level\',i4,'0'i1,'0',i1,'_ull.txt')
              elseif(month_now.lt.10.and.day_now.ge.10)then
                write(line,77772)year_now,month_now,day_now
77772           format('level\',i4,'0',i1,i2,'_ull.txt')
              elseif(month.ge.10.and.day_now.lt.10)then
                write(line,77773)year_now,month_now,day_now
77773           format('level\',i4,i2,'0',i1,'_ull.txt')
              else
                write(line,77774)year_now,month_now,day_now
77774           format('level\',i4,i2,i2,'_ull.txt')
              endif
              fln(99)=line
c              print*,fln(99)(1:40)
c              pause 77788
              continue
              INQUIRE(FILE=fln(99),EXIST=exists)
              if(exists)then
                write(98,*)
     *          'Info: Updating the lake levels with fln',fln(99)(1:40)
                open(unit=99,file=fln(99),status='old',iostat=ios)
                do while(.not.eof(99))
                    read(99,*,iostat=ios)l,junk,lake_elv(l,jz)
                    if(ios.eq.0)then
                      i=ires(l)
                      j=jres(l)
                      n=s(i,j)
                      store2(n)=lake_area(l)*(lake_elv(l,jz)-b7(l))+
     *                                           store_dead(l)
                      store1(n)=store2(n)
                    endif
                end do
               write(98,*)'Info: lake level updated with file:',
     *                              fln(99)(1:72)
              endif

                
              do l=1,noresv
                  lake_release(l,day_now)=-999.0
              end do
!             Read a single day of pre-determined lake outflow                
              if(month_now.lt.10.and.day_now.lt.10)then
                write(line,87771)year_now,month_now,day_now
87771           format('resrl\',i4,'0'i1,'0',i1,'_ufl.txt')
              elseif(month_now.lt.10.and.day_now.ge.10)then
                write(line,87772)year_now,month_now,day_now
87772           format('resrl\',i4,'0',i1,i2,'_ufl.txt')
              elseif(month.ge.10.and.day_now.lt.10)then
                write(line,87773)year_now,month_now,day_now
87773           format('resrl\',i4,i2,'0',i1,'_ufl.txt')
              else
                write(line,87774)year_now,month_now,day_now
87774           format('resrl\',i4,i2,i2,'_ufl.txt')
              endif
              fln(99)=line
              INQUIRE(FILE=fln(99),EXIST=exists)
              if(exists)then
              print*,fln(99)(1:40)
                print*,'NEW <<<<<'
                write(98,*)'Info:
     *              Updating the lake levels with fln',fln(99)(1:40)
                open(unit=99,file=fln(99),status='old',iostat=ios)
                do while(.not.eof(99))
                    read(99,*,iostat=ios)l,junk,lake_release(l,day_now)
                end do
                write(98,*)'Info: lake discharge updated '
              endif
          endif
      firstpass=.false.  
      return      
        
      END SUBROUTINE lake_lvl_update

      