!***********************************************************************
      subroutine wqread
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
! this subroutine reads in the water quality data file
! lflv - oct/98
      
      use area_watflood
	implicit none

      integer :: lpos,ios,i,j,n,l

! open and read *.wqd file (*=file name from fln(1), .wqd=water quality data)
!      lpos=index(fln(1),'.')
!      fln(99)=fln(1)(1:lpos)//'wqd'
!      fln(1)='basin\duff.wqd'
!      open(unit=99,file=fln(1),status='old',iostat=ios)
!      if(ios.ne.0)then
!        print*,' Error opening file name ',fln(1)
!        print* 
!        stop 'Program aborted in WQread.for @ 22'
!      endif


      Print*,' Reading the water quality file  ',fln(26)


! sediment data...
      read(256,*)
      read(256,*)gamma,ro,viskin,grav,a_wq,b_wq
!       changed a & b to a_wq & b_wq

      read(256,*)(gc(i),i=1,classcount)
      read(256,*)(cf(i),i=1,classcount)

!     particle size - d50
      read(256,*)
      do 10 i=ycount,1,-1
         read(256,*)(diam(i,j),j=1,xcount)
   10 continue
   
!     specific weight - spg
      read(256,*)
      do 20 i=ycount,1,-1
         read(256,*)(spew(i,j),j=1,xcount)
   20 continue
   
!     erodibility - erod
      read(256,*)
      do 30 i=ycount,1,-1
         read(256,*)(erodi(i,j),j=1,xcount)
   30 continue
      
!     put into vector format
      do 50 n=1,naa
         i=yyy(n)
         j=xxx(n)
         d50(n)=diam(i,j)
         spg(n)=spew(i,j)
         erod(n)=erodi(i,j)
   50 continue

! nutrient data...
      read(256,*)
      read(256,*)ncrn,ndec,pdec,sdep   
      read(256,*)nscn,ncpw,nrec,nlec
      read(256,*)pscn,pcpw,prec,plec
   
!     fertilizer application

!     initialize values to zero
      do i=1,ycount
        do j=1,xcount
          mnfer(i,j)=0.0
          mpfer(i,j)=0.0
          mnfa(i,j)=0.0
          mpfa(i,j)=0.0
        end do
      end do

      read(256,*)
      read(256,*)nofer
      do 60 l=1,nofer
         read(256,*)i,j,mnfer(i,j),mpfer(i,j),mnfa(i,j),mpfa(i,j)
   60 continue

      print*,' Done reading '

!     put into vector format
      do 70 n=1,naa
         i=yyy(n)
         j=xxx(n)
         nfer(n)=mnfer(i,j)
         pfer(n)=mpfer(i,j)
         nfa(n)=mnfa(i,j)
         pfa(n)=mpfa(i,j)
   70 continue
      
      return

!     water quality data file not found:
99000 write(*,*)' error opening or reading wqd file - file not found'
      pause 'abnormal ending in spl #1 - hit any key to continue'
      stop  'program aborted in wqread - file not found'    

      end
