      subroutine write_modflow(ju,juold,jz)

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


!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

        DIMENSION     :: smc5(16)
        CHARACTER(14) :: date
        CHARACTER(3)  :: eofmark
        CHARACTER(1)  :: lineflg,smok,openflg
        character(20) :: junk
        REAL(4)		:: optlow,time,tot1,qwert,conv,scale,
     *               smc5,tj1,clock,t,thr,dtmin,dtmax,div,aintvl,sintvl,
     *               tot2,e1,tdum,qtemp,diff2,sdlz,
     *               wfo_spec_version_number
        INTEGER       :: rbin,inta,block,no1,classcount1,na1,ycount1,xcount1,
     *                   iallcnt1,n1,ii1,jan,m,ios,iallocate,
     *                   l,ii,juold,jj,lun,nhr,nhf,
     *                   icase,iz,jz,nchr,mz,ju,mon,i,j,n,
     *                   iwest,ieast,isouth,inorth,nocolumns,nno,k,
     *                   nu,igrdshft,minsnwflg,oldjan,flgevp22,
     *                   noresv1,nj,npick            
        CHARACTER(128):: qstr

        INTEGER(kind=2) :: result1,ntest

        CHARACTER(10) :: coordsys
!        INTEGER      :: xcount,ycount
!        REAL         :: xorigin,yorigin,xdelta,ydelta,
        real         :: a66

      real :: ha,fpw,kdn,nratio

!     write output for MODFLOW:  added May, 2000 NK
!     recharge is accumulated over 24 hour in this case
      do n=1,naa
        i=yyy(n)
        j=xxx(n)
        qrgrid(i,j)=rechrg(n)/24.
      end do
      if(ju.ne.juold.and.ju.ge.1)then
!       The grid for modflow might be different so the grid can be extended 
!       by iwest, ieast,isouth and inorth
!       Eventually these can be put into an input file
!       At this point, signs are wrong.
!        iwest=2
!        ieast=3
!        isouth=1
!        inorth=0

        iwest=0
        ieast=0
        isouth=0
        inorth=0
        
        nocolumns=xcount-iwest-ieast
        write(262,4701)ju,-1-isouth+ycount-inorth,nocolumns
!       write(*,4701)ju,-1-isouth+ycount-inorth,nocolumns
        if(nocolumns.le.10)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4710)
     *        (qrgrid(i,j),j=1+iwest,10+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.20)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4720)
     *        (qrgrid(i,j),j=1+iwest,20+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.30)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4730)
     *        (qrgrid(i,j),j=1+iwest,30+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.40)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4740)
     *        (qrgrid(i,j),j=1+iwest,40+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.50)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4750)
     *        (qrgrid(i,j),j=1+iwest,50+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.60)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4760)
     *        (qrgrid(i,j),j=1+iwest,60+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.70)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4770)
     *        (qrgrid(i,j),j=1+iwest,70+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.80)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4780)
     *        (qrgrid(i,j),j=1+iwest,80+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.90)then
          do i=ycount-inorth,1+isouth,-1
            write(262,4790)
     *        (qrgrid(i,j),j=1+iwest,90+iwest),jz,ju-1
          end do
        elseif(nocolumns.le.100)then
          do i=ycount-inorth,1+isouth,-1
            write(262,7900)
     *       (qrgrid(i,j),j=1+iwest,100+iwest),jz,ju-1
          end do
        endif

!       rechrg is a 24 hr value. It is accumulated in runof6 and 
!       then reset here as soon as it is written to the file
        do n=1,naa
          rechrg(n)=0.0
        end do

      endif  !(ju.ne.juold.and.ju.ge.1)

      return

4701  format(' Recharge in mm: ju=',i5,' rows=',i5,' columns=',i5)
4710          format(10f5.1,i5,i5,'  /jz,ju-1')
4720          format(20f5.1,i5,i5,'  /jz,ju-1')
4730          format(30f5.1,i5,i5,'  /jz,ju-1')
4740          format(40f5.1,i5,i5,'  /jz,ju-1')
4750          format(50f5.1,i5,i5,'  /jz,ju-1')
4760          format(60f5.1,i5,i5,'  /jz,ju-1')
4770          format(70f5.1,i5,i5,'  /jz,ju-1')
4780          format(80f5.1,i5,i5,'  /jz,ju-1')
4790          format(90f5.1,i5,i5,'  /jz,ju-1')
7900          format(100f5.1,i5,i5,'  /jz,ju-1')

 7500 format(i5,f10.3,i5,'  /jz,dtmin,ju  runoff in cms')
 7501 format(199f9.5)
 7502 format(i5,f10.3,i5,'  /jz,dtmin,ju  recharge in mm')
 7503 format(i5,f10.3,i5,'  /jz,dtmin,ju  baseflow in cms')




      end subroutine write_modflow