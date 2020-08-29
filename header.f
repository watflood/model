      SUBROUTINE header()
      
!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen 
        
!    WATFLOOD(R) is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    This WATFLOOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.

!    You should have received a copy of the GNU Lesser General Public License
!    along with WATFLOOD.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
      

!********************************************************************

!       THIS SUBROUTINE PRINTS THE BANNER FOR PROGRAM CHARM

!********************************************************************
!     rev. 9.9.54  Jan.  19/15  - NK: Put par & shd file names for 1st event in the headers

      USE area_watflood
      implicit none

      CHARACTER(10) :: time
      CHARACTER(8)  :: cday
      integer  :: n,i,j,iset,istep2,l,ktt,k,ktemp,iasdf,ios,noread
     *             ,nnx
      real*4   :: datemp,qdagrd,qinit

!     SAVES THE LOCAL VARIABLES FROM ONE RUN TO NEXT
      SAVE

      if(abs(dds_flag).ne.1)then
      
!  THIS SECTION GETS THE COMPUTER'S DATE AND TIME
!  THIS INFO USED TO COME FROM HEADER.FI, HOWEVER F90 INTRINSICS
!  AND CALLS ARE DIFFERENT AND THEREFORE IT NEEDED TO BE MODFIED.
       call date_and_time(cday,time)
       
        write(6,6001)
        write(6,6002)
        write(6,6003)
        write(6,6004) 
        write(6,6005) 
        write(6,6006) 
        write(6,6007) 
        write(6,6008) 
        write(6,6009) 
        write(6,6002)
        write(6,6001)

      if(modelflg.eq.'r'.or.modelflg.eq.'l'.or.modelflg.eq.'i')then
        write(6,7001)
        write(6,7002)
        write(6,7003)
        write(6,7004) shd_fln    ! writes the name from event #1
        write(6,7005) par_fln    ! writes the name from event #1
        write(6,7006) fln(3)
        write(6,7007) fln(4)
        write(6,7008) fln(6)
        write(6,7009) fln(7)
        write(6,6010)program_version,program_date,fln(8)
        write(51,6017)
        write(51,6010) ! write version number to the simout/spl.txt file
        write(51,6017)
! NOTE: FORMATS 6011/6012 NEEDED TO BE MODIFIED TO CHARACTER TYPES FROM
!       INTEGERS FOR NEW INTRINSIC.  IF THIS DOESN'T WORK, USE VALUES
!       ARRAY FROM INTRINSIC.
        write(6,7011) time(1:2),time(3:4),time(5:6)   !,fln(10)
        write(6,7012) cday(1:4),cday(5:6),cday(7:8)   !,fln(13)
        write(6,7019) ! fln(15)
        write(6,7013) iopt,snwflg,sedflg,vapflg,smrflg,resinflg,
     *              tbcflg,resumflg,contflg,routeflg,crseflg,
     *              ensimflg,picflg,wetflg,modelflg,shdflg,trcflg,
     *              fln(19)
        write(6,7018) itype,fln(20)
        write(6,7020) fln(21)     
        write(6,7014) 
        write(6,7015) 
        write(6,7016) 
        write(6,7020) 
        write(6,7001) 
        write(6,7017)
      else

        write(6,6051)
        write(6,6052)
        write(6,6053)
        write(6,6054) 
        write(6,6055) 
        write(6,6056) 
        write(6,6057) 
        write(6,6058) 
        write(6,6059) 
        write(6,6052)
        write(6,6051)
        write(6,6060)
        write(6,*)
        write(6,6021)
        write(6,6020)

c
c        write(6,6020)fln(20)
c	  write(6,6022)



! NOTE: FORMATS 6011/6012 NEEDED TO BE MODIFIED TO CHARACTER TYPES FROM
!       INTEGERS FOR NEW INTRINSIC.  IF THIS DOESN'T WORK, USE VALUES
!       ARRAY FROM INTRINSIC.
        write(6,6010)program_version,program_date,shd_fln
        write(6,6011) time(1:2),time(3:4),time(5:6),par_fln
        write(6,6012) cday(1:4),cday(5:6),cday(7:8),fln(3)
        write(6,6019) fln(13)
        write(6,6013) iopt,snwflg,sedflg,vapflg,smrflg,resinflg,
     *              tbcflg,resumflg,contflg,routeflg,crseflg,
     *              ensimflg,picflg,wetflg,modelflg,shdflg,trcflg,
     *              frcflg,initflg,
     *              fln(6)  
        write(6,6018) itype,fln(7)
        
        write(6,6020)fln(36)     
        write(6,6014)fln(8) 
        write(6,6020)fln(10) 
        write(6,6015)fln(15) 
        write(6,6016)fln(62) 
        write(6,6020)fln(20)
        write(6,6025)
        write(6,6022)fln(21)
        write(6,6023)
        write(6,6024)
        write(6,6020)
        write(6,6021)
        write(6,6017)
        
      endif

!     TS - ADDED WATFLOOD WETLAND HEADER MESSAGE:
      if(wetflg.eq.'y'.or.ensimflg.eq.'y'.or.trcflg.eq.'y'.
     *     or.sedflg.eq.'y'.or.frcflg.eq.'y')then
!        print*
	  print*,'******************************************'
      endif

d     print*,'*      >>>>>>DEBUG MODE<<<<<<            *'

      if(ensimflg.eq.'y')then
        if(int(totaltime).gt.ireport_start
     *            .and.int(totaltime).le.ireport_end)then
           print*,'*    Writing a WATFLOOD.WFO file         *'
        endif
	  if(id.eq.1)then  !totaltime not yet defined
           print*,'*    Writing a WATFLOOD.WFO file         *'
        endif
      endif
      
      if(modelflg.eq.'n')then
        if(wetflg.eq.'y')then
          print*,'*      RUNNING WATFLOOD WETLAND          *'
        endif
        if(ruleflg)then
          print*,'*      USING LAKE TARGET LEVELS          *'
        endif
c        if(dlyflg)then
        if(flgevp2.eq.4)then
          print*,'*   RUNNING w/daily temp. differences    *'
        endif
        if(lakeflg.eq.'y')then
          print*,'*    RUNNING WATFLOOD Lake Evapn.        *'
        endif
        if(sedflg.eq.'y')then
          print*,'*      RUNNING WATFLOOD QUALITY          *'
        endif
        if(trcflg.eq.'y')then
          print*,'*    RUNNING WATFLOOD VIRTUAL TRACER      *'
          print*,'*    Tracer ',itrace, '                *'
        endif
        if(frcflg.eq.'y')then
          print*,'*    RUNNING ISOTOPE FRACTIONATION       *'
        endif
        if(ntrlflg.eq.'y')then
          print*,'*      RUNNING WITH NATURAL FLOWS        *'
        endif
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
        if(iceflg.eq.'y')then
          print*,'*  river ice effects enabled  (iceflg=y) *'
        endif
!     rev. 10.1.44 Dec.  02/15  - NK: Reworked icerivflg & icelakeflg
        if(icerivflg.eq.'y')then
          print*,'* river ice effects enabled (icerivflg=y)*'
        endif
        if(icelakeflg.eq.'y')then
          print*,'* lake ice effects enabled (icerivflg=y) *'
        endif
	endif
      if(wetflg.eq.'y'.or.ensimflg.eq.'y'.or.trcflg.eq.'y'.
     *     or.sedflg.eq.'y'.or.frcflg.eq.'y')then
        print*,'******************************************'
!	  print*
      endif
! FORMATSsplx
      
      elseif(id.eq.1)then
          
        write(6,6051)
        write(6,6010)program_version,program_date,shd_fln
        write(6,6051)
          
      endif

 6051 format(1x,'*************^^**********************************')
 6052 format(1x,'*                                               *')
 6053 format(1x,'*    cccc  hh   hh     a     rrrrr  mmm   mmm(R)*')
 6054 format(1x,'*   cc  cc hh   hh    aaa    rr  rr mmmm mmmm   *')  
 6055 format(1x,'*   cc     hh   hh   aa aa   rr  rr mm mmm mm   *')
 6056 format(1x,'*   cc     hhhhhhh  aa   aa  rrrrr  mm  m  mm   *')
 6057 format(1x,'*   cc     hh   hh  aaaaaaa  rrrrr  mm     mm   *')
 6058 format(1x,'*   cc  cc hh   hh aa     aa rr  rr mm     mm   *')
 6059 format(1x,'*    cccc  hh   hh aa     aa rr  rr mm     mm   *')
 6060 format(1x,'     Canadian Hydrological And Routing Model')




 1000 format(' ',3(2a30/))
 6001 format(1x,'*******************************************************
     ***********************')
 6002 format(1x,'*                                                      
     *                     *')
 6003 format(1x,'*  ww  ww         ww  a    tttttttt ffffff ll       ooo
     *      ooo    ddddd   *')
 6004 format(1x,'*   ww  ww       ww  aaa   tttttttt ffffff ll      oooo
     *o    ooooo   dddddd  *')
 6005 format(1x,'*    ww  ww     ww  aa aa     tt    ff     ll     oo   
     *oo  oo   oo  dd   dd *')
 6006 format(1x,'*     ww  ww   ww  aaa aaa    tt    ffff   ll     oo   
     *oo  oo   oo  dd   dd *')
 6007 format(1x,'*      ww  ww ww  aaaaaaaaa   tt    ffff   ll     oo   
     *oo  oo   oo  dd   dd *')
 6008 format(1x,'*       ww  www  aa       aa  tt    ff     llllll  oooo
     *o    ooooo   dddddd  *')
 6009 format(1x,'*        ww  w  aa         aa tt    ff     llllll   ooo
     *      ooo    ddddd   *')
 6010 format(1x,'*      ver=',2a10,'          *',2x,a30)
 6011 format(1x,'*      runtime    ',2(a2,':'),a2,15x,'*',2x,a31)
 6012 format(1x,'*      rundate  ',a4,'-',a2,'-',a2,15x,'*',2x,a30)
 6019 format(1x,'*                                        *',2x,a30)
 6013 format(1x,'*   debug level  ',i2,1x,18a1,    2x' *',2x,a30)
 6018 format(1x,'*   channel type ',i2,' 123456789012345678   *',2x,a30)
 6014 format(1x,'*             WATFLOOD(R)                *',2x,a30)
 6015 format(1x,'*    copyright (c) by n kouwen 1985-2018 *',2x,a30)
 6016 format(1x,'*    university of waterloo,  canada     *',2x,a30)
 6017 format(' ')
 6025 format(1x,'*              Written for               *',2x,a30)
 6022 format(1x,'*   E N V I R O N M E N T  C A N A D A   *',2x,a30)
 6023 format(1x,'*       Open source code under the       *')
 6024 format(1x,'*    GNU Lesser General Public License   *')
 6020 format(1x,'*                                        *',2x,a30)
 6021 format(1x,'******************************************',2x,a30)


 7001 format(1x,'******************************************',2x,a30)
 7002 format(1x,'*                                        * files used')
 7003 format(1x,'*   RRRRRR   TTTTTTTT EEEEEE    999999   *')
 7004 format(1x,'*   RR    R     TT    EE       99    99  *',2x,a30)
 7005 format(1x,'*   RR    R     TT    EE       99    99  *',2x,a30)
 7006 format(1x,'*   RRRRRR      TT    EEEE      9999999  *',2x,a30)
 7007 format(1x,'*   RR  RR      TT    EE             99  *',2x,a30)
 7008 format(1x,'*   RR   RR     TT    EE             99  *',2x,a30)
 7009 format(1x,'*   RR    R     TT    EEEEEEE   999999   *',2x,a30)
!7010 format(1x,'*      ver=x.x.xx   xxx. xx/xx           *',2x,a30)
 7011 format(1x,'*      runtime    ',2(a2,':'),a2,15x,'*',2x,a31)
 7012 format(1x,'*      rundate  ',a4,'-',a2,'-',a2,15x,'*',2x,a30)
 7019 format(1x,'*                                        *',2x,a30)
 7013 format(1x,'*     debug level  ',i2,1x,16a1,    2x' *',2x,a30)
 7018 format(1x,'*     channel type ',i2,' 1234567890123456   *',1x,a30)
 7014 format(1x,'*    copyright (c) by n kouwen 1985-2016 *',2x,a30)
 7015 format(1x,'*    department of civil engineering     *',2x,a30)
 7016 format(1x,'*    university of waterloo,  canada     *',2x,a30)
 7017 format(' ')
 7020 format(1x,'*                                        *',2x,a30)

      RETURN

      END SUBROUTINE header
