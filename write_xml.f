      subroutine write_xml
      
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
     
!     rev. 10.1.64 Jan.  26/17  - NK: Added XML output file 
      
      use area_watflood
      use areacg
      implicit none

!     SAVES THE LOCAL VARIABLES FROM ONE CALL TO NEXT
      SAVE

      integer        :: i,j,k,l,n,jul_day1_id1,ios,jd,strLength
      logical        :: firstpass
      character(80)  :: line,tempname
      character(14)  :: str
      
      data firstpass/.true./

      if(firstpass)then
          
c        jul_day1_id1=jul_day1  
        jul_day1_id1=jul_day_now 
        
        do l=1,no
          if(l.le.9)then
            write(line,10013)l
10013       format('sta_00',i1,'.xml')
          elseif(l.le.99)then
            write(line,10014)l
10014       format('sta_0',i2,'.xml')
          elseif(l.le.999)then
            write(line,10015)l
10015       format('sta_',i3,'.xml')
          endif
          open(unit=1000+l,file=line(1:11),iostat=ios)
        end do  
        do l=1,no
          if(l.le.9)then
            write(line,10033)l
10033       format('header_sta_00',i1,'.xml')
          elseif(l.le.99)then
            write(line,10034)l
10034       format('header_sta_0',i2,'.xml')
          elseif(l.le.999)then
            write(line,10035)l
10035       format('header_sta_',i3,'.xml')
          endif
          open(unit=2000+l,file=line(1:18),iostat=ios)
        end do  
      endif     ! firstpass
      
      do l=1,no
!       jul_day1 is the first jul day this event
        jd=jul_day1
c         do k=24,nl,24
         do k=deltat_report,nl,deltat_report
           hh=mod(k-1,24)  
           if(hh.le.9)then             
           if(qsyn(l,k).lt.0.00000)then
!            missing values               
             write(1000+l,10020)'<event date="',yyyymmdd12(jd),
     *        '" time="0',hh,':00:00" value="',-999.000,'" flag="0"/>'
10020        format(8x,a13,a10,a9,i1,a15,f8.3,a12)       
           elseif(qsyn(l,k).lt.10.00000)then
             write(1000+l,10021)'<event date="',yyyymmdd12(jd),
     *        '" time="0',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
10021        format(8x,a13,a10,a9,i1,a15,f5.3,a12)       
           elseif(qsyn(l,k).lt.100.00000)then
             write(1000+l,10022)'<event date="',yyyymmdd12(jd),
     *        '" time="0',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
10022        format(8x,a13,a10,a9,i1,a15,f6.3,a12)       
           elseif(qsyn(l,k).lt.1000.00000)then
             write(1000+l,10023)'<event date="',yyyymmdd12(jd),
     *        '" time="0',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
10023        format(8x,a13,a10,a9,i1,a15,f7.3,a12)       
           elseif(qsyn(l,k).lt.10000.00000)then
             write(1000+l,10024)'<event date="',yyyymmdd12(jd),
     *        '" time="0',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
10024        format(8x,a13,a10,a9,i1,a15,f8.3,a12)       
           elseif(qsyn(l,k).lt.100000.00000)then
             write(1000+l,10025)'<event date="',yyyymmdd12(jd),
     *        '" time="0',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
10025        format(8x,a13,a10,a9,i1,a15,f9.3,a12)       
           endif
           else   ! HH > 9
           if(qsyn(l,k).lt.0.000)then
!            missing values               
             write(1000+l,20020)'<event date="',yyyymmdd12(jd),
     *         '" time="',hh,':00:00" value="',-999.000,'" flag="0"/>'
20020        format(8x,a13,a10,a8,i2,a15,f8.3,a12)       
           elseif(qsyn(l,k).lt.10.00000)then
             write(1000+l,20021)'<event date="',yyyymmdd12(jd),
     *         '" time="',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
20021        format(8x,a13,a10,a8,i2,a15,f5.3,a12)       
           elseif(qsyn(l,k).lt.100.00000)then
             write(1000+l,20022)'<event date="',yyyymmdd12(jd),
     *         '" time="',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
20022        format(8x,a13,a10,a8,i2,a15,f6.3,a12)       
           elseif(qsyn(l,k).lt.1000.00000)then
             write(1000+l,20023)'<event date="',yyyymmdd12(jd),
     *         '" time="',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
20023        format(8x,a13,a10,a8,i2,a15,f7.3,a12)       
           elseif(qsyn(l,k).lt.10000.00000)then
             write(1000+l,20024)'<event date="',yyyymmdd12(jd),
     *         '" time="',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
20024        format(8x,a13,a10,a8,i2,a15,f8.3,a12)       
           elseif(qsyn(l,k).lt.100000.00000)then
             write(1000+l,20025)'<event date="',yyyymmdd12(jd),
     *         '" time="',hh,':00:00" value="', qsyn(l,k),'" flag="0"/>'
20025        format(8x,a13,a10,a8,i2,a15,f9.3,a12)       
           endif
           endif   ! hh
           
           if(mod(k,24).eq.0)jd=jd+1
         end do
        if(id.eq.ni)then
          write(1000+l,10016)'</series>'
10016     format(4x,a9)
c          if(l.eq.no)then
c              write(1000+l,10018)'</TimeSeries>'
c10018         format((a),$)                   !  $ for no carriage return
c          endif
          close(unit=1000+l,status='keep')
        endif           
      end do

!     Last event     
!     Now that we know the end date we can write the header           
      if(id.eq.ni)then    !  needs to be writte at the end as the end date is needed in the header
!       write the charm.xml file
!       check to see that flows are averaged for the by having dt = reporting time step      
        if(kt.ne.deltat_report)then
          print*
          print*,'WARNING'
          print*,'The reporting time step must be 24 hours to write'
          print*,'the xml file. Please change this in the event file'
          print*,kt,deltat_report
          print*
          print*
        endif
          
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
!       write the header
!       can't do this until we're done
        do l=1,no
c          write(2000+l,10001)'<?xml version="1.0" encoding="UTF-8"?>'

c10001     format(a38)          
c          write(2000+l,10001)'<TimeSeries xmlns:xsi="http://www.w3.org',
c     *                       '/2001/XMLSchema-instance" xmlns:xsd="htt',
c     *                       'p://www.w3.org/2001/XMLSchema" xmlns="ht',
c     *                       'tp://www.wldelft.nl/fews/PI">'
          if(l.eq.1)then
              write(2000+l,10101)'<TimeSeries' 
10101         format(a11)        
              write(2000+l,10102)
     *           'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' 
10102         format(4x,a53)        
              write(2000+l,10103)
     *             'xmlns:xsd="http://www.w3.org/2001/XMLSchema"' 
10103         format(4x,a44)        
              write(2000+l,10104)
     *             'xmlns="http://www.wldelft.nl/fews/PI">'
10104         format(4x,a38)
              write(2000+l,10107)'<timeZone>0.0</timeZone>'
          endif
10107     format(4x,a24)     
          write(2000+l,10105)'<series>'
10105     format(4x,a8)        
          write(2000+l,10106)'<header>'
10106     format(8x,a8)        
          write(2000+l,10004)'<type>instantaneous</type>'
10004     format(12x,a26)
          strLength = LEN_TRIM(gage(l))
          write(2000+l,10005)
     *             '<locationId>',gage(l),'</locationId>'
10005     format(12x,a12,a<strLength>,12a)        
          write(2000+l,10006)'<parameterId>SQIN</parameterId>'    !<<<<fix flow
10006     format(12x,a31)        
          write(2000+l,10007)
     *             '<timeStep unit="second" multiplier="3600"/>'
10007     format(12x,a43)        
          write(2000+l,10008)'<startDate date="',startdateXML,
     *             '" time="00:00:00"/>'
10008     format(12x,a17,a10,a19)        
          write(2000+l,10009)'<endDate date="',date_now,
     *             '" time="00:00:00"/>'        
10009     format(12x,a15,a10,a19)        
          write(2000+l,10010)'<missVal>-999</missVal>'
10010     format(12x,a23)        
          strLength = LEN_TRIM(gage(l))
          write(2000+l,10011)
     *             '<stationName>',gage(l),'</stationName>'
10011     format(12x,a13,a<strLength>,a14)     
          write(2000+l,10012)'<units>CMS</units>'
10012     format(12x,a18)        
          write(2000+l,10044)'</header>'
10044     format(8x,a9)        
          close(unit=2000+l,status='keep')
        end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
        
!     rev. 10.2.10 Nov.  04/17  - NK: Fixed XML file    
        open(unit=99,file='toFews\watflood.xml',status='unknown',
     *                      iostat=ios)
        print*,'Opened  toFews\watflood.xml'
        
        do l=1,no
          if(l.le.9)then  
              
            write(tempname,10200)'header_sta_00',l,'.xml'
10200       format(a13,i1,a4) 
c            print*,tempname(1:18)
            open(unit=100,file=tempname,status='unknown',iostat=ios)
            do while(.not.eof(100))
                read(100,10201,iostat=ios)line      ! header
10201           format(a80)                
                strLength = LEN_trim(line)
                if(ios.eq.0)write(99,10202)line
10202           format(a<strLength>)     
            end do
            close(unit=100,status='keep',iostat=ios)
            strLength = LEN_TRIM(gage(l))
            write(tempname,10203)'sta_00',l,'.xml'
10203       format(a6,i1,a4) 
            open(unit=100,file=tempname,status='unknown',iostat=ios)
            do while(.not.eof(100))
                read(100,10201,iostat=ios)line      ! data
                strLength = LEN_trim(line)
                if(ios.eq.0)write(99,10202)line
            end do
            close(unit=100,status='keep',iostat=ios)
            print*,'Closed ',tempname(1:18)
            
          elseif(l.le.99)then  

            write(tempname,10204)'header_sta_0',l,'.xml'
10204       format(a12,i2,a4) 
c            print*,tempname(1:18)
            open(unit=100,file=tempname,status='unknown',iostat=ios)
            do while(.not.eof(100))
                read(100,10201,iostat=ios)line      ! header
                strLength = LEN_trim(line)
                if(ios.eq.0)write(99,10202)line
            end do
            close(unit=100,status='keep',iostat=ios)
            strLength = LEN_TRIM(gage(l))
            write(tempname,10206)'sta_0',l,'.xml'
10206       format(a5,i2,a4) 
            open(unit=100,file=tempname,status='unknown',iostat=ios)
            do while(.not.eof(100))
                read(100,10201,iostat=ios)line      ! data
                strLength = LEN_trim(line)
                if(ios.eq.0)write(99,10202)line
            end do
            close(unit=100,status='keep',iostat=ios)
            print*,'Closed ',tempname(1:18)
              
          elseif(l.le.999)then  

            write(tempname,10208)'header_sta_',l,'.xml'
10208       format(a11,i3,a4) 
c            print*,tempname(1:18)
            open(unit=100,file=tempname,status='unknown',iostat=ios)
            do while(.not.eof(100))
                read(100,10201,iostat=ios)line      ! header
                strLength = LEN_trim(line)
                if(ios.eq.0)write(99,10202)line
            end do
            close(unit=100,status='keep',iostat=ios)
            strLength = LEN_TRIM(gage(l))
            write(tempname,10210)'sta_',l,'.xml'
10210       format(a4,i3,a4) 
            open(unit=100,file=tempname,status='unknown',iostat=ios)
            do while(.not.eof(100))
                read(100,10201,iostat=ios)line      ! data
                strLength = LEN_trim(line)
                if(ios.eq.0)write(99,10202)line
            end do
            close(unit=100,status='keep',iostat=ios)
            print*,'Closed ',tempname(1:18)
              
          endif   
        end do
        write(99,10018)'</TimeSeries>'
10018   format((a),$)                   !  $ for no carriage return
        close(unit=99,status='keep')
        print*,'Closed    toFews\watflood.xml'
        call execute_command_line("del *.xml")   !, wait=.false.)        
      endif     ! last event
        
      firstpass=.false.
        
      return
      
      end subroutine write_xml     
