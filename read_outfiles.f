      subroutine read_outfiles
      
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
     
!     rev. 10.1.76 Apr.  05/17  - NK: Reorganized the outfiles.* file
      
      use area_watflood
      USE EF_module
      save
      logical exists
      integer i,j,ios
      
      ! THE OUTFILES.TXT FILE READS THE NAMES OF ALL THE OUTPUT FILES
      outfileflg='n'    ! used in sub
      ioflg=0
      INQUIRE(FILE='outfiles.txt',EXIST=exists)
      IF(exists)THEN
	  outfileflg='y'
        open(unit=99,file='outfiles.txt',
     *               status='old',iostat=ios)
        if(ios.eq.0)then
          print*
!         AN OUTFILES.TXT FILE EXISTS AND WE'LL READ IT:
          print*,'reading outfile names'
          read(99,*)  ! header line
          do i=26,30
            ioflg=ioflg+1
            read(99,5001,iostat=ios)j,outfln(i)
d            print*,i,j,outfln(i)(1:50)
            if(ios.ne.0)then 
              print*,'Problems on unit 99'
              print*,'Warning: error reading file name outfiles.txt'
              print*,'possible cause: existing file is read-only'
              print*,'or end of file is reached - need 76 lines'
              print*,'iostat code =',ios
	        print*,'Got as far as line ',ioflg
              STOP 'program aborted in spl.for @ 810'
            endif
          end do
          do i=51,80
            ioflg=ioflg+1
            read(99,5001,iostat=ios)j,outfln(i)
d            print*,i,j,outfln(i)(1:50)
            if(ios.ne.0)then 
              print*,'Problems on unit 99'
              print*,'Warning: error reading file name outfiles.txt'
              print*,'possible cause: existing file is read-only'
              print*,'or end of file is reached'
              print*,'iostat code =',ios
	        print*,'Got as far as line ',ioflg
              STOP 'program aborted in spl.for @ 771'
            endif
          end do
          do i=90,100
            ioflg=ioflg+1
            read(99,5001,iostat=ios)j,outfln(i)
d            print*,i,j,outfln(i)(1:50)
            if(ios.ne.0)then 
              print*,'Problems on unit 99'
              print*,'Warning: error reading file name outfiles.txt'
              print*,'possible cause: existing file is read-only'
              print*,'or end of file is reached - need 76 lines'
              print*,'iostat code =',ios
	        print*,'Got as far as line ',ioflg
              STOP 'program aborted in spl.for @ 797'
            endif
          end do
          do i=901,900+classcount
            ioflg=ioflg+1
            read(99,5001,iostat=ios)j,outfln(i)
d            print*,i,j,outfln(i)(1:50)
            if(ios.ne.0)then 
              print*,'Problems on unit 99'
              print*,'Warning: error reading file name outfiles.txt'
              print*,'possible cause: existing file is read-only'
              print*,'or end of file is reached - need 76 lines'
              print*,'iostat code =',ios
	        print*,'Got as far as line ',ioflg
              STOP 'program aborted in spl.for @ 784'
            endif
          end do
          close(unit=99)
!         opened in watbal.for
          fln(61)=outfln(61)
          fln(67)=outfln(67)           !used in write_r2s
          fln(70)=filename(70)             !used in write_tb0
	    fln(72)=outfln(72)
	    fln(79)=outfln(79)
	    fln(100)=outfln(100)
          filename(96)=outfln(96)
          filename(97)=outfln(97)
          fln(100)=outfln(100)             !used in write_r2c
          print*,'finished reading outfiles.txt'
          print*
        else
          write(98,*)'Error: opening outfiles.txt'
	    write(98,*)'Error: Continuing using default output files'
        endif 
      else
       write(98,*)'Info: outfiles.txt file not found, defaults used'
      endif
 5001 format(i5,a256)

      return
      
      END SUBROUTINE read_outfiles
