  subroutine disaggregate(jz)

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
    
!    First coded March 2019 NK
!    This s/r is used for FEWS only
    
  use area_watflood
  implicit none
  save

  real,    dimension(:,:), allocatable :: unused,pLocal
  integer, dimension(:,:), allocatable :: no_hrs_precip

  integer :: nhdt,nh,n,i,j,jz,iAllocate
  logical :: firstpass
  
  data firstpass/.true./

  if(firstpass)then
      allocate(pLocal(ycount,xcount),stat=iAllocate)
      if(iAllocate.ne.0) STOP'Error with allocation of pLocal in `disaggregate`'
      allocate(no_hrs_precip(ycount,xcount),stat=iAllocate)
      if(iAllocate.ne.0) STOP'Error with allocation of no_hrs_precip in `disaggregate`'
      allocate(unused(ycount,xcount),stat=iAllocate)
      if(iAllocate.ne.0) STOP'Error with allocation of unused in `disaggregate`'
  endif
  
! p(i,j) is the precip read from the input file but it is modified here
! for subsequent time steps  
  
    
!  if(mod((jz-1)/deltaT2,deltaT2).eq.0)then
  if(mod(jz+23,deltaT2).eq.0)then
    do i=1,ycount                                                                                              
        do j=1,xcount                                                                                             
            no_hrs_precip(i,j)=min(int(p(i,j)+1),24)  
            ! assumes 1 mm/hr and 24 hour data time step 
            ! time step to be obtained in read_pcp_nc from the precip file header
            pLocal(i,j)=p(i,j)                                                                                   
            unused(i,j)=p(i,j)  
            nhdt=0
        end do                                                                                                 
    end do  
!    write(777,*)'xx',pLocal(ycount/2,xcount/2)
  endif
  
  deltat=24
  
    nhdt=nhdt+1
    do i=1,ycount                                                                                             
      do j=1,xcount                                                                                          
        if(pLocal(i,j).le.0.0)then                                                                              
          p(i,j)=0.0                   
        elseif(pLocal(i,j).le.smearfactor)then                                                                  
          if(unused(i,j).ge.0.0)then                                                                       
           p(i,j)=unused(i,j)                                                                       
           unused(i,j)=0.0                                                                                 
          endif                                                                                            
        elseif(pLocal(i,j).le.deltaT2)then                                                                       
          if(unused(i,j).gt.0.000001)then                                                                  
            p(i,j)=amin1(smearfactor,unused(i,j))                                                   
            unused(i,j)=unused(i,j)-p(i,j)                                                          
          else                                                                                             
            p(i,j)=0.0                                                                              
          endif                                                                                            
        elseif(pLocal(i,j).gt.deltaT2)then                                                                       
          p(i,j)=pLocal(i,j)/deltaT2                                                                      
        else                                                                                               
          print*,'nhdt,i,j,pLocal(i,j)/',nhdt,i,j,pLocal(i,j)                                                        
          print*,'This should never happen!!!'                                                             
          stop                                                                                             
        endif                                                                                              
      end do                                                                                               
    end do         
    
 !   write(777,*)jz,mod(jz+23,deltaT2),unused(ycount/2,xcount/2),p(ycount/2,xcount/2)

  firstpass=.false.
    
  end subroutine disaggregate