      subroutine errormsg(io_err,SYS_ERR,STAT,UNIT,cond)
      
!***********************************************************************
!    Copyright (C) 1987-2018 by Nicholas Kouwen  
         
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
      

      integer         io_err,sys_err,stat,unit,cond


      write(98,98001)io_err,sys_err,cond,unit,stat

98001 format(' io_err = ',i2,' sys_err = ',i2,' cond =',i2,

     *            ' unit = ',i2,' stat =',i2)


      end subroutine errormsg
