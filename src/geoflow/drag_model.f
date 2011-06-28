C*******************************************************************
C* Copyright (C) 2003 University at Buffalo
C*
C* This software can be redistributed free of charge.  See COPYING
C* file in the top distribution directory for more details.
C*
C* This software is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
C*
C* Author: 
C* Description: 
C*
C*******************************************************************
C* $Id:$
C*

C***********************************************************************
      subroutine calc_drag_force (uvec, vsolid, vfluid, den_solid,
     &             den_fluid, vterminal, drag, tiny)
C***********************************************************************

      implicit none
      integer i, j
      double precision uvec(6), vsolid(2), vfluid(2)
      double precision den_fluid, den_solid
      double precision drag(4), vterminal, tiny

!     local variables and constants
      double precision exponant, temp, volf, denfrac
      double precision delv(2)
      parameter (exponant = 3.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     model in pitman-le paper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      volf = uvec(2)/uvec(1)
      temp = uvec(2)*(1.-volf)**(1.-exponant)/vterminal
      denfrac = den_fluid/den_solid

       
      do 10 i=1,4
10      drag(i)=0.
 
      ! fluid vel - solid vel
      if ( uvec(2).gt.tiny ) then
        do 20 j=1,2
20        delv(j) = vfluid(j)-vsolid(j)

!       compute individual drag-forces
        drag(1) = (1.-denfrac)*temp*delv(1)
        drag(2) = (1.-denfrac)*temp*delv(2)
        drag(3) = (1.-denfrac)*temp*delv(1)/denfrac
        drag(4) = (1.-denfrac)*temp*delv(2)/denfrac
      endif

      end subroutine calc_drag_force

