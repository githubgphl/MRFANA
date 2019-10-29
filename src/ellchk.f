c     **********************************************************************
c
c      Copyright (c) 2017, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein, Ian J. Tickle and Gerard Bricogne
c
c     **********************************************************************

c     check if a reflection is within or outside of defined ellipsoid

      subroutine ellchk(hkl,ival)
c
c i/p: hkl(3) = Miller indices
c      
c o/p: ival = 1 (inside) 0 (outside)
c
c
      implicit none
c
      include 'params.inc'
      include 'cell.inc'
      include 'aniso.inc'
c
      integer hkl(3)
      integer ival
c
      double precision r(3), d, rl(3)
      integer i, j
c
c----+=================================================================
c
c     default = outside
      ival = -1
c
c     get reciprocal lattice vector:
      rl = 0.0d0
      do i = 1, 3
        do j = 1, 3
          rl(j) = rl(j) + hkl(i)*cell_orth_mat_inv(i,j)
        end do
      end do
c
      r = 0.0d0
      do i = 1, 3
        do j = 1, 3
          r(j) = r(j) + EllFit_mat_inv(j,i)*rl(i)
        end do
      end do
      d = 0.0d0
      do i = 1, 3
        d = d + (r(i)/EllFit_res_inv(i))**2
      end do
      if (d.le.1.0d0) then
c       inside
        ival = 1
      end if
c
c----+=================================================================
c
      return
      end
