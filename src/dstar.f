c     **********************************************************************
c
c      Copyright (c) 2017, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     computing resolution using the correct (dataset) information

      double precision function dstar(hkl,idataset)
c
      implicit none
c
      include 'params.inc'
      include 'cell.inc'
c
      integer hkl(3), idataset
c
      double precision d(3)
      integer i, j
c
c----+=================================================================
c
      d = 0.0d0
      do i = 1, 3
        do j = 1, 3
          if (idataset.le.0) then
            d(j) = d(j) + hkl(i)*cell_orth_mat_inv(i,j)
          else
            d(j) = d(j) + hkl(i)*datcell_orth_mat_inv(i,j,idataset)
          end if
        end do
      end do
      dstar = 0.0d0
      do i = 1, 3
        dstar = dstar + d(i)**2
      end do
c
c----+=================================================================
c
      return
      end
