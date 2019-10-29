c     **********************************************************************
c
c      Copyright (c) 2014, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     calculate cell volume

      double precision function cellvol(cell)
c
      implicit none
c
      double precision cell(6)
      double precision d2r, c(6)
      parameter (d2r=0.01745329251994329576d0)
c
      c(1) = cell(1)
      c(2) = cell(2)
      c(3) = cell(3)
      c(4) = cell(4)*d2r
      c(5) = cell(5)*d2r
      c(6) = cell(6)*d2r
      cellvol = c(1)*c(2)*c(3) *
     +  sqrt ( (1.d0 - cos(c(4))**2 - cos(c(5))**2 - cos(c(6))**2) +
     +         2.d0 * (cos(c(4))    * cos(c(5))    * cos(c(6))))
c
      return
      end
