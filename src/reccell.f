c     **********************************************************************
c
c      Copyright (c) 2007, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     reciprocal cell handling

      subroutine reccell(cell,rcell)
c
      implicit none
c
      double precision cell(6),rcell(6)
      double precision d2r
      parameter (d2r=0.01745329251994329576d0)
      double precision c(6),v,cosa,cosb,cosg,sina,sinb,sing
      double precision cellvol
      external cellvol
c
      c = cell
      c(4) = d2r*cell(4)
      c(5) = d2r*cell(5)
      c(6) = d2r*cell(6)
c
      v = cellvol(cell)
c
      rcell(1) = c(2)*c(3)*sin(c(4))/v
      rcell(2) = c(3)*c(1)*sin(c(5))/v
      rcell(3) = c(1)*c(2)*sin(c(6))/v
c
      cosa = (cos(c(6))*cos(c(5))-cos(c(4))) / (sin(c(5))*sin(c(6)))
      sina = sqrt(1.d0-cosa**2)
      cosb = (cos(c(4))*cos(c(6))-cos(c(5))) / (sin(c(4))*sin(c(6)))
      sinb = sqrt(1.d0-cosb**2)
      cosg = (cos(c(4))*cos(c(5))-cos(c(6))) / (sin(c(4))*sin(c(5)))
      sing = sqrt(1.d0-cosg**2)
c
      rcell(4) = atan2(sina,cosa)/d2r
      rcell(5) = atan2(sinb,cosb)/d2r
      rcell(6) = atan2(sing,cosg)/d2r
c
      return
      end
