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

      subroutine setbmat(cell,wave,b)
c
      implicit none
c
      double precision cell(6),wave,b(3,3)
      double precision rcell(6)
      double precision d2r
      parameter (d2r=0.01745329251994329576d0)
c
      call reccell(cell,rcell)
c
      b(1,1) =  rcell(1)
      b(1,2) =  rcell(2) * cos(rcell(6)*d2r)
      b(1,3) =  rcell(3) * cos(rcell(5)*d2r)
      b(2,1) =  0.0d0
      b(2,2) =  rcell(2) * sin(rcell(6)*d2r)
      b(2,3) = -rcell(3) * sin(rcell(5)*d2r) * cos(cell(4)*d2r)
      b(3,1) = 0.0d0
      b(3,2) = 0.0d0
      b(3,3) = wave/cell(3)
c
      return
      end
