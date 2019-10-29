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

      double precision function dotvec(a,b)
c
      implicit none
c
      double precision a(3), b(3)
c
      dotvec = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
c
      return
      end
