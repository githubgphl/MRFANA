c     **********************************************************************
c
c      Copyright (c) 2011, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

      integer stderr, stdin, stdout, lun, pltlun
c
      parameter (stderr=0, stdin=5, stdout=6, pltlun=99)
c
      common /units_i/ lun
