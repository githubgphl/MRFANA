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

c     handling of multi-record MTZ files (including batch information)

c
c     batch header limits
c
      integer    maxrbatch
      parameter (maxrbatch=185)
      integer    maxcbatch
      parameter (maxcbatch=94)
c
      real rbatch(maxrbatch)
      character*1 cbatch(maxcbatch)
      integer ibatch(29)
c
      equivalence (ibatch(1),rbatch(1))
c
      integer i_batch_info_cell
      integer i_batch_info_wave
      integer i_batch_info_phistart
      integer i_batch_info_jumpax
      integer i_batch_info_angle
      integer n_batch_info
c
      parameter (i_batch_info_cell=1)
      parameter (i_batch_info_wave=7)
      parameter (i_batch_info_phistart=8)
      parameter (i_batch_info_jumpax=10)
      parameter (i_batch_info_angle=11)

      parameter (n_batch_info=11)

