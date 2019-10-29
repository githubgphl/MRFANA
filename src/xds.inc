c     **********************************************************************
c
c      Copyright (c) 2015, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     XDS specific

      double precision xds_cell(6), xds_wave
c
      common /xds_d/ xds_cell, xds_wave
c
      integer xds_item_iset
c
      common /xds_i/ xds_item_iset
