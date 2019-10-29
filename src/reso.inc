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

c     resolution information

      double precision reslow, reshigh
      double precision current_reslow, current_reshigh
      double precision BinningByEqualNumber(maxshell)
      double precision BinningByEqualNumber_run(maxshell)
      double precision reslow_ell, reshigh_ell
c
      common /d_reso/ reslow, reshigh, current_reslow, current_reshigh,
     +  BinningByEqualNumber, BinningByEqualNumber_run, reslow_ell,
     +  reshigh_ell
