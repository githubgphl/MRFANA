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

c     various defaults and user-supplied settings

c     strings
      character pgmnam*1024, filnam*1024
c
      common /args_c/ pgmnam, filnam

c     integer
      integer iverb, nshell, nrefbin, nrefbin_run, nIce,
     +  rand_seed,
     +  iset_nshell, iwide_format,
     +  calc_type_r(2,3),
     +  N_BinningByEqualNumber, N_BinningByEqualNumber_run
c
      common /args_i/ iverb, nshell, nrefbin, nrefbin_run, nIce,
     +  rand_seed,
     +  iset_nshell, iwide_format,
     +  calc_type_r,
     +  N_BinningByEqualNumber, N_BinningByEqualNumber_run

c     double precision
      double precision check_part_min
      double precision check_part_max
c
      common /args_d/ check_part_min, check_part_max

c     logical
      logical         l_skip_empty_runs,
     +  L_BinningByEqualNumber, L_BinningByEqualNumber_run
c
      common /args_l/ l_skip_empty_runs,
     +  L_BinningByEqualNumber, L_BinningByEqualNumber_run

