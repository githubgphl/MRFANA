c     **********************************************************************
c
c      Copyright (c) 2010, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     run information (defined via RUN cards on stdin)

      double precision run_reslim(2,0:maxruns)
c
      common /runs_d/ run_reslim
c
      integer numrun,
     +  run_num     (           0:maxruns),
     +  run_i2n     (           0:maxruns),
     +  run_batlim  (         2,0:maxruns),
     +  run_batexcl (maxrunexcl,0:maxruns),
     +  run_nbatexcl(           0:maxruns),
     +  run_set     (           0:maxruns),
     +  run_nref    (           0:maxruns)
c
      common /runs_i/ numrun, run_num, run_i2n,
     +  run_batlim, run_batexcl, run_nbatexcl, run_set, run_nref
c
      character*65 run_dname(0:maxruns)
c
      common /runs_c/ run_dname
