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

      integer symm, nsyml, nshell_run(maxruns), batlim(2)
c
      common /info_i/ symm, nsyml, nshell_run, batlim
c
      character*10 spgrn
c
      common /info_c/ spgrn
c
      double precision cell(6),
     +     resmin , resmax , resstp ,
     +     resmin3, resmax3, resstp3,
     +     resmin_run(maxruns) ,
     +     resmax_run(maxruns) ,
     +     resstp_run(maxruns) ,
     +     resmin3_run(maxruns),
     +     resmax3_run(maxruns),
     +     resstp3_run(maxruns)
c
      common /info_d/ cell,
     +                resmin , resmax , resstp ,
     +                resmin3, resmax3, resstp3,
     +                resmin_run , resmax_run , resstp_run ,
     +                resmin3_run, resmax3_run, resstp3_run
c
      logical lshell3
c
      common /info_l/ lshell3
