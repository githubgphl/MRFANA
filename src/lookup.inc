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

c     setup pointers into reflection array (easier to remember and to
c     read)

      integer icol(maxcol)
      integer i_h, i_k, i_l,
     +        i_i, i_sigi,
     +        i_xobs, i_yobs, i_zobs,
     +        i_lp,
     +        i_xcal, i_ycal, i_zcal,
     +        i_misym, i_batch, i_resol,
     +        i_reso3, i_flag, i_run, i_out,
     +        i_f, i_sigf, i_rand, i_frac, i_iset,
     +        i_ell
      parameter (i_h    = 1)
      parameter (i_k    = 2)
      parameter (i_l    = 3)
      parameter (i_i    = 4)
      parameter (i_sigi = 5)
      parameter (i_xobs = 6)
      parameter (i_yobs = 7)
      parameter (i_zobs = 8)
      parameter (i_lp   = 9)
      parameter (i_xcal =10)
      parameter (i_ycal =11)
      parameter (i_zcal =12)
      parameter (i_misym=13)
      parameter (i_batch=14)
      parameter (i_resol=15)
      parameter (i_reso3=16)
      parameter (i_flag =17)
      parameter (i_run  =18)
      parameter (i_out  =19)
c     if the file was passed through SCALA2SHARP:
      parameter (i_f    =20)
      parameter (i_sigf =21)
c     random number
      parameter (i_rand =22)
c     fraction observed
      parameter (i_frac =23)
c     pointer to SET identifier (needed for XSCALE data)
      parameter (i_iset =24)
c     pointer to ellipse in/out marker
      parameter (i_ell  =25)

c
      common /lookup_i/ icol
