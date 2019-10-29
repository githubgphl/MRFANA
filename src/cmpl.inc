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

c     handling completeness (in a hopefully efficient way within MRFANA)

      integer cmpl_maxshell
      parameter (cmpl_maxshell=maxshell+maxruns*2)
c
      integer cmpl_nshell, cmpl_nref(2,cmpl_maxshell),
     +  cmpl_nref_ell(2,2,cmpl_maxshell)
c
      common /cmpl_i/ cmpl_nshell, cmpl_nref, cmpl_nref_ell
c
      double precision cmpl_resl(2,cmpl_maxshell)
c
      common /cmpl_d/ cmpl_resl
