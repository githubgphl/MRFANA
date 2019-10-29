c     **********************************************************************
c
c      Copyright (c) 2017, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     handling cell (real- and reciprocal cell)

      double precision cell_orth_mat(3,3)
      double precision cell_orth_mat_inv(3,3)
      double precision datcell_orth_mat(3,3,maxdatasets)
      double precision datcell_orth_mat_inv(3,3,maxdatasets)

      common /d_cell/ cell_orth_mat, cell_orth_mat_inv,
     +                datcell_orth_mat, datcell_orth_mat_inv
