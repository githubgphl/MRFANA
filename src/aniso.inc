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

c     Handling description of ellipsoid for computation of ellipsoidal
c     completeness. See STARANISO [1] and autoPROC [2,3] for details.
c
c     [1] Tickle, I.J., Flensburg, C., Keller, P., Paciorek, W., Sharff,
c         A., Vonrhein, C., Bricogne, G. (2018). STARANISO. Cambridge,
c         United Kingdom: Global Phasing Ltd.
c         http://staraniso.globalphasing.org/
c
c     [2] Vonrhein, C., Flensburg, C., Keller, P., Sharff, A., Smart,
c         O., Paciorek, W., Womack, T. & Bricogne, G. (2011). Data
c         processing and analysis with the autoPROC toolbox. Acta
c         Cryst. D67, 293-302.
c         https://www.globalphasing.com/autoproc/

      double precision EllFit_mat(3,3)
      double precision EllFit_mat_inv(3,3)
      double precision EllFit_res(3)
      double precision EllFit_res_inv(3)

      common /d_aniso/ EllFit_mat,
     +                 EllFit_mat_inv,
     +                 EllFit_res,
     +                 EllFit_res_inv

      logical l_ell

      common /l_aniso/ l_ell
