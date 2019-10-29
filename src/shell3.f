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

c     wrappers to get resolution (in A) depending on storage/binning
c     method (with some checks and safety buffer).

      double precision function stol2res(stol)
c
      implicit none
c
      include 'params.inc'
      include 'info.inc'
      include 'reso.inc'
c
      double precision stol
c
c----+=================================================================
c
        if (lshell3) then
          if (stol.gt.0.0) then
            stol2res = 1.d0/(stol**(1./3.))
          else
            stol2res = current_reshigh - 0.001
          end if
        else
          if (stol.gt.0.0) then
            stol2res = 1.d0/sqrt(stol)
          else
            stol2res = current_reshigh - 0.001
          end if
        end if
c
c----+=================================================================
c
      return
      end
c
c----+=================================================================
c
      double precision function res2stol(res)
c
      implicit none
c
      include 'params.inc'
      include 'info.inc'
c
      double precision res
c
c----+=================================================================
c
        if (lshell3) then
          if (res.gt.0.0d0) then
            res2stol = 1./(res**3)
          else
            res2stol = 0.999999d0
          end if
        else
          if (res.gt.0.0d0) then
            res2stol = 1./(res**2)
          else
            res2stol = 0.999999d0
          end if
        end if
c
c----+=================================================================
c
      return
      end
