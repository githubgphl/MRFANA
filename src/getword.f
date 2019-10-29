c     **********************************************************************
c
c      Copyright (c) 2006, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

      subroutine getword(lin,wrd)
c
      implicit none
c
      character*(*) lin,wrd
c
      integer i
c
c----+=================================================================
c
      wrd = ' '
      i = index(lin,' ')
      wrd = lin(1:(i-1))
c
c----+=================================================================
c
      return
      end
