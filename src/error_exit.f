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

      subroutine error_exit(txt)
c
      implicit none
c
      include 'units.inc'
c
      character*(*) txt
      logical isatty
c
c----+=================================================================
c
      write(stdout,'(/,2a,/)') ' ERROR: ',txt(1:len_trim(txt))
      call flush(stdout)
      if ((.not.isatty(stdout)).or.(.not.isatty(stderr)))
     +  write(stderr,'(/,2a,/)') ' ERROR: ',txt(1:len_trim(txt))
      call flush(stderr)
      stop 1
c
c----+=================================================================
c
      end
