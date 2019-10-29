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

c     simple wrapper around file closing function (depends o file
c     format)

      subroutine close_file(filtyp)
c
      implicit none
c
      include 'units.inc'
c
      character filtyp*13
c
c----+=================================================================
c
      if (filtyp.eq.'MTZ') then
        call lrclos(lun)
      else
        close(lun)
      end if
c
c----+=================================================================
c
      return
      end
