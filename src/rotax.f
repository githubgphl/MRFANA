c     **********************************************************************
c
c      Copyright (c) 2014, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c calculates rotation matrix of angle ROT around axis AX (right-hand
c rule)

      subroutine rotax(ax,rot,mat)
c
      implicit none
c
      double precision ax(3), rot, mat(3,3)
c
      double precision c_rot, s_rot, c
c
      double precision d2r
      parameter      ( d2r     =  0.01745329251994329576d0 )

c
c----+=================================================================
c
      c_rot = cos(rot*d2r)
      s_rot = sin(rot*d2r)
      c     = 1.0d0 - c_rot

      mat(1,1) = 1.0d0 + c*(ax(1)*ax(1)-1)
      mat(1,2) = -ax(3)*s_rot+c*ax(1)*ax(2)
      mat(1,3) =  ax(2)*s_rot+c*ax(1)*ax(3)

      mat(2,1) =  ax(3)*s_rot+c*ax(1)*ax(2)
      mat(2,2) = 1.0d0 + c*(ax(2)*ax(2)-1)
      mat(2,3) = -ax(1)*s_rot+c*ax(2)*ax(3)

      mat(3,1) = -ax(2)*s_rot+c*ax(1)*ax(3)
      mat(3,2) =  ax(1)*s_rot+c*ax(2)*ax(3)
      mat(3,3) = 1.0d0 + c*(ax(3)*ax(3)-1)
c
c----+=================================================================
c
      return
      end
