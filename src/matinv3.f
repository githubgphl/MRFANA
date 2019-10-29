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

c     See: http://fortranwiki.org/fortran/show/Matrix+inversion

      subroutine matinv3(a,b)
c
      implicit none
c
      double precision a(3,3), b(3,3)
c
      double precision detinv
c
c----+=================================================================
c

c     Calculate the inverse determinant of the matrix
      detinv = 1.0d0/
     +           (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
     +          - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)
     +          + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

c     Calculate the inverse of the matrix
      B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
      B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
      B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
      B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
      B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
      B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
      B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
      B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
c
c----+=================================================================
c
      return
      end
