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

c     see dtrek2scala_/mkordtr.f

      subroutine jumpax_angle(cell,wave,umat,jumpax,ngonax,axis1,axis2
     .  ,axis3,angles,rotation_axis,angle)
c
      implicit none
c
      double precision cell(6)
      double precision wave
      double precision umat(3,3)
      double precision angle
      double precision axis1(3)
      double precision axis2(3)
      double precision axis3(3)
      double precision angles(3)
      double precision rotation_axis(3)
      integer jumpax, ngonax
c
      double precision bmat(3,3)
      double precision ub(3,3)
      double precision dub(3,3)
      double precision d(3,3)
      double precision rmat(3,3,3)
      double precision v1(3)
      double precision veclen
      double precision a(3)
      double precision ax(3)
      double precision r2d
      parameter (r2d = 57.2957795130823208768d0)
      integer i, j, k
c
      double precision dotvec
      external dotvec
c
c----+=================================================================
c
cv
      k=jumpax
cv
      do j = 1, 3
        do i = 1, 3
          d(i,j) = 0.
          rmat(i,j,1) = 0.
          rmat(i,j,2) = 0.
          rmat(i,j,3) = 0.
        end do
        d(j,j) = 1.
        rmat(j,j,1) = 1.
        rmat(j,j,2) = 1.
        rmat(j,j,3) = 1.
      end do

      call setbmat(cell,wave,bmat)

c     UB = UMAT x BMAT
      ub = matmul(umat,bmat)

      if (ngonax.gt.0) then

        call rotax(axis1,angles(1),rmat(1,1,1))
        if (ngonax.gt.1) then
          call rotax(axis2,angles(2),rmat(1,1,2))
          if (ngonax.gt.2) then
            call rotax(axis3,angles(3),rmat(1,1,3))
          end if
        end if
c
        if (ngonax.eq.1) then
          do j = 1, 3
            do i = 1, 3
              d(i,j) = rmat(i,j,1)
            end do
          end do
        else
          d = matmul(matmul(rmat(:,:,1),rmat(:,:,2)),rmat(:,:,3))
c          call ml3mat(3,rmat(1,1,1),3,rmat(1,1,2),3,rmat(1,1,3),d)
        end if
c
      end if

c     DUB = D x UB
      dub = matmul(d,ub)
c
      do i = 1, 3
        ax(1) = 0.
        ax(2) = 0.
        ax(3) = 0.
        ax(i) = 1.
c
        v1 = matmul(dub,ax)
        veclen = sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
        v1(1) = v1(1)/veclen
        v1(2) = v1(2)/veclen
        v1(3) = v1(3)/veclen
        a(i) = r2d*acos(abs(dotvec(rotation_axis,v1)))

        if (i.eq.1) then
          angle = a(i)
          jumpax = i
        endif
        if (a(i).lt.angle) then
          angle = a(i)
          jumpax = i
        endif

      end do
c
c     should we return the smallest angle or the angle for a given
c     jumpax value? Stick with the latter:
      if (jumpax.ne.k) then
        angle = a(k)
      end if
c
c----+=================================================================
c
      return
      end
