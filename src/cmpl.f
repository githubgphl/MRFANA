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

c     computing completeness

      subroutine cmpl(rlow,rhigh,ityp,nrefpos)
c
      implicit none
c
      include 'params.inc'
      include 'cmpl.inc'
      include 'units.inc'
      include 'info.inc'
      include 'aniso.inc'
c
      double precision rlow, rhigh
      integer nrefpos(6), ityp
c     ityp = 0 : low-resolution limit is excluded
c            1 :                         included
c
      integer i, hklmin(3), hklmax(3), hkl(3), ih, ik, il, isysabs, ic,
     +  ellval
      real epsi
      double precision reso
c
      integer inasu
      double precision dstar
      external inasu
c
      equivalence (ih,hkl(1)),(ik,hkl(2)),(il,hkl(3))
c
c----+=================================================================
c
c     check buffered values:
c
      do i = 1, cmpl_nshell
        if ((cmpl_resl(1,i).eq.rlow ).and.
     +      (cmpl_resl(2,i).eq.rhigh)     ) then
          nrefpos = 0
c         normal
          nrefpos(1) = cmpl_nref(1,i)
c         ano
          nrefpos(2) = cmpl_nref(2,i)
          if (l_ell) then
c           inside ellipsoid
            nrefpos(3) = cmpl_nref_ell(1,1,i)
            nrefpos(4) = cmpl_nref_ell(2,1,i)
c           outside ellipsoid
            nrefpos(5) = cmpl_nref_ell(1,2,i)
            nrefpos(6) = cmpl_nref_ell(2,2,i)
          end if
          return
        else if (cmpl_resl(1,i).lt.rhigh) then
          goto 5
        end if
      end do
 5    continue
c
c     the next part is taken more or less directly from $CCP4/src/unique.f
c
      do i = 1, 3
        hklmax(i) = int(cell(i)/(rhigh-0.01d0))+1
      end do
c
      call hklrange(hklmin(1),hklmax(1),
     +              hklmin(2),hklmax(2),
     +              hklmin(3),hklmax(3))
c
      nrefpos = 0
c
c     loop over reflections
c
      do ih = hklmin(1), hklmax(1)
        do ik = hklmin(2), hklmax(2)
          do il = hklmin(3), hklmax(3)
c
c           since we can't observe the (0,0,0) reflection:
            if ((hkl(1).eq.0).and.
     +          (hkl(2).eq.0).and.
     +          (hkl(3).eq.0)     ) goto 10
c
c           is reflection in asymmetric unit:
            if (inasu(hkl,nsyml) .le. 0) goto 10
c
c           calculate resolution
            reso = 1.0d0/sqrt(dstar(hkl,0))

            if (reso.lt.rhigh) goto 10
            if (ityp.eq.0) then
              if (reso.ge.rlow) goto 10
            else
              if (reso.gt.rlow) goto 10
            end if
c
c           using multiplicity routine to get systematic absent flag
            isysabs = 0
c
            call epslon(hkl,epsi,isysabs)
            if (isysabs.eq.1) goto 10
c
c           all checks passed: count it
            nrefpos(1) = nrefpos(1) + 1

            if (l_ell) then
              call ellchk(hkl,ellval)
              if (ellval.eq.1) then
                nrefpos(3) = nrefpos(3) + 1
              else
                nrefpos(5) = nrefpos(5) + 1
              end if
            end if
            call centr(hkl, ic)
c
c           anomalous differences
            if (ic.eq.0) then
              nrefpos(2) = nrefpos(2) + 1

              if (l_ell) then
                if (ellval.eq.1) then
                  nrefpos(4) = nrefpos(4) + 1
                else
                  nrefpos(6) = nrefpos(6) + 1
                end if
              end if

            end if
c
 10         continue
          end do
        end do
      end do
c
      if (cmpl_nshell.lt.cmpl_maxshell) then
        cmpl_nshell = cmpl_nshell + 1
        cmpl_resl(1,cmpl_nshell) = rlow
        cmpl_resl(2,cmpl_nshell) = rhigh
        cmpl_nref(1,cmpl_nshell) = nrefpos(1)
        cmpl_nref(2,cmpl_nshell) = nrefpos(2)
        cmpl_nref_ell(1,1,cmpl_nshell) = nrefpos(3)
        cmpl_nref_ell(2,1,cmpl_nshell) = nrefpos(4)
        cmpl_nref_ell(1,2,cmpl_nshell) = nrefpos(5)
        cmpl_nref_ell(2,2,cmpl_nshell) = nrefpos(6)
      endif
c
c----+=================================================================
c
      return
      end
