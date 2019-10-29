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

c     sort on H K L M/ISYM

      subroutine sort_refl(refl,ncol,nref)
c
      implicit none
c
      include 'params.inc'
      include 'lookup.inc'
c
      integer ncol, nref
      double precision refl(1:ncol,1:nref)
c
      double precision, allocatable :: reft(:,:)
      real, allocatable :: hkl(:)
      integer, allocatable :: pt(:)
      integer istat, iref, jref, jcol, nH, mH, nK, mK, nL, mL
      double precision Hold, Kold, Lold
c
c----+=================================================================
c
      allocate(reft(1:ncol,1:nref),stat=istat)
      if (istat.ne.0) call error_exit('unable to allocate reft array')
      allocate(hkl(1:nref),stat=istat)
      if (istat.ne.0) call error_exit('unable to allocate hkl array')
      allocate(pt(1:nref),stat=istat)
      if (istat.ne.0) call error_exit('unable to allocate pt array')
c     ------------------------------------------------------
c     sort on H
c
      do iref = 1, nref
        hkl(iref) = sngl(refl(icol(i_h),iref))
      end do
      call quick_sort(nref,hkl,pt)
      do iref = 1, nref
        do jcol = 1, ncol
          reft(jcol,iref) = refl(jcol,pt(iref))
        end do
      end do
      do iref = 1, nref
        do jcol = 1, ncol
          refl(jcol,iref) = reft(jcol,iref)
        end do
      end do
c     ------------------------------------------------------
c     sort on K
c
      Hold = -999999.9d0
      nH = 0
      do iref = 1, nref
        if (refl(icol(i_h),iref).ne.Hold) then
          if (nH.gt.0) then
            mH = 0
            do jref = iref-nH, iref-1
              mH = mH + 1
              hkl(mH) = sngl(refl(icol(i_k),jref))
            end do
            call quick_sort(mH,hkl,pt)
            do jref = 1, mH
              do jcol = 1, ncol
                reft(jcol,jref) = refl(jcol,(pt(jref)+iref-nH-1))
              end do
            end do
            mH = 0
            do jref = iref-nH, iref-1
              mH = mH + 1
              do jcol = 1, ncol
                refl(jcol,jref) = reft(jcol,mH)
              end do
            end do
          end if
          nH = 0
        end if
        nH = nH + 1
        Hold = refl(icol(i_h),iref)
      end do
c
      iref = nref + 1
      if (nH.gt.0) then
        mH = 0
        do jref = iref-nH, iref-1
          mH = mH + 1
          hkl(mH) = sngl(refl(icol(i_k),jref))
        end do
        call quick_sort(mH,hkl,pt)
        do jref = 1, mH
          do jcol = 1, ncol
            reft(jcol,jref) = refl(jcol,(pt(jref)+iref-nH-1))
          end do
        end do
        mH = 0
        do jref = iref-nH, iref-1
          mH = mH + 1
          do jcol = 1, ncol
            refl(jcol,jref) = reft(jcol,mH)
          end do
        end do
      end if
c     ------------------------------------------------------
c     sort on L
c
      Hold = -999999.9d0
      Kold = -999999.9d0
      nK = 0
      do iref = 1, nref
        if ((refl(icol(i_h),iref).ne.Hold).or.
     +      (refl(icol(i_k),iref).ne.Kold)    ) then
          if (nK.gt.0) then
            mK = 0
            do jref = iref-nK, iref-1
              mK = mK + 1
              hkl(mK) = sngl(refl(icol(i_l),jref))
            end do
            call quick_sort(mK,hkl,pt)
            do jref = 1, mK
              do jcol = 1, ncol
                reft(jcol,jref) = refl(jcol,(pt(jref)+iref-nK-1))
              end do
            end do
            mK = 0
            do jref = iref-nK, iref-1
              mK = mK + 1
              do jcol = 1, ncol
                refl(jcol,jref) = reft(jcol,mK)
              end do
            end do
          end if
          nK = 0
        end if
        nK = nK + 1
        Hold = refl(icol(i_h),iref)
        Kold = refl(icol(i_k),iref)
      end do
c
      iref = nref + 1
      if (nK.gt.0) then
        mK = 0
        do jref = iref-nK, iref-1
          mK = mK + 1
          hkl(mK) = sngl(refl(icol(i_l),jref))
        end do
        call quick_sort(mK,hkl,pt)
        do jref = 1, mK
          do jcol = 1, ncol
            reft(jcol,jref) = refl(jcol,(pt(jref)+iref-nK-1))
          end do
        end do
        mK = 0
        do jref = iref-nK, iref-1
          mK = mK + 1
          do jcol = 1, ncol
            refl(jcol,jref) = reft(jcol,mK)
          end do
        end do
      end if
c     ------------------------------------------------------
c     sort on M/ISYM
c
      Hold = -999999.9d0
      Kold = -999999.9d0
      Lold = -999999.9d0
      nL = 0
      do iref = 1, nref
        if ((refl(icol(i_h),iref).ne.Hold).or.
     +      (refl(icol(i_k),iref).ne.Kold).or.
     +      (refl(icol(i_l),iref).ne.Lold)    ) then
          if (nL.gt.0) then
            mL = 0
            do jref = iref-nL, iref-1
              mL = mL + 1
              hkl(mL) = sngl(refl(icol(i_misym),jref))
            end do
            call quick_sort(mL,hkl,pt)
            do jref = 1, mL
              do jcol = 1, ncol
                reft(jcol,jref) = refl(jcol,(pt(jref)+iref-nL-1))
              end do
            end do
            mL = 0
            do jref = iref-nL, iref-1
              mL = mL + 1
              do jcol = 1, ncol
                refl(jcol,jref) = reft(jcol,mL)
              end do
            end do
          end if
          nL = 0
        end if
        nL = nL + 1
        Hold = refl(icol(i_h),iref)
        Kold = refl(icol(i_k),iref)
        Lold = refl(icol(i_l),iref)
      end do
c
      iref = nref + 1
      if (nL.gt.0) then
        mL = 0
        do jref = iref-nL, iref-1
          mL = mL + 1
          hkl(mL) = sngl(refl(icol(i_misym),jref))
        end do
        call quick_sort(mL,hkl,pt)
        do jref = 1, mL
          do jcol = 1, ncol
            reft(jcol,jref) = refl(jcol,(pt(jref)+iref-nL-1))
          end do
        end do
        mL = 0
        do jref = iref-nL, iref-1
          mL = mL + 1
          do jcol = 1, ncol
            refl(jcol,jref) = reft(jcol,mL)
          end do
        end do
      end if
c
      deallocate(reft,stat=istat)
      if (istat.ne.0) call error_exit('unable to deallocate reft array')
      deallocate(hkl,stat=istat)
      if (istat.ne.0) call error_exit('unable to deallocate hkl array')
      deallocate(pt,stat=istat)
      if (istat.ne.0) call error_exit('unable to deallocate pt array')
c
c----+=================================================================
c
      return
      end
