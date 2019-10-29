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

c     See: http://www.ccp4.ac.uk/html/mtzlib.html#orientation

      subroutine read_batch(filtyp,batch_info,batlim)
c
      implicit none
c
      include 'batch.inc'
      include 'units.inc'
      include 'xds.inc'
c
      character filtyp*13
      integer batlim(2)
c     0
c        NWORDS,NINTGR,NREALS,IORTYP,LBCELL(6),MISFLG,
c     11
c        JUMPAX,NCRYST,LCRFLG,LDTYPE,JSCAXS,NBSCAL,NGONAX,LBMFLG,
c     19
c        NDET,LBSETID,INTPAD(8),
c     29
c        CELL(6),UMATR(3,3),PHXYZ(3,2),CRYDAT(12),DATUM(3),
c     65
c        PHISTT,PHIEND,SCANAX(3),TIME1,TIME2,
c     72
c        BSCALE,BBFAC,SDSCAL,SDBB,PHIRANGE,BATPAD(11),
c     88
c        E1(3),E2(3),E3(3),
c     97
c        GONPAD(12),SOURCE(3),S0(3),BEMDAT(25),
c     140
c        DX1,THETA1,DETLM1(2,2),DX2,THETA2,DETLM2(2,2),DETPAD(33)
c
c     1-6 = cell dimensions
c     7   = wavelength [A]
c     8,9 = start & stop values of Phi (degrees) relative to datum
c     10 = jumpax
      double precision batch_info(n_batch_info,batlim(1):batlim(2))
      double precision angle
      double precision dbatch(maxrbatch)
c
      integer ibat,i, j, ibat_skip(2)
c
c----+=================================================================
c
      if (filtyp.eq.'MTZ') then

        ibat_skip(1)=0
        ibat_skip(2)=0

        call lrrewd(lun)
c
        do i = batlim(1), batlim(2)
          do j = 1, n_batch_info
            batch_info(j,i) = 0.0d0
          end do
        end do

        ibat = 0
        do while (ibat.ge.0)
          call lrbat(lun,ibat,rbatch,cbatch,0)
          if (ibat.lt.0) exit

          if ((ibat.lt.batlim(1)).or.(ibat.gt.batlim(2))) then
            if (ibat_skip(2).eq.0) then
              ibat_skip(1) = ibat
            end if
            ibat_skip(2)=ibat
            cycle
c            write(stdout,*) 'ibat,batlim',ibat,batlim
c            call error_exit('batch number outside range')
          else
            if (ibat_skip(2).ne.0) then
              write(stdout,'(/,1x,a,i8,a,i8)') 'skipping BATCHes '
     .          ,ibat_skip(1),'...',ibat_skip(2)
            end if
            ibat_skip(2)=0
          end if

          do i = 1, maxrbatch
            dbatch(i) = dble(rbatch(i))
          end do

          do i = 1,6
            batch_info((i_batch_info_cell+i-1),ibat) = dbatch((i+29))
          end do
          batch_info( i_batch_info_wave       ,ibat) = dbatch(116)
          batch_info( i_batch_info_phistart   ,ibat) = dbatch(66)
          batch_info((i_batch_info_phistart+1),ibat) = dbatch(67)
          batch_info( i_batch_info_jumpax     ,ibat) = float(ibatch(12))
c
c         avoid NaN:
          if ((dbatch(36).eq.dbatch(36)).and.(ibatch(12).gt.0)) then
            call jumpax_angle(
     +        batch_info(i_batch_info_cell,ibat),
     +        batch_info(i_batch_info_wave,ibat),
     +        dbatch(36),
     +        ibatch(12),
     +        ibatch(18),
     +        dbatch(89),
     +        dbatch(92),
     +        dbatch(95),
     +        dbatch(63),
     +        dbatch(68),
     +        angle)
          else
            angle = 0.0d0
          end if
          batch_info(i_batch_info_angle,ibat) = angle
        end do

        if (ibat_skip(2).ne.0) then
          write(stdout,'(/,1x,a,i8,a,i8,/)') 'skipping BATCHes '
     .      ,ibat_skip(1),'...',ibat_skip(2)
        end if
c
      else if (filtyp(1:4).eq.'XDS_') then
c
        do ibat = batlim(1), batlim(2)
          do j = 1, n_batch_info
            batch_info(j,ibat) = 0.0
          end do

          batch_info( i_batch_info_wave       ,ibat) = xds_wave
          do i = 1,6
            batch_info((i_batch_info_cell+i-1),ibat) = xds_cell(i)
          end do

        end do
c
      end if
c
c----+=================================================================
c
      return
      end
