c     **********************************************************************
c
c      Copyright (c) 2007, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein, Claus Flensburg and Gerard Bricogne
c
c     **********************************************************************

c     read reflection data

      subroutine read_refl(filtyp,refl,ncol,ncoltot,ranges,nref,clabs
     +  ,ilabs)
c
      implicit none
c
      include 'args.inc'
      include 'units.inc'
      include 'params.inc'
      include 'lookup.inc'
      include 'info.inc'
      include 'runs.inc'
      include 'reso.inc'
      include 'xds.inc'
      include 'dset.inc'
      include 'xml.inc'
      include 'aniso.inc'
c
      character filtyp*13
      character*30 clabs(maxcol)
      integer nref(2), ncol, ncoltot, ilabs(maxcol)
      double precision refl(1:ncoltot,1:nref(1))
      double precision ranges(2,maxcol)
c
      character*1024 line
      integer iref, i, ihkl(3), jhkl(3), isym, ibat, irun, iexcl,
     +  iset, ibatsetid, nref_5, jrun, run2run(maxruns),
     .  numrun_off, nreset(maxruns)
      real ohkl(3)
      double precision ohkl_dp(3)
      real adata(maxcol), epsi
      double precision reso2, reso3, reso, rand
      logical leof, lsort, luse, lrun2run(maxruns)
c
      integer istat, isysabs
      integer, allocatable :: bat2dsetid(:)
      real, parameter :: huge_real = HUGE(1.0)
      real, parameter :: tiny_real = TINY(1.0)
      real, parameter :: sqrt_huge_real = sqrt(HUGE(1.0))
      double precision, parameter :: huge_dp = HUGE(1.0d0)
      double precision, parameter :: tiny_dp = TINY(1.0d0)
      double precision, parameter :: sqrt_huge_dp = sqrt(HUGE(1.0d0))
c
      double precision dstar
c      real     sthlsq
c      external sthlsq
      double precision calc_resol
      external         calc_resol
c
c----+=================================================================
c
      lsort = .false.
      batlim(1) =  9999999
      batlim(2) = -9999999
c
      nref_5 = nint(0.05*nref(1))
c
c ----------------------------------------------------------------------
c
      if (filtyp.eq.'MTZ') then

        if (icol(i_batch).gt.0) then
          allocate(bat2dsetid(nint(ranges(1,icol(i_batch))):
     +                        nint(ranges(2,icol(i_batch)))),stat=istat)
          if (istat.ne.0)
     +      call error_exit('unable to allocate bat2dsetid array')
        else
          allocate(bat2dsetid(1:1),stat=istat)
          if (istat.ne.0)
     +      call error_exit('unable to allocate bat2dsetid(1:1) array')
        end if
        bat2dsetid = 0

        call lrrewd(lun)
        leof = .false.
        nref(2) = 0
        do while (.not.leof)
          nref(2) = nref(2) + 1
          if (nref(2).le.nref(1)) then
            if (mod(nref(2),nref_5).eq.0) then
              write(stdout,'(a,i10,a,f5.1,a)')
     +          '    reading reflection ',nref(2),
     +          ' ( ',(100.0*nref(2))/nref(1),' %)'
            end if
            call lrrefl(lun,reso2,adata,leof)
c
c           use same routine for resolution calculation throughout
            ihkl(1) = nint(adata(1))
            ihkl(2) = nint(adata(2))
            ihkl(3) = nint(adata(3))
c
            refl(icol(i_ell),nref(2)) = 0.0d0
            if (l_ell) then
              call ellchk(ihkl,i)
              refl(icol(i_ell),nref(2)) = dble(i)
            end if
c
            refl(icol(i_run),nref(2)) = -1.0d0

            if (icol(i_batch).gt.0) then
              ibat = nint(adata(icol(i_batch)))
              if (bat2dsetid(ibat).eq.0) then
                call lrbsetid(lun,ibat,ibatsetid)
                bat2dsetid(ibat) = ibatsetid
              end if
            else
              ibat = 1
              bat2dsetid(ibat) = 1
            end if
            reso2 = dstar(ihkl,bat2dsetid(ibat))

            iset = 0
            reso = 1./sqrt(reso2)
c
c           check overall resolution limit
c
            if ((reso.gt.(reslow +reseps)).or.
     +          (reso.lt.(reshigh-reseps))     ) then
              nref(2) = nref(2) - 1
              goto 20
            end if
c
c           check ellipsoidal resolution limit
c
            if (l_ell.and.(refl(icol(i_ell),nref(2)).gt.0.0d0)) then
              reslow_ell  = max(reslow_ell ,reso)
              reshigh_ell = min(reshigh_ell,reso)
            end if
c
c           check systematic absences
c
            call epslon(ihkl,epsi,isysabs)
            if (isysabs.gt.0) then
              nref(2) = nref(2) - 1
              goto 20
            end if
c
c           check if this reflection is part of a RUN
c
            luse = .true.
            if (numrun.gt.0) then
              luse = .false.
              do irun = 1, numrun
                if ((ibat.ge.run_batlim(1,irun)   ).and.
     +              (ibat.le.run_batlim(2,irun)   ).and.
     +              ( (iset.eq.run_set(irun)).or.
     +                (iset.eq.0)                )      ) then
c
c                 do we exclude this batch?
c
                  do iexcl = 1, run_nbatexcl(irun)
                    if (ibat.eq.run_batexcl(iexcl,irun)) then
                      goto 10
                    end if
                  end do
c
c                 is it inside resolution limits?
c
                  if ((reso.le.(run_reslim(1,irun)+reseps)).and.
     +                (reso.ge.(run_reslim(2,irun)-reseps))     ) then
                    luse = .true.
                    refl(icol(i_run),nref(2)) = dble(irun)
                    run_nref(irun) = run_nref(irun) + 1
                    goto 10
                  end if
                end if
              end do
 10           continue
            end if
            if (.not.luse) then
              nref(2) = nref(2) - 1
              goto 20
            end if
            run_reslim(1,0) = max(run_reslim(1,0),reso)
            run_reslim(2,0) = min(run_reslim(2,0),reso)
c
            batlim(1) = min(batlim(1),ibat)
            batlim(2) = max(batlim(2),ibat)
c
c           check RUN resolution limit
c
            reso3 = sqrt(reso2)**3
c
c     check if we need to sort
c
            if ((.not.lsort).and.(nref(2).gt.1)) then
              if ( (adata(1).lt.ohkl(1))      .or.
     +            ((adata(1).eq.ohkl(1)).and.
     +             (adata(2).lt.ohkl(2))     ).or.
     +            ((adata(1).eq.ohkl(1)).and.
     +             (adata(2).eq.ohkl(2)).and.
     +             (adata(3).lt.ohkl(3))     )    )
     +          then
                lsort = .true.
                if (iverb.gt.1) then
                  write(stdout,'(a,3i5,a,3i5,a,i10)')
     +              ' We will need re-sorting of the reflexions'
     +              //' : old HKL = ',nint(ohkl(1)),nint(ohkl(2))
     +              ,nint(ohkl(3)),' current HKL = ',nint(adata(1))
     +              ,nint(adata(2)),nint(adata(3)),' at record ',nref(2)
                end if
              end if
            end if
            if (.not.lsort) then
              ohkl(1) = adata(1)
              ohkl(2) = adata(2)
              ohkl(3) = adata(3)
              ohkl_dp = dble(ohkl)
            end if
c
            do i = 1, ncol
              refl(i,nref(2)) = dble(adata(i))
            end do

c           test for some ridiculous SigI
            if (refl(icol(i_sigi),nref(2)).lt.0.0d0) then
              write(stdout,'(a,3i5,a,i6)')
     +          ' Something odd with reflection HKL=(',ihkl,') batch =
     +          ',int(refl(i_batch,nref(2)))
              call error_exit
     +          ('sigI value below zero?')
            end if
            if (refl(icol(i_sigi),nref(2)).gt.huge_dp) then
              write(stdout,'(a,3i5,a,i6)')
     +          ' Something odd with reflection HKL=(',ihkl,') batch =
     +          ',int(refl(i_batch,nref(2)))
              call error_exit('sigI larger than allowed as a REAL*8?')
            end if
            if (refl(icol(i_sigi),nref(2)).gt.sqrt_huge_dp) then
              write(stdout,'(a,3i5,a,i6)')
     +          ' Something odd with reflection HKL=(',ihkl,') batch =
     +          ',int(refl(i_batch,nref(2)))
              call error_exit
     +          ('sigI^2 larger than allowed as a REAL*8 - this'/
     +          /' is not going to work!')
            end if

            refl(icol(i_resol),nref(2)) = reso2
            refl(icol(i_reso3),nref(2)) = reso3
            if (icol(i_frac).gt.0) then
              refl(icol(i_frac ),nref(2)) = dble(adata(icol(i_frac)))
            end if

c     update ranges
            if (nref(2).eq.1) then
              ranges(1,icol(i_resol)) = reso2
              ranges(2,icol(i_resol)) = reso2
              ranges(1,icol(i_reso3)) = reso3
              ranges(2,icol(i_reso3)) = reso3
              if (icol(i_frac).gt.0) then
                ranges(1,icol(i_frac )) = refl(icol(i_frac ),nref(2))
                ranges(2,icol(i_frac )) = refl(icol(i_frac ),nref(2))
              end if
            else
              ranges(1,icol(i_resol)) =
     +          min(reso2,ranges(1,icol(i_resol)))
              ranges(2,icol(i_resol)) =
     +          max(reso2,ranges(2,icol(i_resol)))
              ranges(1,icol(i_reso3)) =
     +          min(reso3,ranges(1,icol(i_reso3)))
              ranges(2,icol(i_reso3)) =
     +          max(reso3,ranges(2,icol(i_reso3)))
              if (icol(i_frac).gt.0) then
                ranges(1,icol(i_frac )) = min(ranges(1,icol(i_frac )),
     +            refl(icol(i_frac ),nref(2)))
                ranges(2,icol(i_frac )) = max(ranges(2,icol(i_frac )),
     +            refl(icol(i_frac ),nref(2)))
              end if
            end if
            call random_number(rand)
            refl(icol(i_rand),nref(2)) = rand
          else
            nref(2) = nref(1)
            goto 30
          end if
 20       continue
        end do
 30     call lrrewd(lun)
c
        deallocate(bat2dsetid)
c
c ----------------------------------------------------------------------
c
      else if (filtyp(1:4).eq.'XDS_') then
        rewind(lun)
        nref(2) = 0
 100    read(lun,'(a)',end=110) line
        if (line(1:1).ne.'!') then
          nref(2) = nref(2) + 1
          if (mod(nref(2),nref_5).eq.0) then
            write(stdout,'(a,i10,a,f5.1,a)')
     +        '    reading reflection ',nref(2),
     +        ' ( ',(100.0*nref(2))/nref(1),' %)'
          end if
          read(line,*) (refl(ilabs(i),nref(2)),i=1,ncol)
c
          if ((xds_item_iset.gt.0).and.(xds_item_iset.le.ncol)) then
            iset = nint(refl(ilabs(xds_item_iset),nref(2)))
          else
            iset = 0
          end if
c
c
c         generate some items
c
          ihkl(1) = nint(refl(icol(i_h),nref(2)))
          ihkl(2) = nint(refl(icol(i_k),nref(2)))
          ihkl(3) = nint(refl(icol(i_l),nref(2)))
c
          if (l_ell) then
            call ellchk(ihkl,i)
            refl(icol(i_ell),nref(2)) = dble(i)
          end if
c
c         since we have actual measurement indices: need to put back
c         into asu and get M/ISYM:
c
          call asuput(ihkl,jhkl,isym)
          refl(icol(i_misym),nref(2)) = dble(isym)
          refl(icol(i_h),nref(2)) = dble(jhkl(1))
          refl(icol(i_k),nref(2)) = dble(jhkl(2))
          refl(icol(i_l),nref(2)) = dble(jhkl(3))
c
c         take the Z column as BATCH
          ibat = int(refl(icol(i_zcal),nref(2)))+1
          refl(icol(i_batch),nref(2)) = dble(ibat)
          refl(icol(i_run),nref(2)) = -1.0d0
c          reso2 = 4.*sthlsq(jhkl(1),jhkl(2),jhkl(3))
c          reso = 1./sqrt(reso2)
          reso2 = dstar(jhkl,0)
          reso = 1.0d0/sqrt(reso2)
c
c         check overall resolution limit
c
          if ((reso.gt.(reslow +reseps)).or.
     +        (reso.lt.(reshigh-reseps))     ) then
            nref(2) = nref(2) - 1
            goto 100
          end if
c
c         check ellipsoidal resolution limit
c
          if (l_ell.and.(refl(icol(i_ell),nref(2)).gt.0.0d0)) then
            reslow_ell = min(reslow_ell,reso)
            reshigh_ell = max(reshigh_ell,reso)
          end if
c
c         check systematic absences
c
          call epslon(jhkl,epsi,isysabs)
          if (isysabs.gt.0) then
            nref(2) = nref(2) - 1
            goto 100
          end if
c
          if (icol(i_frac).gt.0) then
            refl(icol(i_frac),nref(2)) = refl(icol(i_frac),nref(2))
     .        /100.0d0
          end if
c
c         check if this reflection is part of a RUN
c
          luse = .true.
          if (numrun.gt.0) then
            luse = .false.
            do irun = 1, numrun

              if ((ibat.ge.run_batlim(1,irun)).and.
     +            (ibat.le.run_batlim(2,irun)).and.
     +            ( (iset.eq.run_set(irun)).or.
     +              (run_set(irun).eq.0   )    )   ) then
c
c               do we exclude this batch?
c
                do iexcl = 1, run_nbatexcl(irun)
                  if (ibat.eq.run_batexcl(iexcl,irun)) then
                    goto 130
                  end if
                end do
c
c               is it inside resolution limits?
c
                if ((reso.le.(run_reslim(1,irun)+reseps)).or.
     +              (reso.ge.(run_reslim(2,irun)-reseps))     ) then
                  luse = .true.
                  refl(icol(i_run),nref(2)) = dble(irun)
                  run_nref(irun) = run_nref(irun) + 1
                  goto 130
                end if
              end if
            end do
 130        continue
          end if
          if (.not.luse) then
            nref(2) = nref(2) - 1
            goto 100
          end if
          run_reslim(1,0) = max(run_reslim(1,0),reso)
          run_reslim(2,0) = min(run_reslim(2,0),reso)
c
          reso3 = sqrt(reso2)**3
c         
          batlim(1) = min(batlim(1),ibat)
          batlim(2) = max(batlim(2),ibat)
          refl(icol(i_resol),nref(2)) = reso2
          refl(icol(i_reso3),nref(2)) = reso3
          call random_number(rand)
          refl(icol(i_rand),nref(2)) = rand
c
c     check if we need to sort
c
          if ((.not.lsort).and.(nref(2).gt.1)) then
            if ( (refl(icol(i_h),nref(2)).lt.ohkl_dp(1))      .or.
     +          ((refl(icol(i_h),nref(2)).eq.ohkl_dp(1)).and.
     +           (refl(icol(i_k),nref(2)).lt.ohkl_dp(2))     ).or.
     +          ((refl(icol(i_h),nref(2)).eq.ohkl_dp(1)).and.
     +           (refl(icol(i_k),nref(2)).eq.ohkl_dp(2)).and.
     +           (refl(icol(i_l),nref(2)).lt.ohkl_dp(3))     )    )
     +        then
              lsort = .true.
              if (iverb.gt.1) then
                write(stdout,'(a,3i5,a,3i5,a,i10)')
     +            ' We will need re-sorting of the reflexions'
     +            //' : old HKL = ',
     +            nint(ohkl(1)),
     +            nint(ohkl(2)),
     +            nint(ohkl(3)),
     +            ' current HKL = ',
     +            nint(refl(icol(i_h),nref(2))),
     +            nint(refl(icol(i_k),nref(2))),
     +            nint(refl(icol(i_l),nref(2))),
     +            ' at record ',nref(2)
              end if
            end if
          end if
          if (.not.lsort) then
            ohkl_dp(1) = refl(icol(i_h),nref(2))
            ohkl_dp(2) = refl(icol(i_k),nref(2))
            ohkl_dp(3) = refl(icol(i_l),nref(2))
            ohkl = sngl(ohkl_dp)
          end if
c     update ranges
          if (nref(2).eq.1) then
            ranges(1,icol(i_misym)) = refl(icol(i_misym),nref(2))
            ranges(2,icol(i_misym)) = refl(icol(i_misym),nref(2))
            ranges(1,icol(i_h))     = refl(icol(i_h),nref(2))
            ranges(2,icol(i_h))     = refl(icol(i_h),nref(2))
            ranges(1,icol(i_k))     = refl(icol(i_k),nref(2))
            ranges(2,icol(i_k))     = refl(icol(i_k),nref(2))
            ranges(1,icol(i_l))     = refl(icol(i_l),nref(2))
            ranges(2,icol(i_l))     = refl(icol(i_l),nref(2))
            ranges(1,icol(i_batch)) = refl(icol(i_batch),nref(2))
            ranges(2,icol(i_batch)) = refl(icol(i_batch),nref(2))
            ranges(1,icol(i_resol)) = refl(icol(i_resol),nref(2))
            ranges(2,icol(i_resol)) = refl(icol(i_resol),nref(2))
            ranges(1,icol(i_reso3)) = refl(icol(i_reso3),nref(2))
            ranges(2,icol(i_reso3)) = refl(icol(i_reso3),nref(2))
          else
            ranges(1,icol(i_misym)) =
     +        min(refl(icol(i_misym),nref(2)),
     +        ranges(1,icol(i_misym)))
            ranges(2,icol(i_misym)) =
     +        max(refl(icol(i_misym),nref(2)),
     +        ranges(2,icol(i_misym)))
            ranges(1,icol(i_h)) =
     +        min(refl(icol(i_h),nref(2)),
     +        ranges(1,icol(i_h)))
            ranges(2,icol(i_h)) =
     +        max(refl(icol(i_h),nref(2)),
     +        ranges(2,icol(i_h)))
            ranges(1,icol(i_k)) =
     +        min(refl(icol(i_k),nref(2)),
     +        ranges(1,icol(i_k)))
            ranges(2,icol(i_k)) =
     +        max(refl(icol(i_k),nref(2)),
     +        ranges(2,icol(i_k)))
            ranges(1,icol(i_l)) =
     +        min(refl(icol(i_l),nref(2)),
     +        ranges(1,icol(i_l)))
            ranges(2,icol(i_l)) =
     +        max(refl(icol(i_l),nref(2)),
     +        ranges(2,icol(i_l)))
            ranges(1,icol(i_batch)) =
     +        min(refl(icol(i_batch),nref(2)),
     +        ranges(1,icol(i_batch)))
            ranges(2,icol(i_batch)) =
     +        max(refl(icol(i_batch),nref(2)),
     +        ranges(2,icol(i_batch)))
            ranges(1,icol(i_resol)) =
     +        min(refl(icol(i_resol),nref(2)),
     +        ranges(1,icol(i_resol)))
            ranges(2,icol(i_resol)) =
     +        max(refl(icol(i_resol),nref(2)),
     +        ranges(2,icol(i_resol)))
            ranges(1,icol(i_reso3)) =
     +        min(refl(icol(i_reso3),nref(2)),
     +        ranges(1,icol(i_reso3)))
            ranges(2,icol(i_reso3)) =
     +        max(refl(icol(i_reso3),nref(2)),
     +        ranges(2,icol(i_reso3)))
          end if
c
c         deal with misfits coming from XDS ( SigI < 0.0 )
c
          refl(icol(i_flag),nref(2)) = 0.0d0
          if (refl(icol(i_sigi),nref(2)).lt.0.0d0) then
            refl(icol(i_sigi),nref(2)) = -refl(icol(i_sigi),nref(2))
            refl(icol(i_flag),nref(2)) = 2.0d0
          end if
c
        end if
        goto 100
 110    rewind(lun)
      end if
c
c     do we want to skip empty runs?
c
      if (l_skip_empty_runs) then
        write(stdout,'(/)')
        do irun = 1, maxruns
          run2run(irun) = irun
          lrun2run(irun) = .false.
          nreset(irun) = 0
        end do
        numrun_off = 0

c       check for empty runs - and remove those:
        irun = 1
        do while (irun.le.numrun)

          if (run_nref(irun).eq.0) then
            do jrun = irun, (numrun - 1)
              run_num(jrun)      = run_num(jrun+1)
              run_i2n(jrun)      = run_i2n(jrun+1)
              run_batlim(1,jrun) = run_batlim(1,jrun+1)
              run_batlim(2,jrun) = run_batlim(2,jrun+1)
              run_nbatexcl(jrun) = run_nbatexcl(jrun+1)
              run_dname(jrun)    = run_dname(jrun+1)
              run_set(jrun)      = run_set(jrun+1)
              run_nref(jrun)     = run_nref(jrun+1)
            end do
            write(stdout,'(a,i4)') ' skipping empty run ',(irun
     .        +numrun_off)
            lrun2run(irun+numrun_off) = .false.
            numrun_off = numrun_off + 1
            numrun     = numrun - 1
          else
            write(stdout,'(a,i4)') ' keeping run ',(irun
     .        +numrun_off)
            run2run((irun+numrun_off)) = irun
            lrun2run(irun+numrun_off) = .true.
            irun = irun + 1
          end if
        end do

c       do we need to reset the run pointers in the reflection array?
        do iref = 1, nref(2)
          irun = int(refl(icol(i_run),iref))
          if (lrun2run(irun)) then
            refl(icol(i_run),iref) =dble(run2run(irun))
            nreset(run2run(irun)) = nreset(run2run(irun)) + 1
          end if
        end do

      end if
c
      write(stdout,'(/)')
      do jrun = 1, numrun
        write(stdout,'(a,i3,a,i8,a,i8,a,i15,a)')
     +    ' run ',                  run_i2n(jrun),
     +    ' contains batches ',run_batlim(1,jrun),
     +    ' to ',              run_batlim(2,jrun),
     +    ' and ',nreset(jrun),' reflections'
      end do
c
c     in case the limits given on a RUN card are actually outside the
c     min/max BATCH number in the reflection file:
      do irun = 1, numrun
        batlim(1) = min(batlim(1),run_batlim(1,irun))
        batlim(2) = max(batlim(2),run_batlim(2,irun))
      end do
c
c     overall/total batlim
      run_batlim(1,0) = batlim(1)
      run_batlim(2,0) = batlim(2)
c
c     calclulate good binning unless explicitly given by user (-n)
c     flag. We are using the value from the -nref flag here.
c
      if ((nshell.eq.0).and.(.not.L_BinningByEqualNumber)) then
        nshell = min(maxbin,max(minbin,nint(float(nref(2))/nrefbin)))
        if (lshell3) then
          write(stdout,'(/,a,i3,a,i6,a,/)') ' binning will be in '
     +      ,nshell,' shells of equal volume (based on ',nrefbin
     +      ,' measured reflections per shell)'
        else
          write(stdout,'(/,a,i3,a,i6,a,/)') ' binning will be in '
     +      ,nshell,' shells (based on ',nrefbin
     +      ,' measured reflections per shell)'
        end if
      else if ((nshell.eq.0).and.(L_BinningByEqualNumber)) then
        write(stdout,'(/,a,a,i6,a,/)') ' binning will be in '
     +    ,' shells of equal number of reflections (based on ',nrefbin
     +    ,' measured unique reflections per shell)'
      else if ((nshell.gt.0).and.(L_BinningByEqualNumber)) then
        write(stdout,'(/,a,i3,a,/)') ' binning will be in '
     +    ,nshell,' shells of equal number of unique reflections'
      end if
c
c     setting up resolution shells
c
      resmin = ranges(1,icol(i_resol))
      resmax = ranges(2,icol(i_resol))
      if (nshell.gt.0) resstp = (resmax - resmin) / nshell
      resmin3= ranges(1,icol(i_reso3))
      resmax3= ranges(2,icol(i_reso3))
      if (nshell.gt.0) resstp3= (resmax3 - resmin3) / nshell
c
c     same for each run
c
      if (l_ell) then
c
c       find highest resolution for ellipsoid
c
        reso = min(EllFit_res(1),min(EllFit_res(2),EllFit_res(3)))
        write(stdout,'(/,a,f9.4,a)')
     +    ' Highest resolution as defined by ellipsoid = ',reso,
     +    ' A'
c
c       take lowest resolution value of the above and the resolution
c       limit actually encountered on input
c
        if (reshigh_ell.lt.reso) then
          reshigh_ell = max(reso,reshigh_ell)
          write(stdout,'(/,2a,f9.4,a)')
     +      ' Highest resolution of available data',
     +      ' within ellipsoid = ',reshigh_ell
     +      ,' A'
        end if
      end if
      do irun = 1, numrun
c
        nshell_run(irun) = min(maxbin,max(minbin,
     +    nint(float(run_nref(irun))/nrefbin_run)))
        if (iset_nshell.eq.1) then
          nshell_run(irun) = min(nshell,nshell_run(irun))
        end if
c
        if (l_ell) then
c
c         limit to high resolution of ellipsoid definition or data
c         within ellipsoid as read:
c
          if (run_reslim(2,irun).lt.reshigh_ell) then
            run_reslim(2,irun) = max(run_reslim(2,irun),reshigh_ell)
            write(stdout,'(/,a,i4,a,f9.4,a)')
     +        ' Limiting high resolution of run ',irun,' to ',
     +        run_reslim(2,irun),' A'
          end if

        end if

        resmin_run(irun) = max(resmin,(1.0d0/run_reslim(1,irun))**2)
        resmax_run(irun) = min(resmax,(1.0d0/run_reslim(2,irun))**2)
        resstp_run(irun) = (resmax_run(irun) - resmin_run(irun))
     +                     / nshell_run(irun)
        resmin3_run(irun) = max(resmin3,(1.0d0/run_reslim(1,irun))**3)
        resmax3_run(irun) = min(resmax3,(1.0d0/run_reslim(2,irun))**3)
        resstp3_run(irun) = (resmax3_run(irun) - resmin3_run(irun))
     +                     / nshell_run(irun)

        if (L_BinningByEqualNumber_run) then
            if (iset_nshell.eq.1) then
              write(stdout,'(a,i3,a,i3,a)') ' binning for run ',
     +          run_num(irun),' will be in ',nshell_run(irun),
     +          ' shells of equal number of unique reflections'
            else
              write(stdout,'(a,i3,a,a,i6,a)') ' binning for run ',
     +          run_num(irun),' will be in',
     +          ' shells of equal number of reflections (based on ',
     +          nrefbin_run,' measured unique reflections per shell)'
            end if
        else
          if (lshell3) then
            if (iset_nshell.eq.1) then
              write(stdout,'(a,i3,a,i3,a)') ' binning for run ',
     +          run_num(irun),' will be in ',nshell_run(irun),
     +          ' shells of equal volume'
            else
              write(stdout,'(a,i3,a,i3,a,i6,a)') ' binning for run ',
     +          run_num(irun),' will be in ',nshell_run(irun),
     +          ' shells of equal volume (based on ',
     +          nrefbin_run,' measured reflections per shell)'
            end if
          else
            if (iset_nshell.eq.1) then
              write(stdout,'(a,i3,a,i3,a)') ' binning for run ',
     +          run_num(irun),' will be in ',nshell_run(irun),
     +          ' shells'
            else
              write(stdout,'(a,i3,a,i3,a,i6,a)') ' binning for run ',
     +          run_num(irun),' will be in ',nshell_run(irun),
     +          ' shells (based on ',
     +          nrefbin_run,' measured reflections per shell)'
            end if
          end if
        end if
      end do
c
      if (iverb.gt.0) then
        write (stdout,'(/,a,i8)')
     +    ' Nref (total) ................. ',nref(1)
        write (stdout,'(a,i8)')
     +    ' Nref (used) .................. ',nref(2)
        write (stdout,'(a,i3,/)')
     +    ' Ncol ......................... ',ncol
        write (stdout,'(/,a,f8.3,a,f8.3,a)')
     +    ' Resolution ................... ',sqrt(1.0d0/resmin),' - '
     +    ,sqrt(1.0d0/resmax),' A'

        if (l_ell) then
          write (stdout,'(a,f8.3,a,f8.3,a)')
     +      ' Resolution (ellipsoid) ....... ',reslow_ell,
     +      ' - '
     +      ,reshigh_ell,' A'
        end if
        write (stdout,'(a,i3,1x,a,a,a)')
     +    ' Symmetry ..................... ',symm,'"',
     +    spgrn(1:len_trim(spgrn)),'"'
        write (stdout,'(a,3f10.3,3f8.3/)')
     +    ' Cell ......................... ',cell
        do i = 1, ncoltot
          write (stdout,'(a,i3.3,a,i4,1x,a)')
     +      ' Col-',i,' ...................... ',ilabs(i),
     +      clabs(i)(1:len_trim(clabs(i)))
          if (iverb.gt.1) then
            write (stdout,'(a,i3.3,a,2g15.5)')
     +        ' Ran-',i,' ...................... ',ranges(1,i),ranges(2
     +        ,i)
          end if
        end do
        write (stdout,'(/)')
      end if
c
      xml_sg = spgrn
c
c     resort reflections
c
      if (lsort) then
        write(stdout,'(/,a,/)')
     +    ' need to re-sort reflections (H=slow, K=medium, L=fast)'
        call sort_refl(refl,ncoltot,nref(2))
      end if
c
      if (iverb.gt.2) then
        write(stdout,'(/,a,i10,a/)') ' Monitoring every ',
     +    max(1,nint(3000.0/iverb)),' reflexion'
        do iref = 1, nref(2), max(1,nint(3000.0/iverb))
          write(stdout,'(3i5,f8.3,i4,i5,2g15.5)')
     +      nint(refl(icol(i_h),iref)),
     +      nint(refl(icol(i_k),iref)),
     +      nint(refl(icol(i_l),iref)),
     +      1.d0/sqrt(refl(icol(i_resol),iref)),
     +      nint(refl(icol(i_misym),iref)),
     +      nint(refl(icol(i_batch),iref)),
     +      refl(icol(i_i),iref),
     +      refl(icol(i_sigi),iref)
        end do
      end if
c
c----+=================================================================
c
      return
      end
