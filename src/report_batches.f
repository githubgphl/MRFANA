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

c     report batch information for run(s) similar to SCALA /AIMLESS

      subroutine report_batches(runs,nrun,batch_info,batch_info_lim)
c
      implicit none
c
      include 'params.inc'
      include 'runs.inc'
      include 'units.inc'
      include 'info.inc'
      include 'batch.inc'
      include 'xml.inc'
c
      integer nrun
      integer runs(nrun), batch_info_lim(2)
      double precision batch_info(n_batch_info
     +  ,batch_info_lim(1):batch_info_lim(2))
c
      integer ibat, ncell, i, j, iex, irun, jrun, istat
      logical ldo
      double precision  cellm(6), cellsd(6), wavelim(2), osclim(2
     +  ,maxruns)
      double precision cell1(6), mean, var
c
      double precision, allocatable :: batch_cell(:,:)
c
      integer nbat, jumpax
      character axes(4)*4
      data axes/'none','a*','b*','c*'/
c
c----+=================================================================
c
      do i = 1, 6
        cell1(i) = 0.0d0
        cellm(i) = 0.0d0
        cellsd(i) = 0.0d0
      end do
      ncell = 0
      wavelim(1) = 999999.9d0
      wavelim(2) = 0.0d0
c
      nbat = 0
      do irun = 1, nrun
        jrun = runs(irun)
        nbat = nbat + ( run_batlim(2,jrun) - run_batlim(1,jrun) + 1 )
        osclim(1,irun) = 999999.9d0
        osclim(2,irun) = 0.0d0
      end do
      allocate(batch_cell(1:6,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit
     .  ('unable to allocate batch_cell array in s/r report_batches')
c
c     loop over all runs
c
      do irun = 1, nrun
        jrun = runs(irun)
c
c     loop over all batches in that run
c
        do ibat = run_batlim(1,jrun), run_batlim(2,jrun)

c     check if this is an excluded batch:
          ldo = .true.
          do iex = 1, run_nbatexcl(jrun)
            if (run_batexcl(iex,jrun).eq.ibat) then
              ldo = .false.
              exit
            end if
          end do

c     not an excluded batch:
          if (ldo) then

c     we have some cell values for it (if there were intermittent
c     batches missing, they won't have a cell value)
            if (batch_info(i_batch_info_cell,ibat).gt.0.) then

c     get cell parameters
              ncell = ncell + 1
              do i = 1, 6
                batch_cell(i,ncell) =
     +            batch_info((i_batch_info_cell+i-1),ibat)
                cell1(i) = cell1(i) +
     +            batch_info((i_batch_info_cell+i-1),ibat)
              end do
c     get wavelength
              if (wavelim(1).gt.wavelim(2)) then
                wavelim(1) = batch_info(i_batch_info_wave,ibat)
                wavelim(2) = batch_info(i_batch_info_wave,ibat)
              end if
              wavelim(1) =
     +          min(wavelim(1),batch_info(i_batch_info_wave,ibat))
              wavelim(2) =
     +          max(wavelim(1),batch_info(i_batch_info_wave,ibat))
c     get osc range
              if (osclim(1,irun).gt.osclim(2,irun)) then
                osclim(1,irun) =
     +            min(batch_info( i_batch_info_phistart   ,ibat),
     +                batch_info((i_batch_info_phistart+1),ibat))
                osclim(2,irun) =
     +            max(batch_info( i_batch_info_phistart   ,ibat),
     +                batch_info((i_batch_info_phistart+1),ibat))
              end if
              osclim(1,irun) = min(osclim(1,irun),
     +          min(batch_info( i_batch_info_phistart   ,ibat),
     +              batch_info((i_batch_info_phistart+1),ibat)))
              osclim(2,irun) = max(osclim(2,irun),
     +          max(batch_info( i_batch_info_phistart   ,ibat),
     +              batch_info((i_batch_info_phistart+1),ibat)))
            end if
          end if
        end do
c
c     report cell dimensions (and sd)
c
      end do
      if (ncell.gt.0) then
        do i = 1, 6
c     get mean
          mean = cell1(i) / ncell
          cell1(i) = 0.0d0
c     sum (X-<X>)**2
          do j = 1, ncell
            cell1(i) = cell1(i) + ( (batch_cell(i,j)-mean) **2 )
          end do
c     variance
          var = cell1(i) / (ncell - 1)
c     store mean and sd
          cellm(i)    = mean
          cellsd(i)   = sqrt(var)
          xml_cell(i) = cellm(i)
        end do
        write(stdout,'(3x,a,6f8.2)') 'Cell:',cellm
        ldo=.false.
        do i = 1, 6
          if (cellsd(i).ne.0.d0) then
            ldo = .true.
            exit
          end if
        end do
        if (ldo) then
          write(6,'(3x,a,6F8.4)') '     ',cellsd
        end if
      end if
      if (wavelim(1).lt.wavelim(2)) then
        write(stdout,'(3x,2(a,F9.5),a)')
     +    'Wavelength:',wavelim(1),' - ',wavelim(2),' A'
        xml_wave = (wavelim(1)+wavelim(2))/2
      else if (wavelim(1).lt.999999.) then
        write(stdout,'(3x,a,F9.5,a)')
     +    'Wavelength:',wavelim(1),' A'
        xml_wave = wavelim(1)
      end if
      do irun = 1, nrun
        jrun = runs(irun)
        write(stdout,'(3x,a,i4,a,i10,a,i10)')
     +    'Run number: ',max(1,run_num(jrun)),' consists of batches ',
     +    run_batlim(1,jrun),' to ',run_batlim(2,jrun)
        if (jrun.gt.0) then
          write(stdout,'(6x,a,2f8.2)')
     +      'Resolution range for run: ',
     +      1./sqrt(resmin_run(jrun)),
     +      1./sqrt(resmax_run(jrun))
        end if
        if (osclim(1,irun).ne.osclim(2,irun)) then
          write(stdout,'(6x,2(a,f19.2))')
     +      'Phi range: ',osclim(1,irun),' to ',osclim(2,irun)
        end if

        jumpax=nint(batch_info(i_batch_info_jumpax,run_batlim(1,jrun)))
        if (jumpax.gt.0) then
          write(stdout,'(6x,3a,f6.1,a)')
     +      'Closest reciprocal axis to spindle: ',axes(jumpax+1)
     .      ,' (angle ',batch_info(i_batch_info_angle,run_batlim(1
     .      ,jrun)),' degrees)'
        else
          write(stdout,'(6x,2a)')
     +      'Closest reciprocal axis to spindle: ',axes(jumpax+1)
        end if
      end do
      write(stdout,'(/)')
c
      deallocate(batch_cell,stat=istat)
      if (istat.ne.0) call error_exit
     .  ('unable to deallocate batch_cell array in s/r report_batches')
c
c----+=================================================================
c
      return
      end
