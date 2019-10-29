c     **********************************************************************
c
c      Copyright (c) 2007, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein, Claus Flensburg, Ian J. Tickle and
c               Gerard Bricogne
c
c     **********************************************************************

      program mrfana

c     analyses reflection files that have unmerged (multiple
c     records per unique HKL) data

      implicit none
c
      include 'args.inc'
      include 'params.inc'
      include 'runs.inc'
      include 'units.inc'
      include 'info.inc'
      include 'batch.inc'
      include 'xml.inc'
c
      character filtyp*13
      character*30 clabs(maxcol)
      character*1024 runs_in_dname
      integer i0_in_dname,i1_in_dname,i2_in_dname
      integer nref(2), ncol, ncoltot, ilabs(maxcol),
     +  i
      double precision ranges(2, maxcol)
      integer i_dname, n_dname, ntotal_d
c
      double precision, allocatable :: refl(:,:)
      double precision, allocatable :: batch_info(:,:)
      integer istat
c
      integer irun(maxruns), jrun
c
c----+=================================================================
c

c     initialisation routine
      call init(nref,ncol,ranges,clabs,ilabs)

c     read command-line arguments (and cards on stdin)
      call get_args()

c     open reflection data file
      write(stdout,*) ' Calling open_file ... '
      call open_file(filtyp)

c     get information about data
      write(stdout,*) ' Calling get_info ... '
      call get_info(filtyp,nref,ncol,ncoltot,ranges,clabs,ilabs)

c     allocate reflection data array
      allocate(refl(1:ncoltot,1:nref(1)),stat=istat)
      if (istat.ne.0) call error_exit('unable to allocate refl array')

c     read reflection data
      write(stdout,*) ' Calling read_refl ... '
      call read_refl(filtyp,refl,ncol,ncoltot,ranges,nref,clabs,ilabs)

c     since we want to be able to handle multi-sweep/multi-crystal data,
c     we need to handle the concept of batches (aka images)
      allocate(batch_info(1:n_batch_info,batlim(1):batlim(2)),stat
     +  =istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate batch_info array')

c 
      write(stdout,*) ' Calling read_batch ... '
      call read_batch(filtyp,batch_info,batlim)
c
c     now we have all reflections read in - let's do something with them
      write(stdout,*) ' Calling stats ... '
c
c     ------------------------------------------------------------------

      call write_info()
c
c     ------------------------------------------------------------------
c
c     loop over all runs
c
      do jrun = 1, numrun
        if (numrun.gt.1) then
          write(stdout,'(/,2a)')  ' ##############################',
     +                             '##############################'
          write(stdout,'(a,i4)')  ' #################### RUN ',
     +      run_num(jrun)

          write(stdout,'(2a,/)')  ' ##############################',
     +                             '##############################'
        end if
        irun(1) = jrun
c
        if (xml_base.ne.' ') then
          if (run_num(jrun).le.9) then
            write(xml_file,'(2a,i1,a)')
     +        trim(xml_base),'_run',run_num(jrun),'.xml'
          else if (run_num(jrun).le.99) then
            write(xml_file,'(2a,i2,a)')
     +        trim(xml_base),'_run',run_num(jrun),'.xml'
          else if (run_num(jrun).le.999) then
            write(xml_file,'(2a,i3,a)')
     +        trim(xml_base),'_run',run_num(jrun),'.xml'
          end if
        end if
c
        call report_batches(irun,1,batch_info,batlim)
c
c       merging stats (as function of resolution):
        call stats(refl,ncoltot,nref,1,irun)
c       various stats (as function of image/batch number):
        call cumcmpl(refl,ncoltot,nref,ranges,1,irun)
c
      end do
c
c     ------------------------------------------------------------------
c
c     several runs: get stats for total(s) depending on DNAME
c
      ntotal_d = 0
c
      if (numrun.ne.1) then
        i_dname = 1
        do while (i_dname.gt.0)
          i_dname = 0
          n_dname = 0
          runs_in_dname = ""
          i0_in_dname = 0
          do jrun = 1, numrun
            if (run_dname(jrun).ne.' ') then
              if (i_dname.eq.0) then
                i_dname = jrun
                n_dname = n_dname + 1
                irun(n_dname) = jrun
c               store run numbers:
                i0_in_dname = i0_in_dname + 1
                i1_in_dname = ( ( i0_in_dname - 1 ) * 4 ) + 1
                i2_in_dname =     i0_in_dname       * 4 
                write(runs_in_dname(i1_in_dname:i2_in_dname),'(1x,i3)')
     +            run_num(irun(n_dname))
              else
                if (run_dname(jrun).eq.run_dname(i_dname)) then
                  n_dname = n_dname + 1
                  irun(n_dname) = jrun
c                 store run numbers:
                  i0_in_dname = i0_in_dname + 1
                  i1_in_dname = ( ( i0_in_dname - 1 ) * 4 ) + 1
                  i2_in_dname =     i0_in_dname       * 4 
                write(runs_in_dname(i1_in_dname:i2_in_dname),'(1x,i3)')
     +              run_num(irun(n_dname))
                end if
              end if
            end if
          end do
c

c
          if (n_dname.gt.0) then
            ntotal_d = ntotal_d + 1
            write(stdout,'(/,2a)')  ' ##############################',
     +                              '##############################'
            write(stdout,'(a,a)')  ' #################### TOTAL: '
     +        ,run_dname(i_dname)(1:len_trim(run_dname(i_dname)))
            write(stdout,'(2a,/)')  ' ##############################',
     +                              '##############################'
c
            write(stdout,'(4a,/)') ' Dataset '
     +        ,run_dname(i_dname)(1:len_trim(run_dname(i_dname)))
     +        ,' consists of runs:'
     +        ,runs_in_dname(1:len_trim(runs_in_dname))
c
            call report_batches(irun,n_dname,batch_info,batlim)
c
            if (xml_base.ne.' ') then
              write(xml_file,'(4a)')
     +          trim(xml_base),'_',
     +          run_dname(i_dname)(1:len_trim(run_dname(i_dname))),
     +          '.xml'
            end if
c
            call stats(refl,ncoltot,nref,n_dname,irun)
            call cumcmpl(refl,ncoltot,nref,ranges,n_dname,irun)
c
            do jrun = 1, n_dname
              run_dname(irun(jrun)) = ' '
            end do
          end if
        end do

        if (ntotal_d.ne.1) then
c
c         only need to run this one if we have multiple sweeps with
c         different DNAME
          write(stdout,'(/,2a)')  ' ##############################',
     +                            '##############################'
          write(stdout,'(a,i4)')  ' #################### TOTAL '
          write(stdout,'(2a,/)')  ' ##############################',
     +                             '##############################'
c
          if (xml_base.ne.' ') then
            write(xml_file,'(2a)')
     +        trim(xml_base),'_total.xml'
          end if
c
          do i = 1, numrun
            irun(i) = i
          end do

          call report_batches(irun,max(numrun,1),batch_info,batlim)
c
          irun(1) = 0
          call stats(refl,ncoltot,nref,1,irun)
          call cumcmpl(refl,ncoltot,nref,ranges,1,irun)
c
        end if
c
      end if
c
c     ------------------------------------------------------------------
c
      write(stdout,*) ' Calling close_file ... '
      call close_file(filtyp)
      deallocate(refl,stat=istat)
      if (istat.ne.0) call error_exit('unable to deallocate refl array')
      deallocate(batch_info,stat=istat)
      if (istat.ne.0) call error_exit
     +  ('unable to deallocate batch_info array')
c
      write(stdout,'(/,a,/)') ' MRFANA: normal termination'
c
c----+=================================================================
c
      end
c
c----+=================================================================
c
