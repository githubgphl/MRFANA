c     **********************************************************************
c
c      Copyright (c) 2015, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein, Claus Flensburg and Gerard Bricogne
c
c     **********************************************************************

c     routine to compute various metrics as a function of batch (i.e
c     image) number

      subroutine cumcmpl(refl,ncoltot,nref,ranges,nrun,irun)
c
      implicit none
c
      include 'params.inc'
      include 'cmpl.inc'
      include 'lookup.inc'
      include 'reso.inc'
      include 'units.inc'
      include 'aniso.inc'
      include 'args.inc'
c
      integer ncoltot,nref(2),nrun,irun(nrun)
      double precision refl(ncoltot,nref(1))
      double precision ranges(2,maxcol)
c
      integer i, iref, jref, this_bat, this_set, j, istat,
     .  isym,isym_sum, isym_prod, ipm, ipm1, ipm2, ipm3, jpm, k, iwgt
     .  ,jnpos(12), nrefpos(6), nreftot(2)
c
      integer ibat_min, ibat_max, nbat, hkl(3), ic, hklold(3), ic_old,
     +  nhkl(3), num_intensity_mrg(3),
     +  iset_min, iset_max, isysabs, isysabs_old, ellval
      real epsi
      double precision w, AbsDiff_intensity(2), wtmp
      double precision intensity_mrg(3), sig_intensity_mrg(3),
     .  intensity(3,maxmul), sig_intensity(3,maxmul), sum_intensity(3),
     .  wgt(2), var_intensity_mrg(3), var_intensity(3,maxmul),
     .  i_mrg, s_mrg
      integer jbatch(3,maxmul), jset(3,maxmul)
      double precision compl_ano_p, compl_ano_m,
     +                 compl_ano_ell_p, compl_ano_ell_m
c      logical l_flag
c
      integer, allocatable :: bat2i(:,:)
      integer, allocatable :: i2bat(:)
      integer, allocatable :: i2set(:)
      integer, allocatable :: inpos(:,:)
      double precision, allocatable :: store(:,:)
c     Rmerge, Rmeas, Rpim:
      double precision, allocatable :: rval(:,:,:,:)
      double precision rvaltmp(3,2)
      double precision ivaltmp(2)
c     <I>, <I/sigI>
      double precision, allocatable :: ival(:,:,:)
      double precision, allocatable ::     intensity_mrg_batch(:,:)
      double precision, allocatable :: sig_intensity_mrg_batch(:,:)
      double precision, allocatable :: var_intensity_mrg_batch(:,:)
      integer         , allocatable :: num_intensity_mrg_batch(:,:)

c     pointers - to make it better readable
      integer iAll, iPlus, iMinus
      parameter (iAll=1, iPlus=2, iMinus=3)
c
      integer njref
      integer, allocatable :: jref2iref(:)
      integer :: jref2ibin(1)
      integer n_jref2ibin
      double precision :: jref2reso(1:2,1:(maxshell+1)), dp_tmp(1)
c
      character*13 c_cmpl(2)
c
c----+=================================================================
c
      if (l_ell) then
        c_cmpl(1) = 'Compl. Spher.'
      else
        c_cmpl(1) = 'Completeness '
      end if
      c_cmpl(2) = 'Compl. Ellip.'
c
      allocate(jref2iref(1:(nref(2)+1)),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate jref2iref array')

      if ((nrun.eq.1).and.(irun(1).gt.0)) then
        call refl_test_active(ncoltot,nref,refl,nrun,irun,
     +    jref2iref,njref,
     +    .false.,0,dp_tmp,nrefbin_run, jref2ibin, n_jref2ibin,
     +    jref2reso)
      else
        call refl_test_active(ncoltot,nref,refl,nrun,irun,
     +    jref2iref,njref,
     +    .false.,0,dp_tmp,nrefbin    , jref2ibin, n_jref2ibin,
     +    jref2reso)
      end if
c
c     these are for testing if we have both I+ and I- measured:
c     ISYM values are odd/even for I+/I- and we should have then
c     ( http://www.ccp4.ac.uk/html/mtzformat.html#itemnames )
c
c       isym_sum   > 0
c       isym_prod == 0
c
c     when using mod(isym,2) for those calculations. If we have at least
c     one I+ measurement, isym_sum will be larger than zero. And if we
c     have at least one I- measurement isym_prod will be zero.

      isym_sum  = 0
      isym_prod = 1
      nhkl = 0
      num_intensity_mrg = 0
          intensity_mrg = 0.0d0
      var_intensity_mrg = 0.0d0
      sig_intensity_mrg = 0.0d0

c
      if (icol(i_iset).gt.0) then
        allocate(bat2i(
     +    nint(ranges(1,icol(i_batch))):
     +    nint(ranges(2,icol(i_batch))) ,
     +    nint(ranges(1,icol(i_iset))):
     +    nint(ranges(2,icol(i_iset))) ),stat=istat)
      else
        allocate(bat2i(
     +    nint(ranges(1,icol(i_batch))):
     +    nint(ranges(2,icol(i_batch))) ,
     +    0:0),stat=istat)
      endif
      if (istat.ne.0)
     +  call error_exit('unable to allocate bat2i array')
c
      bat2i = 0
c      do i = nint(ranges(1,icol(i_batch))), nint(ranges(2
c     .  ,icol(i_batch)))
c        bat2i(i) = 0
c      end do

      this_set = 0
c     loop over all reflections
      j = 0
      do jref = 1, njref
        iref = jref2iref(jref)

        this_bat = nint(refl(icol(i_batch),iref))
        if (icol(i_iset).gt.0) then
          this_set = nint(refl(icol(i_iset),iref))
        endif

        j = j + 1
        if (j.eq.1) then
          ibat_min = this_bat
          ibat_max = this_bat
          iset_min = this_set
          iset_max = this_set
        end if
        ibat_min = min(this_bat,ibat_min)
        ibat_max = max(this_bat,ibat_max)
        bat2i(this_bat,this_set) = bat2i(this_bat,this_set) + 1
        iset_min = min(this_set,iset_min)
        iset_max = max(this_set,iset_max)
      end do
c
c
      nbat = 0
      do j = iset_min,iset_max
        do i = ibat_min,ibat_max
          if (bat2i(i,j).gt.0) nbat=nbat + 1
        end do
      end do

      allocate(i2bat(1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate i2bat array')
      i2bat = 0
      allocate(i2set(1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate i2set array')
      i2set = 0
      allocate(inpos(1:8,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate inpos array')
      inpos = 0

      nbat = 0
      do j = iset_min,iset_max
        do i = ibat_min,ibat_max
          if (bat2i(i,j).gt.0) then
            nbat=nbat + 1
            i2bat(nbat) = i
            i2set(nbat) = j
            bat2i(i,j) = nbat
          end if
        end do
      end do
c
c     now we have arrays to go from index to batch number and back
      allocate(store(1:2,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate store array')
      store = 0.0d0

      allocate(rval(1:5,1:2,1:2,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate rval array')
      rval = 0.0d0

      allocate(ival(1:3,1:2,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate ival array')
      ival = 0.0d0

      allocate(intensity_mrg_batch(1:3,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate intensity_mrg_batch array')
      intensity_mrg_batch = 0.0d0

      allocate(var_intensity_mrg_batch(1:3,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit(
     +  'unable to allocate var_intensity_mrg_batch array')
      var_intensity_mrg_batch = 0.0d0

      allocate(sig_intensity_mrg_batch(1:3,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit(
     +  'unable to allocate sig_intensity_mrg_batch array')
      sig_intensity_mrg_batch = 0.0d0

      allocate(num_intensity_mrg_batch(1:3,1:nbat),stat=istat)
      if (istat.ne.0)
     +  call error_exit(
     +  'unable to allocate num_intensity_mrg_batch array')
      num_intensity_mrg_batch = 0

c     =================================================================
c     loop over all reflections
      njref = njref + 1
      jref2iref(njref) = nref(2) + 1
      do jref = 1, njref
        iref = jref2iref(jref)

        if (iref.le.nref(2)) then
          this_bat = nint(refl(icol(i_batch),iref))
          if (icol(i_iset).gt.0) then
            this_set = nint(refl(icol(i_iset),iref))
          end if

          j = bat2i(this_bat,this_set)

c         things we can store on a per-measurement basis
c         <I/sigI>
          if (refl(icol(i_sigi),iref).gt.0.d0) then
            store(1,j) = store(1,j) + 1
            store(2,j) = store(2,j) +
     +        sngl(refl(icol(i_i),iref)/refl(icol(i_sigi),iref))
          end if
c
        end if
c
c       use this loop also for the final reflection ... ugly, but easier
c
        if (iref.gt.nref(2)) then
c         force popping of stash by setting 'new' reflection indices
          hkl(1) = 0
          hkl(2) = 0
          hkl(3) = 0
        else
c
c         additional information about this reflection:
c
          hkl(1) = nint(refl(icol(i_h),iref))
          hkl(2) = nint(refl(icol(i_k),iref))
          hkl(3) = nint(refl(icol(i_l),iref))
c
          call centr(hkl,ic)
c
c         check for systematic absences:
          isysabs = 0
          call epslon(hkl,epsi,isysabs)

        end if
c
c       --------------------------------------------------------------
c       is this a new unique reflection?
c       --------------------------------------------------------------
c
        if ( ( (hkl(1).ne.hklold(1)) .or.
     +         (hkl(2).ne.hklold(2)) .or.
     +         (hkl(3).ne.hklold(3))     ).and.
     +       (iref.gt.1                  )     )then
c
c         ============================================================
c         first process the previous one
c         ============================================================
c
          if ((nhkl(iAll).gt.0).and.(isysabs_old.eq.0)) then
c
c           is this a centric one or not?
            call centr(hklold,ic_old)
c           1 = all observations = iAll
c           2 = I+ observations  = iPlus
c           3 = I- observations  = iMinus
c
            if (l_ell) call ellchk(hklold,ellval)

c           get merged intensities (and sigmas) for all
            do ipm = iAll, iMinus
              if (num_intensity_mrg(ipm).gt.0) then
                intensity_mrg(ipm) = 
     +          intensity_mrg(ipm) / var_intensity_mrg(ipm)
                sig_intensity_mrg(ipm) =
     +          sqrt(1.d0/var_intensity_mrg(ipm))
              else
                intensity_mrg(ipm)     = 0.0d0
                sig_intensity_mrg(ipm) = 0.0d0
                var_intensity_mrg(ipm) = 0.0d0
              end if
            end do

c           for cumulative completeness:
            jnpos(1) = 0
            jnpos(2) = 0
            do j = 1, nhkl(iAll)
              k = bat2i(jbatch(iAll,j),jset(iAll,j))
              if (j.eq.1) then
                jnpos(1) = k
                jnpos(5) = k
              endif
              jnpos(1) = min(jnpos(1),k)
              jnpos(5) = max(jnpos(5),k)
            end do
            if (jnpos(1).gt.0) then
              inpos(1,jnpos(1)) = inpos(1,jnpos(1)) + 1
              inpos(3,jnpos(5)) = inpos(3,jnpos(5)) + 1
              jnpos(3) = jnpos(3) + 1
              if (l_ell.and.(ellval.eq.1)) then
                inpos(5,jnpos(1)) = inpos(5,jnpos(1)) + 1
                inpos(7,jnpos(5)) = inpos(7,jnpos(5)) + 1
                jnpos(9) = jnpos(9) + 1
              end if
            end if

c           anomalous
            if ( (ic_old.eq.0)      .and.
     +           (nhkl(iPlus).gt.0) .and.
     +           (nhkl(iMinus).gt.0)     ) then
              jnpos(1) = 0
              jnpos(2) = 0
              jnpos(5) = 0
              jnpos(6) = 0
              do j = 1, nhkl(iPlus)
                k = bat2i(jbatch(iPlus,j),jset(iPlus,j))
                if (j.eq.1) then
                  jnpos(1) = k
                  jnpos(5) = k
                end if
                jnpos(1) = min(jnpos(1),k)
                jnpos(5) = max(jnpos(5),k)
              end do
              do j = 1, nhkl(iMinus)
                k = bat2i(jbatch(iMinus,j),jset(iMinus,j))
                if (j.eq.1) then
                  jnpos(2) = k
                  jnpos(6) = k
                end if
                jnpos(2) = min(jnpos(2),k)
                jnpos(6) = max(jnpos(6),k)
              end do
              if ((jnpos(1).gt.0).and.(jnpos(2).gt.0)) then
                inpos(2,max(jnpos(1),jnpos(2))) =
     +          inpos(2,max(jnpos(1),jnpos(2))) + 1
                inpos(4,min(jnpos(5),jnpos(6))) =
     +          inpos(4,min(jnpos(5),jnpos(6))) + 1
                jnpos(4) = jnpos(4) + 1
                if (l_ell.and.(ellval.eq.1)) then
                  inpos(6,max(jnpos(1),jnpos(2))) =
     +            inpos(6,max(jnpos(1),jnpos(2))) + 1
                  inpos(8,min(jnpos(5),jnpos(6))) =
     +            inpos(8,min(jnpos(5),jnpos(6))) + 1
                  jnpos(10) = jnpos(10) + 1
                end if
              end if
            end if

            if (nhkl(iAll).gt.0) then
c
c             1 = normal statistics, treating I+ and I- identical
c             2 = treating I+ and I- separate
c
              do ipm = iAll, iMinus
c
c               <I> and <I/sigI> stats
c
                if (nhkl(ipm).gt.0) then

c                 collaps both of I+ and I- into one set of stats - with
c                 the only difference that we treat them separately:
                  jpm = min(ipm,2)
c
c                 loop over all measurements for that class:
                  do j = 1, nhkl(ipm)
c
c                   pointer to batch index
c
                    k = bat2i(jbatch(ipm,j),jset(ipm,j))

c                   get merged stats for all measurements (of this type)
c                   that fall into that batch:
                    if (num_intensity_mrg_batch(ipm,k).gt.0) then
                      i_mrg = intensity_mrg_batch(ipm,k) /
     +                    var_intensity_mrg_batch(ipm,k)
                      s_mrg = sqrt(1.0d0/var_intensity_mrg_batch(ipm,k))

                      ival(1,jpm,k) = ival(1,jpm,k) + i_mrg
                      ival(2,jpm,k) = ival(2,jpm,k) + i_mrg/s_mrg
                      ival(3,jpm,k) = ival(3,jpm,k) + 1

c                     mark this as being done - so we don't do it twice!
                      num_intensity_mrg_batch(ipm,k) = 0
                          intensity_mrg_batch(ipm,k) = 0.0d0
                      var_intensity_mrg_batch(ipm,k) = 0.0d0
                      sig_intensity_mrg_batch(ipm,k) = 0.0d0
                    else
c                     just in case:
                      num_intensity_mrg_batch(ipm,k) = 0
                          intensity_mrg_batch(ipm,k) = 0.0d0
                      var_intensity_mrg_batch(ipm,k) = 0.0d0
                      sig_intensity_mrg_batch(ipm,k) = 0.0d0
                    end if
c
                  end do
c
                end if

c               we can only compute R-values if we have more than one
c               measurement:
                if (nhkl(ipm).gt.1) then

c                 collaps both of I+ and I- into one set of stats - with
c                 the only difference that we treat them separately:
                  jpm = min(ipm,2)
c
c                 loop over all measurements for that class:
                  do j = 1, nhkl(ipm)
c
c                   pointer to batch index
c
                    k = bat2i(jbatch(ipm,j),jset(ipm,j))
c
c                   sums for Rmerge, Rmeas and Rpim
c
c                   unweighted,weighted
                    wgt(1) = 1.0d0
                    wgt(2) = 1.0d0/var_intensity(ipm,j)

c                   loop over unweighted, weighted
                    do iwgt = 1, 2
                      AbsDiff_intensity(iwgt) = wgt(iwgt) *
     +                  abs(intensity(ipm,j) - intensity_mrg(ipm))

c                     Rmerge
                      rval(1,iwgt,jpm,k) = 
     +                rval(1,iwgt,jpm,k) + AbsDiff_intensity(iwgt)

c                     Rmeas
                      wtmp = sqrt(float(nhkl(ipm))/(nhkl(ipm)-1))
                      rval(2,iwgt,jpm,k) = 
     +                rval(2,iwgt,jpm,k) + AbsDiff_intensity(iwgt)*wtmp

c                     Rpim
                      wtmp = sqrt(              1./(nhkl(ipm)-1))
                      rval(3,iwgt,jpm,k) = 
     +                rval(3,iwgt,jpm,k) + AbsDiff_intensity(iwgt)*wtmp

c                     denominator
                      rval(4,iwgt,jpm,k) = 
     +                rval(4,iwgt,jpm,k) + wgt(iwgt)*intensity_mrg(ipm)

c                     number
                      rval(5,iwgt,jpm,k) = 
     +                rval(5,iwgt,jpm,k) + 1

                    end do ! end loop over unweighted,weighted

                  end do ! end loop over all measurements for that class
c
                end if

              end do ! end loop over iGeneral, iAnomalous
c
            end if
c
          end if
c
c         ============================================================
c         finished processing the previous one
c         ============================================================
c
c         we are after the last reflection - exit loop
          if (i.gt.nref(2)) cycle
c
          isym_sum  = 0
          isym_prod = 1
          nhkl = 0
          num_intensity_mrg = 0
              intensity_mrg = 0.0d0
          sig_intensity_mrg = 0.0d0
          var_intensity_mrg = 0.0d0
          num_intensity_mrg_batch = 0
              intensity_mrg_batch = 0.0d0
          sig_intensity_mrg_batch = 0.0d0
          var_intensity_mrg_batch = 0.0d0
c
        end if
c
c       ============================================================
c       start processing the new one
c       ============================================================
c
c       we're already past the last reflection:
        if ( (hkl(1).eq.0).and.
     +       (hkl(2).eq.0).and.
     +       (hkl(3).eq.0)     ) cycle

c
c       get ISYM:
c
c       odd  values are I(+)
c       even values are I(-)
c
c       each friedel pair is consecutive
c       i.e. isym=(h,k,l) and isym+1=(-h,-k,-l)
c
c       m    = refl(icol(i_misym),iref)/256
c       isym = refl(icol(i_misym),iref) - m * 256
c
        isym = nint(refl(icol(i_misym),iref))
        isym_sum  = isym_sum  + mod(isym,2)
        isym_prod = isym_prod * mod(isym,2)

c
c         process this measurement (for a given unique reflection)
c
c         1 = all observations = iAll
c         2 = I+ observations  = iPlus
c         3 = I- observations  = iMinus

c         define loop (to store in iAll and - if applicable - also iPlus
c         and/or iMinus)

          ipm1 = iAll
          ipm2 = iAll
          if (ic.eq.0) then
c           pointer to
c           2 = I(+)    or
c           3 = I(-)
            ipm2 = mod(isym+1,2)+2
          end if
          ipm3 = ipm2 - ipm1
          if (ipm3.lt.1) ipm3 = 1
c
          do ipm = ipm1, ipm2, ipm3

c           increment number of measurements:
            nhkl(ipm) = nhkl(ipm) + 1

            jbatch(ipm,nhkl(ipm)) = nint(refl(icol(i_batch),iref))
            if (icol(i_iset).gt.0) then
              jset(ipm,nhkl(ipm)) = nint(refl(icol(i_iset ),iref))
            else
              jset(ipm,nhkl(ipm)) = 0
            end if

c           fetch intensity and sigma
            intensity    (ipm,nhkl(ipm)) = refl(icol(i_i   ),iref)
            sig_intensity(ipm,nhkl(ipm)) = refl(icol(i_sigi),iref)
            var_intensity(ipm,nhkl(ipm)) =
     +      sig_intensity(ipm,nhkl(ipm))**2

c           sum intensities
            sum_intensity(ipm) =
     +      sum_intensity(ipm) + intensity(ipm,nhkl(ipm))

            w = 1.0d0/var_intensity(ipm,nhkl(ipm))

                intensity_mrg(ipm) = 
     +          intensity_mrg(ipm) + w*intensity(ipm,nhkl(ipm))
            var_intensity_mrg(ipm) = 
     +      var_intensity_mrg(ipm) + w
            num_intensity_mrg(ipm) = 
     +      num_intensity_mrg(ipm) + 1

            k = bat2i(jbatch(ipm,nhkl(ipm)),jset(ipm,nhkl(ipm)))
                intensity_mrg_batch(ipm,k) = 
     +          intensity_mrg_batch(ipm,k) + w*intensity(ipm,nhkl(ipm))
            var_intensity_mrg_batch(ipm,k) = 
     +      var_intensity_mrg_batch(ipm,k) + w
            num_intensity_mrg_batch(ipm,k) = 
     +      num_intensity_mrg_batch(ipm,k) + 1

          end do

c
        do j = 1, 3
          hklold(j) = hkl(j)
        end do
        isysabs_old   = isysabs
c
      end do

      cmpl_nshell = 0
      call cmpl(current_reslow,current_reshigh,0,nrefpos)

      write(stdout,'(/,a,/,a,/)')
     +  ' W0 = unweighted',
     +  ' W1 = 1/(sigI**2) weighted'

      if ((iset_max-iset_min).gt.0) then
        write(stdout,'(3a)',advance='no')
     +    ' -----------------------------------------------------',
     +    '-----------------------------------------------------',
     +    '----------------------------------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '------------------------------------'
        endif
        write(stdout,'(/)',advance='no')

        write(stdout,'(7a)',advance='no')
     +    '                                ',
     +    '       merged        ',
     +    '        ',
     +    '         Rmerge ',
     +    '           Rmeas  ',
     +    '            Rpim',
     +    '                  '//c_cmpl(1)
        if (l_ell) then
          write(stdout,'(19x,a)',advance='no') c_cmpl(2)//'  '
        endif
        write(stdout,'(/)',advance='no')

        write(stdout,'(8a)',advance='no')
     +    '                           ',
     +    ' --------------------',
     +    '        ',
     +    '    ----------------',
     +    '  ----------------',
     +    '  ----------------',
     +    '  ----------------',
     +    '------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '  ----------------------------------'
        endif
        write(stdout,'(/)',advance='no')
        write(stdout,'(7a)',advance='no')
     +    '    BATCH ISET     Nref <I/sigI>',
     +    '      <I>    <I/sigI>     #ref',
     +    '      W0      W1',
     +    '        W0      W1',
     +    '        W0      W1  ',
     +    '       All     Ano',
     +    '      -All    -Ano  '
        if (l_ell) then
          write(stdout,'(2a)',advance='no')
     +    '     All     Ano',
     +    '      -All    -Ano  '
        endif
        write(stdout,'(/)',advance='no')

        write(stdout,'(3a)',advance='no')
     +    ' -----------------------------------------------------',
     +    '-----------------------------------------------------',
     +    '----------------------------------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '------------------------------------'
        endif
        write(stdout,'(/)',advance='no')

      else

        write(stdout,'(3a)',advance='no')
     +    ' -----------------------------------------------------',
     +    '-----------------------------------------------------',
     +    '-----------------------------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '------------------------------------'
        endif
        write(stdout,'(/)',advance='no')

        write(stdout,'(7a)',advance='no')
     +    '                           ',
     +    '       merged        ',
     +    '        ',
     +    '         Rmerge ',
     +    '           Rmeas  ',
     +    '            Rpim',
     +    '                  '//c_cmpl(1)
        if (l_ell) then
          write(stdout,'(18x,a)',advance='no') c_cmpl(2)//'  '
        endif
        write(stdout,'(/)',advance='no')

        write(stdout,'(8a)',advance='no')
     +    '                           ',
     +    ' --------------------',
     +    '        ',
     +    '    ----------------',
     +    '  ----------------',
     +    '  ----------------',
     +    '  ----------------',
     +    '------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '  ----------------------------------'
        endif
        write(stdout,'(/)',advance='no')

        write(stdout,'(7a)',advance='no')
     +    '    BATCH     Nref <I/sigI>',
     +    '      <I>    <I/sigI>     #ref',
     +    '      W0      W1',
     +    '        W0      W1',
     +    '        W0      W1  ',
     +    '       All     Ano',
     +    '      -All    -Ano  '
        if (l_ell) then
          write(stdout,'(2a)',advance='no')
     +    '     All     Ano',
     +    '      -All    -Ano  '
        endif
        write(stdout,'(/)',advance='no')

        write(stdout,'(3a)',advance='no')
     +    ' -----------------------------------------------------',
     +    '-----------------------------------------------------',
     +    '-----------------------------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '------------------------------------'
        endif
        write(stdout,'(/)',advance='no')

      endif

      jnpos(1) = 0
      jnpos(2) = 0
      jnpos(7) = 0
      jnpos(8) = 0
      nreftot(1) = 0
      nreftot(2) = 0
      do i = 1, nbat

        jnpos(1) = jnpos(1) + inpos(1,i)
        jnpos(2) = jnpos(2) + inpos(2,i)
        if (l_ell) then
          jnpos(7) = jnpos(7) + inpos(5,i)
          jnpos(8) = jnpos(8) + inpos(6,i)
        end if

c       we have to double the anomalous number of reflections, since we
c       only store the batch were we get the second part of an I+/I-
c       pair (and nrefpos(2) counts I+ and I- separately).

        compl_ano_p = 0.0
        compl_ano_m = 0.0
        if(nrefpos(2) .ne. 0) then
          compl_ano_p = float(jnpos(2))/nrefpos(2)
          compl_ano_m = float(jnpos(4))/nrefpos(2)
        end if

        compl_ano_ell_p = 0.0
        compl_ano_ell_m = 0.0
        if (l_ell) then
          if(nrefpos(4) .ne. 0) then
            compl_ano_ell_p = float(jnpos( 8))/nrefpos(4)
            compl_ano_ell_m = float(jnpos(10))/nrefpos(4)
          endif
        end if

        rvaltmp = 0.0d0
        if (rval(4,1,1,i).gt.0.0d0) then
          rvaltmp(1,1) = min((rval(1,1,1,i)/rval(4,1,1,i)),99.9999)
          rvaltmp(2,1) = min((rval(2,1,1,i)/rval(4,1,1,i)),99.9999)
          rvaltmp(3,1) = min((rval(3,1,1,i)/rval(4,1,1,i)),99.9999)
        endif
        if (rval(4,2,1,i).gt.0.0d0) then
          rvaltmp(1,2) = min((rval(1,2,1,i)/rval(4,2,1,i)),99.9999)
          rvaltmp(2,2) = min((rval(2,2,1,i)/rval(4,2,1,i)),99.9999)
          rvaltmp(3,2) = min((rval(3,2,1,i)/rval(4,2,1,i)),99.9999)
        endif

        ivaltmp = 0.0d0
        if (ival(3,1,i).gt.0.0d0) then
c         <I>
          ivaltmp(1) = ival(1,1,i)/ival(3,1,i)
c         <I/sigI>
          ivaltmp(2) = ival(2,1,i)/ival(3,1,i)
        end if

        if ((iset_max-iset_min).gt.0) then
          write(stdout,
     +      101,
     +      advance='no')
     +      i2bat(i),
     +      i2set(i),
     +      nint(store(1,i)),
     +      store(2,i)/store(1,i),
     +      ivaltmp(1),ivaltmp(2),
     +      nint(rval(5,1,1,i)),
     +      rvaltmp(1,1),rvaltmp(1,2),
     +      rvaltmp(2,1),rvaltmp(2,2),
     +      rvaltmp(3,1),rvaltmp(3,2),
     +      float(jnpos(1))/nrefpos(1),
     +      compl_ano_p,
     +      float(jnpos(3))/nrefpos(1),
     +      compl_ano_m
          if (l_ell) then
            write(stdout,'(2(2x,2(1x,f7.4)))',advance='no')
     +        float(jnpos(7))/nrefpos(3),
     +        compl_ano_ell_p,
     +        float(jnpos(9))/nrefpos(3),
     +        compl_ano_ell_m
          end if
          write(stdout,'(/)',advance='no')
        else
          write(stdout,
     +      102,
     +      advance='no')
     +      i2bat(i),
     +      nint(store(1,i)),
     +      store(2,i)/store(1,i),
     +      ivaltmp(1),ivaltmp(2),
     +      nint(rval(5,1,1,i)),
     +      rvaltmp(1,1),rvaltmp(1,2),
     +      rvaltmp(2,1),rvaltmp(2,2),
     +      rvaltmp(3,1),rvaltmp(3,2),
     +      float(jnpos(1))/nrefpos(1),
     +      compl_ano_p,
     +      float(jnpos(3))/nrefpos(1),
     +      compl_ano_m
          if (l_ell) then
            write(stdout,'(2(2x,2(1x,f7.4)))',advance='no')
     +        float(jnpos(7))/nrefpos(3),
     +        compl_ano_ell_p,
     +        float(jnpos(9))/nrefpos(3),
     +        compl_ano_ell_m
          end if
          write(stdout,'(/)',advance='no')
        end if

        jnpos(3)  = jnpos(3)  - inpos(3,i)
        jnpos(4)  = jnpos(4)  - inpos(4,i)
        jnpos(9)  = jnpos(9)  - inpos(7,i)
        jnpos(10) = jnpos(10) - inpos(8,i)

        nreftot(1) = nreftot(1) + nint(store(1,i))
        nreftot(2) = nreftot(2) + nint(rval(5,1,1,i))

      end do

      compl_ano_p = 0.0
      if(nrefpos(2) .ne. 0) then
        compl_ano_p = float(jnpos(2))/nrefpos(2)
      end if
      if (l_ell) then
        compl_ano_ell_p = 0.0
        if(nrefpos(4) .ne. 0) then
          compl_ano_ell_p = float(jnpos(8))/nrefpos(4)
        end if
      end if

      if ((iset_max-iset_min).gt.0) then
        write(stdout,'(3a)',advance='no')
     +    ' -----------------------------------------------------',
     +    '-----------------------------------------------------',
     +    '----------------------------------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '------------------------------------'
        endif
        write(stdout,'(/)',advance='no')
        write(stdout,'(a8,6x,   i9,    30x,i9   ,56x,2(1x,f7.4))',
     +    advance='no')
     +    'Total',nreftot,
     +    float(jnpos(1))/nrefpos(1),
     +    compl_ano_p
        if (l_ell) then
          write(stdout,'(20x,2(1x,f7.4))',advance='no')
     +      float(jnpos(7))/nrefpos(3),
     +      compl_ano_ell_p
        end if
        write(stdout,'(//)',advance='no')
      else
        write(stdout,'(3a)',advance='no')
     +    ' -----------------------------------------------------',
     +    '-----------------------------------------------------',
     +    '-----------------------------------------'
        if (l_ell) then
          write(stdout,'(a)',advance='no')
     +      '------------------------------------'
        endif
        write(stdout,'(/)',advance='no')
        write(stdout,'(a8,1x,   i9,    30x,i9   ,56x,2(1x,f7.4))',
     +    advance='no')
     +    'Total',nreftot,
     +    float(jnpos(1))/nrefpos(1),
     +    compl_ano_p
        if (l_ell) then
          write(stdout,'(20x,2(1x,f7.4))',advance='no')
     +      float(jnpos(7))/nrefpos(3),
     +      compl_ano_ell_p
        end if
        write(stdout,'(//)',advance='no')
      endif

      deallocate(jref2iref,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate jref2iref array')
      deallocate(bat2i,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate bat2i array')
      deallocate(i2bat,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate i2bat array')
      deallocate(i2set,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate i2set array')
      deallocate(inpos,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate inpos array')
      deallocate(store,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate store array')
      deallocate(rval,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate rval array')
c
 101  format(1x,i8,1x,i4,1x,i8,1x,f8.3,1x,es11.3,1x,f8.3,1x,i8
     .  ,5(2x,2(1x,f7.4)))
 102  format(1x,i8      ,1x,i8,1x,f8.3,1x,es11.3,1x,f8.3,1x,i8
     .  ,5(2x,2(1x,f7.4)))
c
c----+=================================================================
c
      return
      end
