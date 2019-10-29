c     **********************************************************************
c
c      Copyright (c) 2010, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein, Claus Flensburg and Gerard Bricogne
c
c     **********************************************************************

c     compute quality metrics in shells/bins

      subroutine stats(refl,ncoltot,nref,nrun,irun)
c
c     Rmerge: http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_reflns_shell.Rmerge_I_obs.html
c             http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_reflns.Rmerge_F_all.html
c             http://www.iucr.org/__data/iucr/cifdic_html/2/cif_mm.dic/Ireflns.Rmerge_F_all.html
c
c       The value of Rmerge(I) for reflections classified as 'observed'
c       (see _reflns.observed_criterion) in a given shell.
c
c                      sum_i ( sum_j | Ij - <I> | )
c       Rmerge(I) = --------------------------------
c                      sum_i ( sum_j        <I>   )
c
c       Ij = the intensity of the jth observation of reflection i
c       <I> = the mean of the intensities of all observations of reflection i
c
c       sum_i is taken over all reflections
c       sum_j is taken over all observations of each reflection.
c
c     Rsym (J. Drenth, p. 289)
c
c                   sum_i ( sum_j | Ij - <I> | )
c       Rsym(I) = --------------------------------
c                   sum_i ( sum_j   Ij         )
c
c     Rpim (http://www.iucr.org/__data/iucr/cifdic_html/2/mmcif_pdbx.dic/Ireflns.pdbx_Rpim_I_all.html)
c
c      The precision-indicating merging R factor value Rpim,
c      for merging all intensities in this data set.
c   
c             sum_i [1/(Ni - 1)]^1/2 sum_j | Ij - <Ii> |
c      Rpim = --------------------------------------------------
c                           sum_i ( sum_j Ij )
c   
c      Ij   = the intensity of the jth observation of reflection i
c      <Ii> = the mean of the intensities of all observations
c               of reflection i
c      Ni   = the redundancy (the number of times reflection i
c               has been measured).
c   
c      sum_i is taken over all reflections
c      sum_j is taken over all observations of each reflection.
c   
c      Ref: Diederichs, K. & Karplus, P. A. (1997). Nature Struct.
c           Biol. 4, 269-275.
c           Weiss, M. S. & Hilgenfeld, R. (1997). J. Appl. Cryst.
c           30, 203-205.
c           Weiss, M. S. (2001). J. Appl. Cryst. 34, 130-135.
c   
c     Rmeas/Rrim
c
c       The redundancy-independent merging R factor value Rrim, also
c       denoted Rmeas, for merging all intensities in a given shell.
c
c                sum_i [Ni /( Ni - 1)]^(1/2) sum_j | Ij - <Ii> |
c       Rrim = --------------------------------------------------------
c                                 sum_i ( sum_j Ij )
c
c       Ij   = the intensity of the jth observation of reflection i
c       <Ii> = the mean of the intensities of all observations of
c                reflection i
c       Ni   = the redundancy (the number of times reflection i
c                has been measured).
c
c       sum_i is taken over all reflections
c       sum_j is taken over all observations of each reflection.
c
c       Ref: Diederichs, K. & Karplus, P. A. (1997). Nature Struct.
c            Biol. 4, 269-275.
c            Weiss, M. S. & Hilgenfeld, R. (1997). J. Appl. Cryst.
c            30, 203-205.
c            Weiss, M. S. (2001). J. Appl. Cryst. 34, 130-135.
c

      implicit none
c
      include 'args.inc'
      include 'params.inc'
      include 'lookup.inc'
      include 'info.inc'
      include 'units.inc'
      include 'runs.inc'
      include 'reso.inc'
      include 'cmpl.inc'
      include 'aniso.inc'
c
      integer ncoltot,nref(2),nrun,irun(nrun)
      double precision refl(ncoltot,nref(1))
c
c     local
c
      integer    maxrej
      parameter (maxrej=9)
c
c     pointers - to make it better readable
      integer iAll, iPlus, iMinus
      parameter (iAll=1, iPlus=2, iMinus=3)
      integer iGeneral, iAnomalous
      parameter (iGeneral=1, iAnomalous=2)
      integer iBoth, iHalfset1, iHalfset2
      parameter (iBoth=1, iHalfset1=2, iHalfset2=3)
      integer iMin, iMax
      parameter (iMin=1, iMax=2)
c
      integer hklold(3), hkl(3)
      integer i, j, ic, ic_old, ishell, ishellold, kshell, isym,
     +  nhkl(3),nreject(maxrej),nrefpos(6),
     +  isym_sum, isym_prod, ipm,jpm, jrun,iResoLimits, mshell,
     +  nrejmsg,
     +  ibat(3,maxmul), misym(3,maxmul),krun(3,maxmul), iCC(3),iref,
     +  ipm1, ipm2, ipm3, num_intensity_mrg(3,3),jref,k,
     +  ellval, i_elladd,
     +  isysabs, isysabs_old
      real epsi
      double precision intensity(3,maxmul), sum_intensity(3),
     +  mean_intensity,
     +  AbsDiff_intensity(3), sig_intensity(3,maxmul), r1, r2,
     +  r1max, r2min, r1stol, r2stol, r2min_active,
     +  cmpl_max, runs_reslim(2),
     +  r2LastPopulatedBin,
     +  dano_mrg, sano_mrg,
     +  dpm(3), w, intSet1, intSet2,
     +  intensity_mrg(3,3), sig_intensity_mrg(3,3),
     +  sum_intensity_over_sig(3), stol2res, weight_r(3), res2stol
c
c     binning by resolution shell
      double precision qesmin, qesmax, qesstp, qesmin3, qesmax3, qesstp3
c
      logical  l_write_rej
      external stol2res, res2stol
c
      double precision Rval(2,3,maxshell),Rval1(2,3,maxshell),
     +           Rval2(2,3,maxshell)
      integer NumRval(2,3,maxshell)
      double precision IsigI(maxshell)
      integer NumIsigI(maxshell)
      double precision IsigImrg(maxshell)
      integer NumIsigImrg(maxshell)
      integer NumUnique(2,maxshell),NumUniqueEll(2,maxshell)
      integer NumOverall(2,maxshell)
      integer NumReject(maxshell,maxrej)
      double precision Cmpltness(2,maxshell),CmpltnessEll(2,maxshell)
      integer NumCmpltness(2,maxshell),NumCmpltnessEll(2,maxshell)
      double precision Multiplicity(2,maxshell)
      integer NumMultiplicity(2,maxshell)
      double precision ResoLimits(2,maxshell)
      double precision ShellRes(2,maxshell)
      double precision ShellStol(2,maxshell)
      double precision SigAno(maxshell)
      integer NumSigAno(maxshell)
      integer NumCc(2,maxshell)
c
      double precision CC_Sx   (2,maxshell)
      double precision CC_Sy   (2,maxshell)
      double precision CC_Sxx  (2,maxshell)
      double precision CC_Syy  (2,maxshell)
      double precision CC_Sxy  (2,maxshell)
      integer          CC_n    (2,maxshell)
      double precision CC      (2,maxshell)
      double precision CC_print(2,maxshell)
c
      logical ActiveShell(maxshell)
c
      character*13 c_cmpl(2)
c
      character*25 pre
      character*34 reject_reason(maxrej)
      character*10 reject_line_fmt
      data reject_reason /
     +  'BGratio too large                 ',
     +  'PKratio too large (or XDS misfits)',
     +  'Negative < 5sigma                 ',
     +  'BG gradient too large             ',
     +  'Profile-fitted overloads          ',
     +  'Spots on edge                     ',
     +  'Partial                           ',
     +  'unused in R (one observation only)',
     +  'Systematic absence                '/
c
      integer istat, njref
      integer, allocatable :: jref2iref(:)
      integer, allocatable :: jref2ibin(:)
      integer n_jref2ibin
      double precision :: jref2reso(1:2,1:maxshell)
      character*10 restmp
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
      n_jref2ibin = 0

      allocate(jref2iref(1:(nref(2)+1)),stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to allocate jref2iref array')

      if ((nrun.eq.1).and.(irun(1).gt.0)) then
        if (L_BinningByEqualNumber_run) then
          allocate(jref2ibin(1:(nref(2)+1)),stat=istat)
        else
          allocate(jref2ibin(1:1),stat=istat)
        end if
        if (istat.ne.0)
     +    call error_exit('unable to allocate jref2ibin array')
        call refl_test_active(ncoltot,nref,refl,nrun,irun,
     +    jref2iref,njref,
     +    L_BinningByEqualNumber_run,
     +    N_BinningByEqualNumber_run,
     +      BinningByEqualNumber_run,
     +    nrefbin_run, jref2ibin, n_jref2ibin, jref2reso)
      else
        if (L_BinningByEqualNumber) then
          allocate(jref2ibin(1:(nref(2)+1)),stat=istat)
        else
          allocate(jref2ibin(1:1),stat=istat)
        end if
        if (istat.ne.0)
     +    call error_exit('unable to allocate jref2ibin array')
        call refl_test_active(ncoltot,nref,refl,nrun,irun,
     +    jref2iref,njref,
     +    L_BinningByEqualNumber    ,
     +    N_BinningByEqualNumber    ,
     +      BinningByEqualNumber    ,
     +    nrefbin    , jref2ibin, n_jref2ibin, jref2reso)
      end if

c
c     type of CC(1/2) and CC(anom) calculation:
c
c     as done in original Karplus&Diederichs 2012 paper (I think) and
c     also in AIMLESS source code (I think): based on the random number
c     of the first measurement, assign the remaining measurements
c     alternating to the two half-sets. Then calculate two (weighted?)
c     merged intensities and use those for the CC calculation.
c
      if ((nrun.eq.1).and.(irun(1).gt.0).and.(n_jref2ibin.eq.0)) then
c
c       take RUN-specific binning
c
        mshell  = nshell_run(irun(1))
        qesmin  = resmin_run(irun(1))
        qesmax  = resmax_run(irun(1))
        qesstp  = resstp_run(irun(1))
        qesmin3 = resmin3_run(irun(1))
        qesmax3 = resmax3_run(irun(1))
        qesstp3 = resstp3_run(irun(1))

      else

        if (n_jref2ibin.gt.0) then
c
c         take pre-defined binning (either from user or via automatic
c         binning by equal numbers: see s/r refl_test_active)
c
          mshell  = n_jref2ibin
        else
c
c         take overall binning
c
          mshell  = nshell
          qesmin  = resmin
          qesmax  = resmax
          qesstp  = resstp
          qesmin3 = resmin3
          qesmax3 = resmax3
          qesstp3 = resstp3
        end if

      end if

      write(stdout,'(a,i4,a,/)') ' We will be using a total of ',
     +  mshell,' bins for statistics as a function of resolution'

      if (n_jref2ibin.gt.0) then
        write(stdout,'(a)',advance='no')
     +    ' Resolution bins defined as : '
        write(restmp,'(f10.5)') jref2reso(2,1)
        write(stdout,'(a)',advance='no') trim(adjustl(restmp))

        do i = 2, n_jref2ibin - 1
          write(restmp,'(f10.5)') jref2reso(2,i)
          write(stdout,'(a)',advance='no') ','//trim(adjustl(restmp))
        end do

        write(stdout,'(a,/)') ' '
      end if
c
c     ----------------------------------------------------------
c     initialisation
      Rval = 0.0d0
      Rval1 = 0.0d0
      Rval2 = 0.0d0
      NumRval = 0
      do kshell = 1, mshell + 1
        IsigI(kshell) = 0.0d0
        NumIsigI(kshell) = 0
        IsigImrg(kshell) = 0.0d0
        NumIsigImrg(kshell) = 0
        SigAno(kshell) = 0.0d0
        NumSigAno(kshell) = 0
        do j = iGeneral, iAnomalous
          NumCc(j,kshell) = 0
          NumUnique(j,kshell) = 0
          NumUniqueEll(j,kshell) = 0
          NumOverall(j,kshell) = 0
          Cmpltness(j,kshell) = 0.0d0
          CmpltnessEll(j,kshell) = 0.0d0
          NumCmpltness(j,kshell) = 0
          NumCmpltnessEll(j,kshell) = 0
          ResoLimits(j,kshell) = 0.
          Multiplicity(j,kshell) = 0.0d0
          NumMultiplicity(j,kshell) = 0
          CC_Sx (j,kshell) = 0.0d0
          CC_Sy (j,kshell) = 0.0d0
          CC_Sxx(j,kshell) = 0.0d0
          CC_Syy(j,kshell) = 0.0d0
          CC_Sxy(j,kshell) = 0.0d0
          CC_n  (j,kshell) = 0
        end do
        do j = 1, maxrej
          NumReject(kshell,j) = 0
        end do
      end do
c
      do j = 1, 3
        hklold(j) = -99999
      end do
      do ipm = iAll, iMinus
        nhkl(ipm) = 0
        do j = iBoth, iHalfset2
          num_intensity_mrg(j,ipm) = 0
          intensity_mrg(j,ipm) = 0.0d0
          sig_intensity_mrg(j,ipm) = 0.0d0
        end do
      end do
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
      do i = 1, maxrej
        nreject(i) = 0
      end do
      do ipm = iAll, iMinus
        sum_intensity(ipm)= 0.0d0
      end do
c
c     ----------------------------------------------------------
c
      ellval = 0
c
c     -------------------------------------------------------------------------
c     loop over all reflections
c     -------------------------------------------------------------------------
c
      njref = njref + 1
      jref2iref(njref) = nref(2) + 1
      do jref = 1, njref
        iref = jref2iref(jref)
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

          call centr(hkl,ic)
c
c         check for systematic absences:
          isysabs = 0
          call epslon(hkl,epsi,isysabs)
c
c         which resolution shell are we in
c
          if (n_jref2ibin.gt.0) then
            ishell = jref2ibin(jref)
          else
            if (lshell3) then
              ishell = int((refl(icol(i_reso3),iref)-qesmin3)/qesstp3)+1
              ishell = max(1,min(ishell,mshell))
            else
              ishell = int( (refl(icol(i_resol),iref)-qesmin)/qesstp )+1
              ishell = max(1,min(ishell,mshell))
            end if
          end if

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
          if ((nhkl(iAll).gt.0).and.(isysabs_old.gt.0)) then
            nreject(9) = nreject(9) + 1

          else if (nhkl(iAll).gt.0) then
c
c           is this a centric one or not?
            call centr(hklold,ic_old)

c           ellval = 1 (inside)
c                    0 (outside)
            if (l_ell) call ellchk(hklold,ellval)
c
c           1 = all observations = iAll
c           2 = I+ observations  = iPlus
c           3 = I- observations  = iMinus
c
c           get merged intensities (and sigmas) for all (full and
c           halfsets):
            do ipm = iAll, iMinus
              do j = iBoth, iHalfset2
                if (num_intensity_mrg(j,ipm).gt.0) then
                  intensity_mrg(j,ipm) = 
     +            intensity_mrg(j,ipm) / sig_intensity_mrg(j,ipm)
                  sig_intensity_mrg(j,ipm) =
     +            sqrt(1.d0/sig_intensity_mrg(j,ipm))
                else
                  intensity_mrg(j,ipm) = 0.0d0
                  sig_intensity_mrg(j,ipm) = 0.0d0
                end if
              end do
            end do

c           store in current shell and TOTAL:
            do kshell = ishellold,(mshell+1),(mshell+1-ishellold)

c             multiplicity
              NumUnique      (iGeneral,kshell) =
     +        NumUnique      (iGeneral,kshell) + 1
              if (l_ell.and.(ellval.eq.1)) then
                NumUniqueEll   (iGeneral,kshell) =
     +          NumUniqueEll   (iGeneral,kshell) + 1
              end if

              NumOverall     (iGeneral,kshell) =
     +        NumOverall     (iGeneral,kshell) + nhkl(iAll)
              NumMultiplicity(iGeneral,kshell) =
     +        NumMultiplicity(iGeneral,kshell) + 1
              Multiplicity   (iGeneral,kshell) =
     +        Multiplicity   (iGeneral,kshell) + nhkl(iAll)

c             I/sigI
              IsigI(kshell) = IsigI(kshell) + 
     +          sum_intensity_over_sig(iAll)
              NumIsigI(kshell) = NumIsigI(kshell) + nhkl(iAll)

c             Imrg/sigImrg
              IsigImrg(kshell)      =
     +        IsigImrg(kshell)      +
     +        intensity_mrg(iBoth,iAll)/sig_intensity_mrg(iBoth,iAll)
              NumIsigImrg(kshell)   =
     +        NumIsigImrg(kshell)   + 1

            end do

c
c           CC(1/2)
c
            if ((num_intensity_mrg(iHalfset1,iAll).gt.0).and.
     +          (num_intensity_mrg(iHalfset2,iAll).gt.0)     ) then

              intSet1 = intensity_mrg(iHalfset1,iAll)
              intSet2 = intensity_mrg(iHalfset2,iAll)

              do kshell = ishellold, (mshell+1), (mshell+1-ishellold)

                CC_Sx (iGeneral,kshell) =
     +          CC_Sx (iGeneral,kshell) + intSet1
                CC_Sy (iGeneral,kshell) =
     +          CC_Sy (iGeneral,kshell) + intSet2
                CC_Sxx(iGeneral,kshell) =
     +          CC_Sxx(iGeneral,kshell) + intSet1**2
                CC_Syy(iGeneral,kshell) =
     +          CC_Syy(iGeneral,kshell) + intSet2**2
                CC_Sxy(iGeneral,kshell) =
     +          CC_Sxy(iGeneral,kshell) + intSet1*intSet2
                CC_n  (iGeneral,kshell) =
     +          CC_n  (iGeneral,kshell) + 1
              end do
            end if

c           anomalous data:
            if ( (ic_old.eq.0)      .and.
     +           (nhkl(iPlus).gt.0) .and.
     +           (nhkl(iMinus).gt.0)     ) then

c
c             get anomalous difference (of merged intensities) and its
c             sigma:

              dano_mrg = intensity_mrg(iBoth,iPlus) -
     +                   intensity_mrg(iBoth,iMinus)
              sano_mrg = sqrt(sig_intensity_mrg(iBoth,iPlus )**2 +
     +                        sig_intensity_mrg(iBoth,iMinus)**2  )

c             store in current shell and TOTAL:
              do kshell = ishellold,(mshell+1),(mshell+1-ishellold)

c               multiplicity - and do we have both measurements?
                if ((isym_sum.gt.0).and.(isym_prod.eq.0)) then
                  NumUnique      (iAnomalous,kshell) =
     +            NumUnique      (iAnomalous,kshell) + 1
c
                  if (l_ell.and.(ellval.eq.1)) then
                    NumUniqueEll   (iAnomalous,kshell) =
     +              NumUniqueEll   (iAnomalous,kshell) + 1
                  end if
c
                  NumOverall     (iAnomalous,kshell) =
     +            NumOverall     (iAnomalous,kshell) + nhkl(iAll)

c                 original method:
                  NumMultiplicity(iAnomalous,kshell) =
     +            NumMultiplicity(iAnomalous,kshell) + 2
                  Multiplicity   (iAnomalous,kshell) =
     +            Multiplicity   (iAnomalous,kshell) + nhkl(iAll)

                end if

c               DANO|/SANO
                SigAno(kshell) =
     +          SigAno(kshell) + abs(dano_mrg)/sano_mrg
                NumSigAno(kshell) =
     +          NumSigAno(kshell) + 1
c
              end do
c
c             according to AIMLESS source code, the anomalous
c             correlation between half-sets is calculated as follows:
c
c             * divide all I+ measurements into two random half-sets
c               and get average I for each:
c                 I+1  I+2
c
c             * do the same for the I- measurements:
c                 I-1  I-2
c
c             * use CC( (I+1 - I-1), (I+2 - I-2) )
c
              if ((num_intensity_mrg(iHalfset1,iPlus ).gt.0).and.
     +            (num_intensity_mrg(iHalfset2,iPlus ).gt.0).and.
     +            (num_intensity_mrg(iHalfset1,iMinus).gt.0).and.
     +            (num_intensity_mrg(iHalfset2,iMinus).gt.0)     ) then

                dpm(iHalfset1) = intensity_mrg(iHalfset1,iPlus) - 
     +                           intensity_mrg(iHalfset1,iMinus)
                dpm(iHalfset2) = intensity_mrg(iHalfset2,iPlus) - 
     +                           intensity_mrg(iHalfset2,iMinus)
c
                do kshell = ishellold, (mshell+1), (mshell+1-ishellold)
                  CC_Sx (iAnomalous,kshell) =
     +            CC_Sx (iAnomalous,kshell) + dpm(iHalfset1)
                  CC_Sxx(iAnomalous,kshell) =
     +            CC_Sxx(iAnomalous,kshell) + dpm(iHalfset1)**2
                  CC_Sy (iAnomalous,kshell) =
     +            CC_Sy (iAnomalous,kshell) + dpm(iHalfset2)
                  CC_Syy(iAnomalous,kshell) =
     +            CC_Syy(iAnomalous,kshell) + dpm(iHalfset2)**2
                  CC_Sxy(iAnomalous,kshell) =
     +            CC_Sxy(iAnomalous,kshell) +
     +              dpm(iHalfset1)*dpm(iHalfset2)
                  CC_n  (iAnomalous,kshell) =
     +            CC_n  (iAnomalous,kshell) + 1
                end do

              end if

            end if
c
c           ============================================================
c           end of simple processing for completely fetched, previous
c           reflection (ie. anything we can do with stats accumulated on
c           -the-fly while reading all measurements). Now we want to do
c           something with the processed data for this reflection where
c           we need to loop over the measurements a second time.
c           ============================================================

            if (nhkl(iAll).gt.1) then
c
c             1 = normal statistics, treating I+ and I- identical
c             2 = treating I+ and I- separate
c
              do ipm = iAll, iMinus
c               we can only compute R-values if we have more than one
c               measurement:
                if (nhkl(ipm).gt.1) then

c                 collaps both of I+ and I- into one set of stats - with
c                 the only difference that we treat them separately:
                  jpm = min(ipm,2)
c
c                 get mean intensity for that class (iAll, iMinus or
c                 iPlus)
                  mean_intensity = sum_intensity(ipm) / nhkl(ipm)

c

c                 Rmerge
                  weight_r(1) =
     +              1.0
c                 Rmeas
                  weight_r(2) =
     +              sqrt(float(nhkl(ipm))/(nhkl(ipm)-1))
c                 Rpim
                  weight_r(3) =
     +              sqrt(             1.0/(nhkl(ipm)-1))
c
c                 loop over all measurements for that class:
                  do j = 1, nhkl(ipm)

c                   Rmerge Rmeas Rpim
                    do k = 1, 3

                      if (calc_type_r(2,k).eq.1) then
                        AbsDiff_intensity(k) = weight_r(k) *
     +                    abs(intensity(ipm,j) - mean_intensity)
                      else
                        AbsDiff_intensity(k) = weight_r(k) *
     +                    abs(intensity(ipm,j) - intensity_mrg(1,ipm))
                      end if
c
c                     sums for Rmerge Rmeas Rpim
c
                      do kshell = ishellold, (mshell+1),
     +                                       (mshell+1-ishellold)
                        Rval1  (jpm,k,kshell) =
     +                  Rval1  (jpm,k,kshell) + AbsDiff_intensity(k)

                        if (calc_type_r(1,k).eq.1) then
                          Rval2  (jpm,k,kshell) =
     +                    Rval2  (jpm,k,kshell) + intensity(ipm,j)
                        else
                          Rval2  (jpm,k,kshell) =
     +                    Rval2  (jpm,k,kshell) + intensity_mrg(1,ipm)
                        end if
                        NumRval(jpm,k,kshell) =
     +                  NumRval(jpm,k,kshell) + 1
                      end do

                    end do

                  end do

                end if

              end do

            else if (nhkl(iAll).eq.1) then
c             only one measurement for this reflection: can't be used
c             for R-values:
              if (i.le.nref(2)) nreject(8) = nreject(8) + 1
            end if
c
          end if

c         ============================================================
c         finished processing the previous one
c         ============================================================

          if (i.le.nref(2)) then
            do j = 1, maxrej
              do kshell = ishell, (mshell+1), (mshell+1-ishell)
                NumReject (kshell,j) =
     +          NumReject (kshell,j) + nreject(j)
              end do
              nreject(j) = 0
            end do
          end if

c         we are after the last reflection - exit loop
          if (i.gt.nref(2)) then
            goto 10
          end if

c         clear-out some arrays - to be ready for next unique
c         reflection:
          do ipm = iAll, iMinus
            nhkl(ipm) = 0
            sum_intensity(ipm) = 0.0d0
            sum_intensity_over_sig(ipm) = 0.0d0
            iCC(ipm) = 0
            do j = iBoth, iHalfset2
              intensity_mrg(j,ipm) = 0.0d0
              sig_intensity_mrg(j,ipm) = 0.0d0
              num_intensity_mrg(j,ipm) = 0
            end do
          end do
          isym_sum  = 0
          isym_prod = 1

        end if
c
c       ============================================================
c       start processing the new one
c       ============================================================
c
c       we're already past the last reflection:
        if ( (hkl(1).eq.0).and.
     +       (hkl(2).eq.0).and.
     +       (hkl(3).eq.0)     ) then
          goto 10
        end if

c
c       get ISYM:
c
c       odd  values are I(+)
c       even values are I(-)
c
c       each friedel pair is consecutive
c       i.e. isym=(h,k,l) and isym+1=(-h,-k,-l)
c
        isym      = nint(refl(icol(i_misym),iref))
        isym_sum  = isym_sum  + mod(isym,2)
        isym_prod = isym_prod * mod(isym,2)

        ellval = nint(refl(icol(i_ell),iref))
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
c
            if (nhkl(ipm).gt.maxmul) then
              write(stdout,'(/,a,i6,/)')
     +          ' maximum allowed multiplicity = ',maxmul
              call error_exit('multiplicity too large')
            end if

c           fetch intensity and sigma
            intensity    (ipm,nhkl(ipm)) = refl(icol(i_i   ),iref)
            sig_intensity(ipm,nhkl(ipm)) = refl(icol(i_sigi),iref)

c           sum intensities
            sum_intensity(ipm) =
     +      sum_intensity(ipm) + intensity(ipm,nhkl(ipm))

c           sum intensity/sigma
            sum_intensity_over_sig(ipm) =
     +      sum_intensity_over_sig(ipm) +
     +          intensity(ipm,nhkl(ipm))/sig_intensity(ipm,nhkl(ipm))

c           get batch number ...
            ibat (ipm,nhkl(ipm)) = nint(refl(icol(i_batch),iref))

c           and m/isym:
            misym(ipm,nhkl(ipm)) = nint(refl(icol(i_misym),iref))
            krun (ipm,nhkl(ipm)) = nint(refl(icol(i_run),iref))
c
c           assign to halfsets:
c
c           this way we take the random number for the first
c           measurement and then alternate:
            if (iCC(ipm).eq.0) then
c             iBOTH=1 iHalfset1=2 iHalfset2=3:
              iCC(ipm) = int((refl(icol(i_rand),iref))*2.0)+2
            else
              iCC(ipm) = mod((iCC(ipm)-1),2) + 2
            end if
c
c           sum for merged intensity/sigma - variance-weighted:
            w = 1.d0/(sig_intensity(ipm,nhkl(ipm))**2)

c           for CC(1/2):
                intensity_mrg(iBoth,ipm) = 
     +          intensity_mrg(iBoth,ipm) + w*intensity(ipm,nhkl(ipm))
            sig_intensity_mrg(iBoth,ipm) = 
     +      sig_intensity_mrg(iBoth,ipm) + w
            num_intensity_mrg(iBoth,ipm) = 
     +      num_intensity_mrg(iBoth,ipm) + 1

c           for CC(anom):
                intensity_mrg(iCC(ipm),ipm) = 
     +          intensity_mrg(iCC(ipm),ipm) + w*intensity(ipm,nhkl(ipm))
            sig_intensity_mrg(iCC(ipm),ipm) = 
     +      sig_intensity_mrg(iCC(ipm),ipm) + w
            num_intensity_mrg(iCC(ipm),ipm) = 
     +      num_intensity_mrg(iCC(ipm),ipm) + 1

          end do

        do j = 1, 3
          hklold(j) = hkl(j)
        end do
        ishellold     = ishell
        isysabs_old   = isysabs
c        
 10     continue
c
      end do
c
      if (NumUnique(1,(mshell+1)).eq.0) goto 999
c
c     ---------------------------------------------------------------------
c     finalise stats
c
      do jrun = 1, nrun
        if (jrun.eq.1) then
          runs_reslim(1) = run_reslim(1,irun(jrun))
          runs_reslim(2) = run_reslim(2,irun(jrun))
        else
          runs_reslim(1) = max(runs_reslim(1),run_reslim(1,irun(jrun)))
          runs_reslim(2) = min(runs_reslim(2),run_reslim(2,irun(jrun)))
        end if
      end do
c
      cmpl_max = -1.0d0
c
      do kshell = 1, mshell + 1
        do j = 1, 2
          do k = 1, 3
            if ((NumRval(j,k,kshell)   .gt.0    ).and.
     +          (abs(Rval2(j,k,kshell)).gt.0.001)     ) then
              Rval(j,k,kshell) = Rval1(j,k,kshell)/Rval2(j,k,kshell)
            end if
          end do
        end do
        if (NumIsigI(kshell).gt.0) then
          IsigI(kshell) = IsigI(kshell)/NumIsigI(kshell)
        end if
        if (NumIsigImrg(kshell).gt.0) then
          IsigImrg(kshell) = IsigImrg(kshell)/NumIsigImrg(kshell)
        end if
        if (NumSigAno(kshell).gt.0) then
          SigAno(kshell) = SigAno(kshell)/NumSigAno(kshell)
        end if

        do j = iGeneral, iAnomalous
c         default to 1.01 - so that empty bins are ignored
          CC(j,kshell) = 1.01d0
          CC_print(j,kshell) = 0.0d0
c         need to have at least 3 pairs for doing a correlation
          if (CC_n(j,kshell).gt.2) then
            CC(j,kshell) = (
     +        CC_Sxy(j,kshell) - (
     +          CC_Sx(j,kshell)*CC_Sy(j,kshell)/CC_n(j,kshell)
     +        ) ) /
     +        sqrt(
     +        (CC_Sxx(j,kshell)-(CC_Sx(j,kshell)**2/CC_n(j,kshell)))*
     +        (CC_Syy(j,kshell)-(CC_Sy(j,kshell)**2/CC_n(j,kshell)))
     +        )
            CC_print(j,kshell) = CC(j,kshell)
            NumCc(j,kshell) = NumCc(j,kshell) + CC_n(j,kshell)
c
          end if
        end do
c
c       get resolution limits
c
        if (n_jref2ibin.gt.0) then
          r1stol = res2stol(jref2reso(1,kshell))
          r2stol = res2stol(jref2reso(2,kshell))
        else
          if (lshell3) then
            r1stol = qesmin3 + (kshell-1)*qesstp3
            r2stol = r1stol + qesstp3
            if (kshell.ge.mshell) r2stol = qesmax3
            if (kshell.gt.mshell) r1stol = qesmin3
          else
            r1stol = qesmin + (kshell-1)*qesstp
            r2stol = r1stol + qesstp
            if (kshell.ge.mshell) r2stol = qesmax
            if (kshell.gt.mshell) r1stol = qesmin
          end if
        end if
        r1 = stol2res(r1stol)
        r2 = stol2res(r2stol)
c
        if (kshell.le.mshell) then

          ShellRes(1,kshell)  = r1
          ShellRes(2,kshell)  = r2
          ShellStol(1,kshell) = r1stol
          ShellStol(2,kshell) = r2stol

        end if
c
        if (n_jref2ibin.gt.0) then
          ResoLimits(1,kshell) = jref2reso(1,kshell)
          ResoLimits(2,kshell) = jref2reso(2,kshell)
        else
          ResoLimits(1,kshell) =
     +      max(min(r1,runs_reslim(1)),runs_reslim(2))
          ResoLimits(2,kshell) =
     +      min(max(r2,runs_reslim(2)),runs_reslim(1))
        end if
c         
c       limit to resolution limits of run(s)
c
        r1 = min(runs_reslim(1),r1)
        if (r1.gt.runs_reslim(2)) then
          r2=max(runs_reslim(2),r2)
c
c         calculate number of possible reflections in that shell
c
          if (kshell.eq.(mshell+1)) then
            call cmpl(r1,r2LastPopulatedBin,1,nrefpos)
          else
            if (NumOverall(iGeneral,kshell).gt.0) then
              r2LastPopulatedBin = r2
            end if
            call cmpl(r1,r2,0,nrefpos)
          end if
c
c         calculate completeness
c
          do j = iGeneral, iAnomalous
            NumCmpltness(j,kshell) = nrefpos(j)
            if (NumCmpltness(j,kshell).ne.0) then
              Cmpltness(j,kshell) =
     +          float(NumUnique(j,kshell))/NumCmpltness(j,kshell)
            end if
            if (l_ell) then
              NumCmpltnessEll(j,kshell) = nrefpos((2+j))
              if (NumCmpltnessEll(j,kshell).ne.0) then
                CmpltnessEll(j,kshell) =
     +          float(NumUniqueEll(j,kshell))/NumCmpltnessEll(j,kshell)
              end if
            end if
          end do
          if (kshell.lt.(mshell+1)) then
            cmpl_max = max(cmpl_max,Cmpltness(iGeneral,kshell))
          end if
        end if
c
c       multiplicity
        do j = iGeneral, iAnomalous
          if (NumMultiplicity(j,kshell).gt.0) then
            Multiplicity(j,kshell) =
     +      Multiplicity(j,kshell)/NumMultiplicity(j,kshell)
          end if
        end do
      end do
c
c     report some stats
c
      l_write_rej = .false.
      if (iwide_format.eq.1) then
        nrejmsg = 0
        do i = 1, maxrej
          if (NumReject((mshell+1),i).gt.0) then
            if (.not.l_write_rej) then
              write(stdout,'(/,a)') ' Rejection   reason'
              write(stdout,'(a)')
     +          ' ------------------------------------------------'
              l_write_rej = .true.
            endif
            write(stdout,'(a,i1,2a)')
     +        '  [',i,']        ',reject_reason(i)
            nrejmsg = nrejmsg + 1
          end if
        end do
        if (l_write_rej) write(stdout,'(/)')
      end if

      i_elladd = 0
      if (l_ell) then
        i_elladd = 17
      end if

      if (iwide_format.eq.1) then
c                                                                                                   1
c         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8
c123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
c                                                  Rmerge          Rmeas           Rpim                        Completeness   Multiplicity                                  Rejections
c                                             --------------  --------------  --------------                  --------------  ------------                                 -----------
c            Resolution      #uniq     #Rfac     all     ano     all     ano     all     ano   #Isig  I/sigI     all     ano    all    ano CC(1/2)  #CcAno CC(ano) #SigAno      [8] 
c -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        write(reject_line_fmt,'(a,i3,a)') '(',nrejmsg*8,'a)'
        if (l_write_rej) then
          write(stdout,101,advance='no') ('-',i=1,181)
          if (i_elladd.gt.0) then
            do i = 1, i_elladd
              write(stdout,'(a1)',advance='no') '-'
            end do
          endif
          write(stdout,reject_line_fmt) ('-',i=1,(nrejmsg*8))
          write(stdout,102,advance='no')
     +      'Rmerge','Rmeas','Rpim',c_cmpl(1),
     +                   'Multiplicity'
          if (l_ell) then
            write(stdout,'(a)',advance='no') '   '//c_cmpl(2)//'  '
          end if
          write(stdout,'(a11)') 'Rejections'
          write(stdout,103,advance='no')
     +      '--------------','--------------',
     +      '--------------','--------------',
     +      '------------'
          if (l_ell) then
            write(stdout,'(a)',advance='no') '  --------------- '
          end if
          write(stdout,'(a11)') '-----------'

        else
c ------------------------------------------------------------------------------------------------------------------------------------------------
c                                                                                      Completeness   Multiplicity
c                                                                                     --------------  ------------
c            Resolution      #uniq     #Rfac  Rmerge   Rmeas    Rpim   #Isig  I/sigI     all     ano    all    ano CC(1/2)  #CCAno CC(ano)  SigAno
c ------------------------------------------------------------------------------------------------------------------------------------------------
          write(stdout,201,advance='no') ('-',i=1,168)
          if (i_elladd.gt.0) then
            do i = 1, i_elladd
              write(stdout,'(a1)',advance='no') '-'
            end do
          end if
          write(stdout,'(/)',advance='no')
          write(stdout,202,advance='no')
     +      'Rmerge','Rmeas','Rpim',c_cmpl(1),
     +      'Multiplicity'
          if (l_ell) then
            write(stdout,'(a)') '   '//c_cmpl(2)//' '
          else
            write(stdout,'(/)',advance='no')
          end if
          write(stdout,203,advance='no')
     +      '--------------','--------------',
     +      '--------------','--------------',
     +      '------------'
          if (l_ell) then
            write(stdout,'(a)') '  --------------- '
          else
            write(stdout,'(/)',advance='no')
          end if
        end if
        write(stdout,104,advance='no')
     +    'Resolution',
     +    '#uniq',
     +    '#Rfac',
     +    'all','ano',
     +    'all','ano',
     +    'all','ano',
     +    '#Isig','I/sigI',
     +    'all','ano',
     +    'all','ano',
     +    'CC(1/2)',
     +    '#CcAno','CC(ano)','#SigAno','SigAno'
        if (l_ell) write(stdout,'(a)',advance='no') '  all     ano   '
        if (l_write_rej) then
          do j = 1, maxrej
            if (NumReject((mshell+1),j).gt.0) then
              write(stdout,'(a,i1,a)',advance='no') '    [',j,'] '
            end if
          end do
        end if
        write(stdout,'(/)',advance='no')
        if (l_write_rej) then
          write(stdout,101,advance='no') ('-',i=1,181)
          if (i_elladd.gt.0) then
            do i = 1, i_elladd
              write(stdout,'(a1)',advance='no') '-'
            end do
          end if
          write(stdout,reject_line_fmt) ('-',i=1,(nrejmsg*8))
        else
          write(stdout,201,advance='no') ('-',i=1,176)
          if (i_elladd.gt.0) then
            do i = 1, i_elladd
              write(stdout,'(a1)',advance='no') '-'
            end do
          end if
          write(stdout,'(/)',advance='no')
        end if

c       narrow format:
      else
c                                                                                                   1
c         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8
c123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
c ------------------------------------------------------------------------------------------------------------------------------------------------
c                                                                                      Completeness   Multiplicity                                
c                                                                                     --------------  ------------                                
c            Resolution        #uniq   #Rfac   Rmerge  Rmeas    Rpim   #Isig  I/sigI     all     ano    all    ano CC(1/2)  #CcAno CC(ano) SigAno
c ------------------------------------------------------------------------------------------------------------------------------------------------
c
        write(stdout,301,advance='no') ('-',i=1,144)
        if (i_elladd.gt.0) then
          do i = 1, i_elladd
            write(stdout,'(a1)',advance='no') '-'
          end do
        end if
        write(stdout,'(/)',advance='no')
        write(stdout,302,advance='no')
     +    c_cmpl(1),
     +    'Multiplicity'
        if (l_ell) then
          write(stdout,'(32x,a)') '   '//c_cmpl(2)//' '
        else
          write(stdout,'(/)',advance='no')
        end if
        write(stdout,303,advance='no')
     +    '--------------',
     +    '------------'
        if (l_ell) then
          write(stdout,'(32x,a)') '  --------------- '
        else
          write(stdout,'(/)',advance='no')
        end if
        write(stdout,304,advance='no')
     +    'Resolution',
     +    '#uniq',
     +    '#Rfac',
     +    'Rmerge',
     +    'Rmeas',
     +    'Rpim',
     +    '#Isig','I/sigI',
     +    'all','ano',
     +    'all','ano',
     +    'CC(1/2)',
     +    '#CCAno','CC(ano)','SigAno'
        if (l_ell) write(stdout,'(a)',advance='no') '  all     ano'
        write(stdout,'(/)',advance='no')
        write(stdout,301,advance='no') ('-',i=1,144)
        if (i_elladd.gt.0) then
          do i = 1, i_elladd
            write(stdout,'(a1)',advance='no') '-'
          end do
        end if
        write(stdout,'(/)',advance='no')
c
      end if
c
c
      r1max = -9999.9
      r2min =  9999.9

c     avoid writing empty shells towards the end:
      ActiveShell = .true.
      r2min_active = r2min
      do kshell = mshell, 1, -1
        if (NumOverall(1,kshell).gt.0) then
          exit
        end if
        if (NumOverall(1,kshell).eq.0) ActiveShell(kshell) = .false.
      end do
c
      do kshell = 1, mshell
c
        if (.not.ActiveShell(kshell)) cycle
c
c       get resolution limits
c
        r1 = ShellRes(1,kshell)
        r2 = ShellRes(2,kshell)
        r1stol = ShellStol(1,kshell)
        r2stol = ShellStol(2,kshell)
c
c       limit to resolution limits of run(s)
c
        r1 = min(runs_reslim(1),r1)
        if (r1.gt.runs_reslim(2)) then
c
c          if (NumOverall(1,kshell).gt.0) then
c
          r2    = max(runs_reslim(2),r2)
          r1max = max(r1,r1max)
          r2min = min(r2,r2min)
          r2min_active = r2min
c
          if (NumOverall(1,kshell).gt.0) iResoLimits = kshell
          write(pre,'(8x,f7.3,a3,f7.3)') r1,' - ',r2
c
          if (iwide_format.eq.1) then
            if (NumOverall(1,kshell).gt.0) then
              write(stdout,1001,advance='no') pre,
     +          NumUnique(1,kshell),
     +          NumRval(1,1,kshell),
     +          min(99.999,max(0.,Rval(iGeneral,1,kshell))),
     +          min(99.999,max(0.,Rval(iAnomalous,1,kshell))),
     +          min(99.999,max(0.,Rval(iGeneral,2,kshell))),
     +          min(99.999,max(0.,Rval(iAnomalous,2,kshell))),
     +          min(99.999,max(0.,Rval(iGeneral,3,kshell))),
     +          min(99.999,max(0.,Rval(iAnomalous,3,kshell))),
     +          NumIsigImrg(kshell),IsigImrg(kshell),
     +          min(1.0,Cmpltness(iGeneral,kshell)),
     +          min(1.0,Cmpltness(iAnomalous,kshell)),
     +          Multiplicity(iGeneral,kshell),
     +          Multiplicity(iAnomalous,kshell),
     +          CC_print(iGeneral,kshell),
     +          NumCc(iAnomalous,kshell),CC_print(iAnomalous,kshell),
     +          NumSigAno(kshell),SigAno(kshell)
              if (l_ell) then
                write(stdout,'(2f8.4)',advance='no')
     +            min(1.0,CmpltnessEll(iGeneral,kshell)),
     +            min(1.0,CmpltnessEll(iAnomalous,kshell))
              end if
              if(l_write_rej) then
                do j = 1, maxrej
                  if (NumReject((mshell+1),j).gt.0) then
                    write(stdout,'(i8)',advance='no')
     +                NumReject(kshell,j)
                  end if
                end do
              end if
              write(stdout,'(/)',advance='no')
            else
              write(stdout,'(a25,1x,i7)') pre,0
            end if

c           narrow format
          else
            if (NumOverall(1,kshell).gt.0) then
              write(stdout,2001,advance='no') pre,
     +          NumUnique(1,kshell),
     +          NumRval(1,1,kshell),
     +          min(99.999,max(0.,Rval(iGeneral,1,kshell))),
     +          min(99.999,max(0.,Rval(iGeneral,2,kshell))),
     +          min(99.999,max(0.,Rval(iGeneral,3,kshell))),
     +          NumIsigImrg(kshell),IsigImrg(kshell),
     +          min(1.0,Cmpltness(iGeneral,kshell)),
     +          min(1.0,Cmpltness(iAnomalous,kshell)),
     +          Multiplicity(iGeneral,kshell),
     +          Multiplicity(iAnomalous,kshell),
     +          CC_print(iGeneral,kshell),
     +          NumCc(iAnomalous,kshell),CC_print(iAnomalous,kshell),
     +          SigAno(kshell)
              if (l_ell) then
                write(stdout,'(2f8.4)',advance='no')
     +            min(1.0,CmpltnessEll(iGeneral,kshell)),
     +            min(1.0,CmpltnessEll(iAnomalous,kshell))
              end if
              write(stdout,'(/)',advance='no')
            else
              write(stdout,'(a25,1x,i7)') pre,0
            end if
          end if
c
        end if
c
      end do
c
      if (iwide_format.eq.1) then
        if (l_write_rej) then
          write(stdout,101,advance='no') ('-',i=1,173)
          write(stdout,reject_line_fmt,advance='no')
     +      ('-',i=1,(nrejmsg*8))
        else
          write(stdout,201,advance='no') ('-',i=1,168)
        end if

c       narrow format
      else
        write(stdout,301,advance='no') ('-',i=1,144)
      end if
      if (i_elladd.gt.0) then
        do i = 1, i_elladd
          write(stdout,'(a1)',advance='no') '-'
        end do
      end if
      write(stdout,'(/)',advance='no')
c
c     Total
c
      ResoLimits(1,(mshell+1)) = r1max
      ResoLimits(2,(mshell+1)) = r2min_active
c
      write(pre,'(1x,a6,1x,f7.3,a3,f7.3)') 'Total:',r1max,' - ',
     +  r2min_active
      kshell = mshell + 1
      if (iwide_format.eq.1) then
        write(stdout,1001,advance='no') pre,
     +    NumUnique(1,kshell),
     +    NumRval(1,1,kshell),
     +    Rval(iGeneral,1,kshell),
     +    Rval(iAnomalous,1,kshell),
     +    Rval(iGeneral,2,kshell),
     +    Rval(iAnomalous,2,kshell),
     +    Rval(iGeneral,3,kshell),
     +    Rval(iAnomalous,3,kshell),
     +    NumIsigImrg(kshell),
     +    IsigImrg(kshell),
     +    min(1.0,Cmpltness(iGeneral,kshell)),
     +    min(1.0,Cmpltness(iAnomalous,kshell)),
     +    Multiplicity(iGeneral,kshell),
     +    Multiplicity(iAnomalous,kshell),
     +    CC_print(iGeneral,kshell),
     +    NumCc(iAnomalous,kshell),CC_print(iAnomalous,kshell),
     +    NumSigAno(kshell),SigAno(kshell)
        if (l_ell) then
          write(stdout,'(2f8.4)',advance='no')
     +      min(1.0,CmpltnessEll(iGeneral,kshell)),
     +      min(1.0,CmpltnessEll(iAnomalous,kshell))
        end if
        if(l_write_rej) then
          do j = 1, maxrej
            if (NumReject((mshell+1),j).gt.0) then
              write(stdout,'(i8)',advance='no') NumReject(kshell,j)
            end if
          end do
        end if
        write(stdout,'(/)',advance='no')

c       narrow format
      else
        write(stdout,2001,advance='no') pre,
     +    NumUnique(1,kshell),
     +    NumRval(1,1,kshell),
     +    Rval(iGeneral,1,kshell),
     +    Rval(iGeneral,2,kshell),
     +    Rval(iGeneral,3,kshell),
     +    NumIsigImrg(kshell),
     +    IsigImrg(kshell),
     +    min(1.0,Cmpltness(iGeneral,kshell)),
     +    min(1.0,Cmpltness(iAnomalous,kshell)),
     +    Multiplicity(iGeneral,kshell),
     +    Multiplicity(iAnomalous,kshell),
     +    CC_print(iGeneral,kshell),
     +    NumCc(iAnomalous,kshell),CC_print(iAnomalous,kshell),
     +    SigAno(kshell)
        if (l_ell) then
          write(stdout,'(2f8.4)',advance='no')
     +      min(1.0,CmpltnessEll(iGeneral,kshell)),
     +      min(1.0,CmpltnessEll(iAnomalous,kshell))
        end if
        write(stdout,'(/)',advance='no')
      end if
c
      write(stdout,'(/)')
c
      call write_table1(maxshell,mshell,iResoLimits,
     +  ResoLimits,Rval,NumIsigI,NumUnique,
     +  IsigImrg,Cmpltness,CmpltnessEll,Multiplicity,CC,SigAno,nrun)
c
 999  continue
c
c     resetting
c
      deallocate(jref2iref,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate jref2iref array')
      deallocate(jref2ibin,stat=istat)
      if (istat.ne.0)
     +  call error_exit('unable to deallocate jref2ibin array')
c
c----+=================================================================
c
      return
 101  format(1x,24a1,  8a1,  10a1,    8a1,    8a1,    8a1,    8a1,
     +                   8a1,    8a1,  8a1,    8a1,    8a1,    8a1,
     +                   7a1,    7a1,    8a1,
     +                   8a1,    8a1,    8a1,    8a1, 5a1)
 201  format(1x,24a1,  8a1,  10a1,    8a1,    8a1,    8a1,    8a1,
     +                   8a1,    8a1,  8a1,    8a1,    8a1,    8a1,
     +                   7a1,    7a1,    8a1,
     +                   8a1,    8a1,    8a1,    8a1)
 301  format(1x,24a1,  8a1,  10a1,    8a1,            8a1,        
     +                   8a1,          8a1,    8a1,    8a1,    8a1,
     +                   7a1,    7a1,    8a1,
     +                   8a1,    8a1,    8a1)
 102  format(    25x,   8x,   10x,7x,  a5,     4x,7x,  a5,     4x,
     +               7x,  a4,     5x,   8x,     8x,2x, a13,1x,    
     +               1x,         a13,     8x,
     +               8x,     8x,     8x,     8x,
     +               1x)
 202  format(    25x,   8x,   10x,7x,  a5,     4x,7x,  a5,     4x,
     +               7x,  a4,     5x,   8x,     8x,2x, a13,1x,    
     +                1x,         a13)
 302  format(    25x,   8x,   10x,     8x,     8x,     8x,     8x,
     +                    8x,2x, a13,1x,2x,    a12)
 103  format(    25x,   8x,   10x,1x,         a15,1x,         a15,
     +               1x,         a15,   8x,     8x,1x,         a15,
     +               1x,         a13,     8x,
     +               8x,     8x,     8x,     8x,
     +               1x)
 203  format(    25x,   8x,   10x,1x,         a15,1x,         a15,
     +               1x,         a15,   8x,     8x,1x,         a15,
     +               1x,         a13)
 303  format(    25x,   8x,   10x,     8x,             8x,        
     +                    8x,           8x,     8x,1x,         a15,
     +               1x,         a13)
 104  format( a22,3x,1x,a7,1x,a9,1x,  a7,1x,  a7,1x,  a7,1x,  a7,
     +               1x,  a7,1x,  a7,1x,a7,1x,  a7,1x,  a7,1x,  a7,
     +               1x,  a6,1x,  a6,1x,  a7,
     +               1x,a7,1x, a7 ,1x,a7,1x, a7 ,2x)
 304  format( a22,3x,1x,a7,1x,a9,1x,  a7,        1x,  a7,        
     +               1x,  a7,        1x,a7,1x,  a7,1x,  a7,1x,  a7,
     +               1x,  a6,1x,  a6,1x,  a7,
     +               1x,a7,1x, a7 ,1x, a7 ,2x)
 1001 format(    a25,1x,i7,1x,i9,1x,f7.3,1x,f7.3,1x,f7.3,1x,f7.3,
     +               1x,f7.3,1x,f7.3,1x,i7,1x,f7.3,1x,f7.4,1x,f7.4,
     +               1x,f6.2,1x,f6.2,1x,f7.4,
     +               1x,i7,1x,f7.4,1x,i7,1x,f7.3,1x)
 2001 format(    a25,1x,i7,1x,i9,1x,f7.3,        1x,f7.3,        
     +               1x,f7.3,        1x,i7,1x,f7.3,1x,f7.4,1x,f7.4,
     +               1x,f6.2,1x,f6.2,1x,f7.4,
     +               1x,i7,1x,f7.4,1x,f7.3,1x)
      end
