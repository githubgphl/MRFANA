c     **********************************************************************
c
c      Copyright (c) 2015, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     multiple tests to see if that reflection is active (not rejected,
c     inside limits etc)

      subroutine refl_test_active(ncoltot,nref,refl,nrun,irun,
     +  jref2iref,njref,L_Bin,N_Bin,Bin,nrefbin_use,jref2ibin,
     +  n_jref2ibin,jref2reso)
c
      implicit none
c
      include 'params.inc'
      include 'lookup.inc'
      include 'reso.inc'
      include 'aniso.inc'
      include 'args.inc'
c
      integer ncoltot,nref(2),nrun
      double precision refl(ncoltot,nref(1))
      integer irun(nrun)
c
      integer jref2iref((nref(2)+1))
      integer jref2ibin((nref(2)+1))
      integer njref, N_Bin, nrefbin_use, n_jref2ibin
      logical L_Bin
c
      integer iref, jrun, nrej(7), i, istat, hklold(3), hkl(3), j, k,
     +  jref
      character*30 crej(7)
      data crej /'not in current run','not in any run'
     +  ,'outside resolution','partial reflection'
     +  ,'outside partiality-check','FLAGed','outside ellipsoid'/
      double precision reso, jref2reso(1:2,1:maxshell), Bin(N_bin),
     .  jf
      logical l_run
c
      real, allocatable :: resoref(:)
      integer, allocatable :: reso2ref(:,:)
      integer, allocatable :: isort(:)
      integer nreso, ibin, NbeyondNBin
c
c----+=================================================================
c
      n_jref2ibin = 0
      NbeyondNBin = 0

      if (L_bin.and.(N_Bin.eq.0)) then
        allocate(resoref(1:(nref(2)+1)),stat=istat)
        if (istat.ne.0)
     +    call error_exit('unable to allocate resoref array')
        allocate(reso2ref(1:2,1:(nref(2)+1)),stat=istat)
        if (istat.ne.0)
     +    call error_exit('unable to allocate reso2ref array')
        allocate(isort(1:(nref(2)+1)),stat=istat)
        if (istat.ne.0)
     +    call error_exit('unable to allocate isort array')
      end if

      jref2iref = 0
      njref = 0
      nreso = 0

      do i = 1, 6
        nrej(i) = 0
      end do

      current_reslow  = 0.0d0
      current_reshigh = 9999.0d0
      hklold = 0

      do iref = 1, nref(2)

        if (irun(1).gt.0) then

          l_run = .false.
          do jrun = 1, nrun
            if (nint(refl(icol(i_run),iref)).eq.irun(jrun)) then
              l_run = .true.
              exit
            end if
          end do

c         not in current run
          if (.not.l_run) then
            nrej(1) = nrej(1) + 1
            cycle
          end if

        end if

c       not in any run:
        if (nint(refl(icol(i_run),iref)).eq.0) then
          nrej(2) = nrej(2) + 1
          cycle
        end if

c       outside overall resolution limits:
        reso = 1.d0/sqrt(refl(icol(i_resol),iref))
        if ((reso.gt.(reslow +reseps)).or.
     +      (reso.lt.(reshigh-reseps))    ) then
          nrej(3) = nrej(3) + 1
          cycle
        end if

c       outside ellipsoidal resolution limit:
        if (l_ell) then
          if ((reso.lt.(reshigh_ell-reseps))    ) then
            nrej(7) = nrej(7) + 1
            cycle
          end if
        end if

c       skip partial ones
        if (nint(refl(icol(i_misym),iref)).gt.256) then
          nrej(4) = nrej(4) + 1
          cycle
        end if

c       outside partiality-check
        if (icol(i_frac).gt.0) then
          if ( (refl(icol(i_frac),iref).gt.0.0d0)              .and.
     +        ((refl(icol(i_frac),iref).lt.check_part_min).or.
     +         (refl(icol(i_frac),iref).gt.check_part_max)    )) then
            nrej(5) = nrej(5) + 1
            cycle
          end if
        end if

c  http://www.ccp4.ac.uk/html/scala.html#files
c        FLAG            error flag (packed bits) from Mosflm (v6.2.3
c                        or later). By default, if this column is present, 
c                        observations with a non-zero FLAG will be
c                        omitted. They may be conditionally accepted
c                        using the ACCEPT command (qv)
c                        Bit flags:
c                              1         BGRATIO too large
c                              2         PKRATIO too large
c                              4         Negative > 5*sigma
c                              8         BG Gradient too high
c                             16         Profile fitted overload
c                             32         Profile fitted "edge" reflection
c       flagged as bad
        if (icol(i_flag).gt.0) then
          if (nint(refl(icol(i_flag),iref)).ne.0) then
            nrej(6) = nrej(6) + 1
            cycle
          end if
        end if

        njref = njref + 1
        jref2iref(njref) = iref

        current_reslow  = max(reso,current_reslow)
        current_reshigh = min(reso,current_reshigh)
c
c       setup binning by equal numbers
c
        if (L_bin) then
          hkl(1) = nint(refl(icol(i_h),iref))
          hkl(2) = nint(refl(icol(i_k),iref))
          hkl(3) = nint(refl(icol(i_l),iref))
          if (N_bin.eq.0) then
            if ((hkl(1).ne.hklold(1)).or.
     +          (hkl(2).ne.hklold(2)).or.
     +          (hkl(3).ne.hklold(3))    ) then
              nreso = nreso + 1
              resoref(nreso)  = real(reso,4)
              reso2ref(1,nreso) = njref
              isort(nreso) = nreso
            end if
            reso2ref(2,nreso) = njref
          else
            if ((hkl(1).ne.hklold(1)).or.
     +          (hkl(2).ne.hklold(2)).or.
     +          (hkl(3).ne.hklold(3))    ) then
              ibin = 0
              do i = 1, N_Bin
                if (reso.ge.Bin(i)) then
                  ibin = i
                  exit
                end if
              end do
              if (ibin.eq.0) then
                ibin = N_Bin + 1
                NbeyondNBin = NbeyondNBin + 1
              end if
            end if
            jref2ibin(njref) = ibin
          end if
          hklold = hkl
        end if

      end do

      if (njref.eq.0) then
        do i = 1, 7
          if (nrej(i).eq.nref(2)) then
            call error_exit
     +        ('no active reflections available - '
     +        //'all rejected because '''//trim(crej(i))//'''')
          end if
        end do
      end if

      if (L_Bin) then

        if ((N_Bin.eq.0).and.(nreso.gt.0)) then

          call rsortd(resoref,nreso,isort)
c
          if (iset_nshell.eq.0) then
c           how many bins would we get:
            j = nint(real(nreso)/nrefbin_use)
c           does that fit:
            if (j.gt.(maxshell-1)) then
c             no - so restrict to maxshell
              jf = real(nreso)/(maxshell-1)
            else
c             yes - so use it
              jf = real(nreso)/j
            end if
          else
c           have -n and -nref -<N> argument: so want to do equal-number
c           binning with a maximum number of bins of at least a certain
c           number of reflections
            jf = real(nreso)/nshell
            if (jf.lt.nrefbin_use) then
              j = nint(real(nreso)/nrefbin_use)
              jf = real(nreso)/j
            end if

          end if

          k = 0
          n_jref2ibin = 1
          jref2reso(1,n_jref2ibin) = resoref(1)
          do i = 1, nreso
            k = k + 1
            if (k.gt.nint(jf*n_jref2ibin)) then
              n_jref2ibin = n_jref2ibin + 1
              jref2reso(1,n_jref2ibin) = resoref(i)
            end if

c           mark all 
            do jref = reso2ref(1,isort(i)), reso2ref(2,isort(i))
              jref2ibin(jref) = n_jref2ibin
            end do

            jref2reso(2,n_jref2ibin) = resoref(i)

          end do
          jref2reso(1,(n_jref2ibin+1)) = jref2reso(1,1)
          jref2reso(2,(n_jref2ibin+1)) = jref2reso(2,n_jref2ibin)
          deallocate(resoref)
          deallocate(reso2ref)
          deallocate(isort)

c
        else if (N_Bin.gt.0) then

          jref2reso(1,1) = current_reslow
          jref2reso(2,1) = Bin(1)
          do i = 2, N_Bin
            jref2reso(1,i) = Bin((i-1)) - 0.000001
            jref2reso(2,i) = Bin(i)
          end do
          n_jref2ibin = N_Bin
          if (jref2reso(2,n_jref2ibin).gt.current_reshigh) then
            n_jref2ibin = n_jref2ibin + 1
            jref2reso(1,n_jref2ibin) = Bin((n_jref2ibin-1)) - 0.000001
            jref2reso(2,n_jref2ibin) = current_reshigh
          else
            if (NbeyondNBin.gt.0) then
              call error_exit('NbeyondNBin>0')
            end if
          end if

          jref2reso(1,(n_jref2ibin+1)) = jref2reso(1,1)
          jref2reso(2,(n_jref2ibin+1)) = jref2reso(2,n_jref2ibin)

        end if
      end if
c
c----+=================================================================
c
      return
      end
