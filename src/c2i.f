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
c
c     Multi-record files will come in different formats with different
c     meanings of colums. Here we establish a common system of pointers
c     to simplify working with the reflection array later.
c
c INTEGRATE.HKL (See: http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#INTEGRATE.HKL)
c
c  H     : h-index of the reflection.
c
c  K     : k-index of the reflection.
c
c  L     : l-index of the reflection.
c
c  IOBS  : intensity of the reflection obtained by profile fitting. The
c          intensity is already LP-corrected assuming an unpolarized
c          incident beam. The final polarization correction is carried
c          out in the CORRECT step.
c
c  SIGMA : e.s.d. of IOBS as obtained from profile fitting.
c
c  XCAL  : calculated detector X-coordinate of the reflection.
c
c  YCAL  : calculated detector Y-coordinate of the reflection.
c
c  ZCAL  : calculated image number of reflection at diffraction maximum.
c
c  RLP   : Reciprocal Lorentz-Polarization factor computed for unpolarized
c          incident beam.
c
c  PEAK  : Percentage of observed reflection intensity.
c
c  CORR  : Correlation factor between observed and expected reflection
c          profile.
c
c  MAXC  : Largest image pixel value contributing to the reflection.
c
c  XOBS  : observed detector X-coordinate of the reflection. 0 if
c          unobserved.
c
c  YOBS  : observed detector Y-coordinate of the reflection. 0 if
c          unobserved.
c
c  ZOBS  : observed centroid of image numbers of reflection. 0 if
c          unobserved.
c
c  ALF0  : incident beam direction described by polar coordinates
c          (degrees) as sin(bet0)cos(alf0), sin(bet0)sin(alf0),
c          cos(bet0).
c
c  BET0  : incident beam direction described by polar coordinate as
c          sin(bet0)cos(alf0), sin(bet0)sin(alf0), cos(bet0).
c
c  ALF1  : diffracted beam direction described by polar coordinates as
c          sin(bet1)cos(alf1), sin(bet1)sin(alf1), cos(bet1).
c
c  BET1  : diffracted beam direction described by polar coordinates as
c          sin(bet1)cos(alf1), sin(bet1)sin(alf1), cos(bet1).
c
c  PSI   : angle (degrees) as defined by Schwarzenbach and Flack that
c          completely specifies diffraction geometry with respect to the
c          crystal cell axes.
c
c XDS_ASCII.HKL (See: http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#XDS_ASCII.HKL)
c
c  H           : reflection index h.
c
c  K           : reflection index k.
c
c  L           : reflection index l.
c
c  IOBS        : reflection intensity.
c
c  SIGMA(IOBS) : error of the intensity. A negative sign is attached to
c                SIGMA(IOBS) to indicate a MISFIT (IOBS is incompatible
c                with symmetry equivalent reflection intensities). At
c                present, such MISFITs are ignored in the subsequent
c                processing by xscale.
c
c  XD          : X-coordinate (pixels) of reflection on detector.
c
c  YD          : Y-coordinate (pixels) of reflection on detector.
c
c  ZD          : centroid of image numbers that recorded the Bragg peak.
c
c  RLP         : reciprocal LP-correction factor that has been applied to this
c                reflection.
c
c  PEAK        : percentage of observed reflection intensity. A value less than
c                100.0 indicates that the reflection was either incompletely recorded
c                or overlapping with neighbouring reflections or included untrusted
c                pixels in the profile region. Note, that the reflection intensity
c                given in ITEM_IOBS always refers to the complete reflection
c                regardless of the value of ITEM_PEAK.
c
c  CORR        : percentage of correlation between observed and expected
c                reflection profile.
c
c  PSI         : value of Schwarzenbach and Flack's psi-angle (degrees) which
c                completely specifies diffraction geometry with respect to the unit
c                cell axes.
c
c XSCALE (See: http://xds.mpimf-heidelberg.mpg.de/html_doc/xscale_parameters.html#OUTPUT_FILE=)
c
c  H           : reflection index h.
c
c  K           : reflection index k.
c
c  L           : reflection index l.
c
c  IOBS        : reflection intensity.
c
c  SIGMA(IOBS) : error of the intensity
c
c  MERGE=FALSE:
c
c    XD        : X-coordinate (pixels) of reflection on detector.
c
c    YD        : Y-coordinate (pixels) of reflection on detector.
c
c    ZD        : centroid of image numbers that recorded the Bragg peak.
c
c    PSI       : value of Schwarzenbach and Flack's psi-angle (degrees) which
c                completely specifies diffraction geometry with respect to the unit
c                cell axes.
c
c    ISET      : identifying number of the input data set the reflection came from.
c
c
c    DCY       : decay-factor b(h) for this reflection used for 0-dose correction
c                (this item will only be present if a value has been specified for
c                CRYSTAL_NAME=)
c
c MTZ: (See: http://www.ccp4.ac.uk/html/mtzformat.html)
c
c Compulsory columns:
c
c         H K L           indices
c         M/ISYM          partial flag, symmetry number
c         BATCH           batch number
c         I               intensity  (integrated intensity)
c         SIGI            sd(intensity)   (integrated intensity)
c
c Optional columns:
c
c         XDET YDET       position on detector of this reflection: these
c                         may be in any units (e.g. mm or pixels), but the
c                         range of values must be specified in the
c                         orientation data block for each batch. If
c                         these columns are absent, the scale may not be
c                         varied across the detector (i.e. only SCALES
c                         DETECTOR 1 is valid)
c         ROT             rotation angle of this reflection ("Phi"). If
c                         this column is absent, only SCALES BATCH is valid.
c         IPR             intensity  (profile-fitted intensity)     
c         SIGIPR          sd(intensity)   (profile-fitted intensity)
c         SCALE           previously calculated scale factor (e.g. from
c                         previous run of Scala). This will be applied
c                         on input
c         SIGSCALE        sd(SCALE)
c         TIME            time for B-factor variation (if this is
c                         missing, ROT is used instead)
c         MPART           partial flag from Mosflm
c         FRACTIONCALC    calculated fraction, required to SCALE PARTIALS
c         LP              Lorentz/polarization correction (already applied)
c         FLAG            error flag (packed bits) from Mosflm (v6.2.3
c                         or later). By default, if this column is present, 
c                         observations with a non-zero FLAG will be
c                         omitted. They may be conditionally accepted
c                         using the ACCEPT command (qv)
c                         Bit flags:
c                               1         BGRATIO too large
c                               2         PKRATIO too large
c                               4         Negative > 5*sigma
c                               8         BG Gradient too high
c                              16         Profile fitted overload
c                              32         Profile fitted "edge" reflection
c         BGPKRATIOS      packed background & peak ratios, & background
c                         gradient, from Mosflm, to go with FLAG
c
c #----------------------------------------------------------------------
c
c   Column         MTZ            XDS_ASCII.HKL   INTEGRATE.HKL
c   H              H              H               H
c   K              K              K               K
c   L              L              L               L
c   M/ISYM         M/ISYM
c   BATCH          BATCH          
c   ZOBS                          ZD              ZOBS
c   ZCAL                          ZD              ZCAL
c   I              I              IOBS            IOBS
c   SIGI           SIGI           SIGMA(IOBS)     SIGMA
c   XOBS           XDET           XD              XOBS
c   YOBS           YDET           YD              YOBS
c   XCAL           XDET           XD              XCAL
c   YCAL           YDET           YD              YCAL
c   ROT            ROT
c   IPR            IPR
c   SIGIPR         SIGIPR
c   SCALE          SCALE
c   SIGSCALE       SIGSCALE
c   TIME           TIME
c   MPART          MPART
c   FRACTIONCALC   FRACTIONCALC
c   LP             LP             RLP             RLP
c   FLAG           FLAG
c   BGPKRATIOS     BGPKRATIOS
c
      integer function c2i(filtyp,c,ncol,clabs,ilabs)
c
      implicit none
c
      include 'params.inc'
c
      character*(*) filtyp, c
      character*30 cc,clabs(maxcol)
      integer ncol, ilabs(maxcol)
c
      integer i
c
c----+=================================================================
c
      cc = ' '
      if (filtyp(1:4).eq.'XDS_') then

c       we support XDS_ASCII and XSCALE

        if (c.eq.'I') then
          cc = 'IOBS'
        else if (c.eq.'SIGI') then
          if ((filtyp(1: 9).eq.'XDS_ASCII' ).or.
     +        (filtyp(1:10).eq.'XDS_XSCALE')    ) then
            cc = 'SIGMA(IOBS)'
          else
            cc = 'SIGMA'
          end if
        else if (c.eq.'XOBS') then
          if ((filtyp(1: 9).eq.'XDS_ASCII' ).or.
     +        (filtyp(1:10).eq.'XDS_XSCALE')    ) then
            cc = 'XD'
          else
            cc = 'XOBS'
          end if
        else if (c.eq.'XCAL') then
          if ((filtyp(1: 9).eq.'XDS_ASCII' ).or.
     +        (filtyp(1:10).eq.'XDS_XSCALE')    ) then
            cc = 'XD'
          else
            cc = 'XCAL'
          end if
        else if (c.eq.'YOBS') then
          if ((filtyp(1: 9).eq.'XDS_ASCII' ).or.
     +        (filtyp(1:10).eq.'XDS_XSCALE')    ) then
            cc = 'YD'
          else
            cc = 'YOBS'
          end if
        else if (c.eq.'YCAL') then
          if ((filtyp(1: 9).eq.'XDS_ASCII' ).or.
     +        (filtyp(1:10).eq.'XDS_XSCALE')    ) then
            cc = 'YD'
          else
            cc = 'YCAL'
          end if
        else if (c.eq.'ZOBS') then
          if ((filtyp(1: 9).eq.'XDS_ASCII' ).or.
     +        (filtyp(1:10).eq.'XDS_XSCALE')    ) then
            cc = 'ZD'
          else
            cc = 'ZOBS'
          end if
        else if (c.eq.'ZCAL') then
          if ((filtyp(1: 9).eq.'XDS_ASCII' ).or.
     +        (filtyp(1:10).eq.'XDS_XSCALE')    ) then
            cc = 'ZD'
          else
            cc = 'ZCAL'
          end if
        else if (c.eq.'LP') then
          cc = 'RLP'
        else if (c.eq.'FRACTIONCALC') then
          cc = 'PEAK'
        end if
      end if
c
      if (cc.eq.' ') cc = c
c
      c2i = 0
      do i = 1, ncol
        if (clabs(i).eq.cc) then
          c2i = ilabs(i)
          return
        end if
      end do
c
c----+=================================================================
c
      return
      end
