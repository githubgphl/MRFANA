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

      subroutine write_info()
c
      implicit none
c
      include 'units.inc'
c
c----+=================================================================
c

      write(stdout,'(/,2a)')
     +  ' ####--------------------------------------',
     +  '--------------------------'
      write(stdout,'(1a)') ' #### NOTES'
      write(stdout,'(2a,/)')
     +  ' ####--------------------------------------',
     +  '--------------------------'


      write(stdout,'(2a)')
     +  ' ####--------------------------------------',
     +  '--------------------------'
      write(stdout,'(1a)') ' #### Rmerge:'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  According to [1] and Evans (2006) this ',
     +  'is defined as'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '                   sum_i ( sum_j | Ij - ',
     +  '<I> | )'
      write(stdout,'(2a)') '    Rmerge(I) = ------------------------',
     +  '--------                 (1)'
      write(stdout,'(2a)') '                   sum_i ( sum_j        ',
     +  '<I>   )'
      write(stdout,'(1a)') ' '
      write(stdout,'(1a)') '  with'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '    Ij  = the intensity of the jth ',
     +  'observation of reflection i'
      write(stdout,'(2a)') '    <I> = the mean of the intensities of ',
     +  'all observations of reflection i'
      write(stdout,'(1a)') ' '
      write(stdout,'(1a)') '    sum_i is taken over all reflections'
      write(stdout,'(2a)') '    sum_j is taken over all observations ',
     +  'of each reflection.'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  Diederichs & Karplus (1997), Weiss & ',
     +  'Hilgenfeld (1997) and'
      write(stdout,'(2a)') '  Weiss (2001) define this slightly ',
     +  'different as'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '                   sum_i ( sum_j | Ij - ',
     +  '<I> | )'
      write(stdout,'(2a)') '    Rmerge(I) = ------------------------',
     +  '--------                 (2)'
      write(stdout,'(2a)') '                   sum_i ( sum_j   Ij   ',
     +  '      )'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  which is identically to (1) if <I> ',
     +  'is the mean (or average)'
      write(stdout,'(2a)') '  intensity. However, it is the weighted ',
     +  '(by variance) mean of'
      write(stdout,'(2a)') '  the intensities we are computing ',
     +  'during merging:'
      write(stdout,'(1a)') ' '
      write(stdout,'(1a)') '           sum_j ( Wj * Ij )'
      write(stdout,'(2a)') '    <I> = -------------------            ',
     +  '                        (3)'
      write(stdout,'(1a)') '           sum_j ( Wj      )'
      write(stdout,'(1a)') ' '
      write(stdout,'(1a)') '                    1     '
      write(stdout,'(2a)') '    with:  Wj = ---------                ',
     +  '                        (4)'
      write(stdout,'(1a)') '                sd(Ij)**2'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  Then equations (1) and (2) become ',
     +  'different.'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  Drenth (1994) calls (1) Rmerge and ',
     +  '(2) Rsym.'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)')
     +  ' ####--------------------------------------',
     +  '--------------------------'
      write(stdout,'(1a)') ' #### Rmeas (also called Rrim):'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  The redundancy-independent merging R ',
     +  'factor value Rrim, also'
      write(stdout,'(2a)') '  denoted Rmeas, for merging all ',
     +  'intensities in a given shell.'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '           sum_i [Ni /( Ni - 1)]^(1/2) ',
     +  'sum_j | Ij - <Ii> |'
      write(stdout,'(2a)') '    Rrim = -------------------------',
     +  '----------------------       (5)'
      write(stdout,'(2a)') '           sum_i (                  ',
     +  '   sum_j   Ij        )'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '    Ij   = the intensity of the jth ',
     +  'observation of reflection i'
      write(stdout,'(2a)') '    <Ii> = the mean of the intensities ',
     +  'of all observations of reflection i'
      write(stdout,'(2a)') '    Ni   = the redundancy (the number ',
     +  'of times reflection i has been measured).'
      write(stdout,'(1a)') ' '
      write(stdout,'(1a)') '    sum_i is taken over all reflections'
      write(stdout,'(2a)') '    sum_j is taken over all observations ',
     +  'of each reflection.'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  See Diederichs & Karplus (1997), ',
     +  'Weiss & Hilgenfeld (1997), Weiss'
      write(stdout,'(1a)') '  (2001) and [2].'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  Evans (2006) defines this differently, ',
     +  'with the denominator'
      write(stdout,'(2a)') '  summing over <I> similar to the Rmerge ',
     +  'equation (1) above:'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '           sum_i [Ni /( Ni - 1)]^(1/2) ',
     +  'sum_j | Ij - <Ii> |'
      write(stdout,'(2a)') '    Rrim = ----------------------------',
     +  '-------------------       (6)'
      write(stdout,'(2a)') '           sum_i (                    ',
     +  ' sum_j        <Ii> )'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)')
     +  ' ####--------------------------------------',
     +  '--------------------------'
      write(stdout,'(1a)') ' #### Rpim:'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  The precision-indicating merging R ',
     +  'factor value Rpim,'
      write(stdout,'(2a)') '  for merging all intensities in this ',
     +  'data set is defined in [3],'
      write(stdout,'(2a)') '  Diederichs & Karplus (1997), Weiss & ',
     +  'Hilgenfeld (1997) and'
      write(stdout,'(1a)') '  Weiss (2001) as'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '           sum_i [  1/(Ni - 1)]^(1/2) ',
     +  'sum_j | Ij - <Ii> |'
      write(stdout,'(2a)') '    Rpim = -------------------------',
     +  '---------------------        (7)'
      write(stdout,'(2a)') '           sum_i (                  ',
     +  '  sum_j   Ij        )'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '  As above, Evans (2006) defines this ',
     +  'differently as:'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') '           sum_i [  1/(Ni - 1)]^(1/2) ',
     +  'sum_j | Ij - <Ii> |'
      write(stdout,'(2a)') '    Rpim = -------------------------',
     +  '---------------------        (8)'
      write(stdout,'(2a)') '           sum_i (                  ',
     +  '  sum_j        <Ii> )'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') ' To add more complexity, all of the ',
     +  'above statistics (Rmerge, Rmeas'
      write(stdout,'(2a)') ' and Rpim) can be computed over all ',
     +  'measurements of a unique reflection'
      write(stdout,'(2a)') ' or separately for I+ and I-. The latter ',
     +  'has the advantage that an'
      write(stdout,'(2a)') ' existing anomalous signal will not ',
     +  'result in an increase of R-value'
      write(stdout,'(2a)') ' because of this. However, the original ',
     +  'definitions above are not'
      write(stdout,'(2a)') ' taking anomalous signal into account ',
     +  'this way.'
      write(stdout,'(1a)') ' '
      write(stdout,'(1a)') ' ### References:'
      write(stdout,'(1a)') ' '
      write(stdout,'(2a)') ' [1] http://mmcif.wwpdb.org/dictionaries',
     +  '/mmcif_std.dic/Items/_reflns_shell.Rmerge_I_obs.html'
      write(stdout,'(2a)') ' [2] http://www.iucr.org/__data/iucr',
     +  '/cifdic_html/2/mmcif_pdbx.dic/Ireflns.pdbx_Rrim_I_all.html'
      write(stdout,'(2a)') ' [3] http://www.iucr.org/__data/iucr',
     +  '/cifdic_html/2/mmcif_pdbx.dic/Ireflns.pdbx_Rpim_I_all.html'
      write(stdout,'(2a)') ' [4] Diederichs, K. & Karplus, P. A. ',
     +  '(1997). Nature Struct. Biol. 4, 269-275.'
      write(stdout,'(3a)') ' [5] Drenth, J. (1994). Principles ',
     +  'of Protein X-ray Crystallography. Springer-Verlag, ',
     +  'New York. p289-290.'
      write(stdout,'(1a)') ' [6] Evans, P. (2006). Acta D62, 72-82.'
      write(stdout,'(2a)') ' [7] Weiss, M. S. & Hilgenfeld, R. ',
     +  '(1997). J. Appl. Cryst. 30, 203-205.'
      write(stdout,'(2a)') ' [8] Weiss, M. S. (2001). J. Appl. ',
     +  'Cryst. 34, 130-135.'
      write(stdout,'(1a)') ' '
c
c----+=================================================================
c
      return
      end
