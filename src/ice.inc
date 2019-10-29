c     **********************************************************************
c
c      Copyright (c) 2010, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     handling of ice-rings

      integer    NumIceRings
      parameter (NumIceRings = 21)
      integer    MaxNumIceRings
      parameter (MaxNumIceRings = 100)
c
c 20100705:
c
c   experimentally determined ice-ring ranges from deposited
c   diffraction data (mainly JCSG):
c
c     (1) looking at all non-indexed spots when processing with autoPROC
c
c     (2) fitting Gaussian via fityk
c
c     +- 3sigma            Center  Height   Area    FWHM    parameters...
c   
c    1.713 - 1.721 A  :   1.71697  117.84  0.364 0.00290464    117.84 1.71697 0.00145232
c    1.872 - 1.888 A  :   1.87994  442.04  2.393 0.00508574    442.04 1.87994 0.00254287
c    1.906 - 1.922 A  :   1.91419 1493.04  8.486 0.00533977   1493.04 1.91419 0.00266989
c    1.936 - 1.955 A  :   1.94543  370.77  2.575 0.00652518    370.77 1.94543 0.00326259
c    2.056 - 2.076 A  :   2.06608 1970.19 13.623 0.00649583   1970.19 2.06608 0.00324792
c    2.235 - 2.257 A  :   2.24585 2867.87 22.691 0.00743305   2867.87 2.24585 0.00371653
c    2.274 - 2.299 A  :   2.28658  185.37  1.618 0.00819938    185.37 2.28658 0.00409969
c    2.495 - 2.516 A  :   2.50555  324.42  2.506 0.00725785    324.42 2.50555 0.00362892
c    2.649 - 2.681 A  :   2.66487 1210.62 13.955 0.01082920   1210.62 2.66487 0.00541458
c    2.786 - 2.816 A  :   2.80123  155.91  1.663 0.01002100    155.91 2.80123 0.00501048
c    3.219 - 3.241 A  :   3.23026  102.01  0.809 0.00744714    102.01 3.23026 0.00372357
c    3.420 - 3.442 A  :   3.43134 2298.85 17.764 0.00725936   2298.85 3.43134 0.00362968
c    3.644 - 3.664 A  :   3.65415 2734.10 19.685 0.00676383   2734.10 3.65415 0.00338191
c    3.873 - 3.900 A  :   3.88620 1808.54 17.102 0.00888359   1808.54 3.88620 0.00444179
c
c 20161125: added higher resolution rings:
c
c    Center     FWHM        Sigma         +- 3sigma
c | 1.52130 | 0.00206048 | 0.00087494 | 1.524 - 1.519 |
c | 1.47139 | 0.00123968 | 0.00052640 | 1.473 - 1.470 |
c | 1.44217 | 0.00175535 | 0.00074537 | 1.444 - 1.440 |
c | 1.37013 | 0.00129781 | 0.00055109 | 1.372 - 1.368 |
c | 1.36498 | 0.00160928 | 0.00068335 | 1.367 - 1.363 |
c | 1.29755 | 0.00144315 | 0.00061280 | 1.299 - 1.296 |
c | 1.27455 | 0.00072355 | 0.00030724 | 1.275 - 1.274 |
c | 1.26013 | 0.00098256 | 0.00041722 | 1.261 - 1.259 |
c | 1.22289 | 0.00081936 | 0.00034792 | 1.224 - 1.222 |
c | 1.16960 | 0.00092685 | 0.00039357 | 1.171 - 1.168 |
c | 1.12296 | 0.00101736 | 0.00043200 | 1.124 - 1.122 |
c
      double precision IceRings(2,NumIceRings)
      data IceRings /
     +  3.906d0,3.873d0,
     +  3.671d0,3.642d0,
     +  3.450d0,3.420d0,
     +  2.681d0,2.649d0,
     +  2.258d0,2.235d0,
     +  2.076d0,2.056d0,
     +  1.955d0,1.936d0,
     +  1.922d0,1.906d0,
     +  1.888d0,1.872d0,
     +  1.721d0,1.713d0,
     +  1.524d0,1.519d0,
     +  1.473d0,1.470d0,
     +  1.444d0,1.440d0,
     +  1.372d0,1.368d0,
     +  1.367d0,1.363d0,
     +  1.299d0,1.296d0,
     +  1.275d0,1.274d0,
     +  1.261d0,1.259d0,
     +  1.224d0,1.222d0,
     +  1.171d0,1.168d0,
     +  1.124d0,1.122d0/

      double precision IceRingsUsed(2,MaxNumIceRings)
      common /d_ice/ IceRingsUsed

      logical consider_ice
      common /l_ice/ consider_ice

      integer IceRingsNref(MaxNumIceRings)
      common /i_ice/ IceRingsNref
