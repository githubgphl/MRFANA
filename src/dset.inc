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

c     dataset model following the MTZ data model
c
c       http://www.ccp4.ac.uk/html/mtzformat.html#datamodel

      character pname(maxdatasets), xname(maxdatasets),
     +  dname(maxdatasets)
      common /dset_c/ pname, xname, dname

      integer isets(maxdatasets), ndatasets
      common /dset_i/isets, ndatasets

      double precision datcell(6,maxdatasets)
      common /dset_d/ datcell
