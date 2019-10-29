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

c     maximum number of columns (e.g. in MTZ files)
      integer    maxcol
      parameter (maxcol=500)

c     maximum number of datasets
      integer    maxdatasets
      parameter (maxdatasets=100)

c     maximum number of shells/bins we support (remember that we need 1
c     for the overall values)
      integer    maxshell
      parameter (maxshell=201)

c     maximum number of runs
      integer    maxruns
      parameter (maxruns=1000)

c     how many EXCLUDE cards can we have on a RUN definition:
      integer    maxrunexcl
      parameter (maxrunexcl=20)

c     maximum multiplicity we will support
      integer    maxmul
      parameter (maxmul=9999)
c
      integer    minbin,   maxbin
      parameter (minbin=1, maxbin=100)
c
      double precision reseps
      parameter (reseps=0.0001d0)
