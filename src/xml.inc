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

c     XML related

      character xml_sg*10, xml_timestamp*19,
     +          xml_base*1024, xml_file*1024
      common /xml_c/ xml_sg, xml_timestamp, xml_base, xml_file
c
      double precision xml_cell(6), xml_wave
      common /xml_d/ xml_cell, xml_wave
