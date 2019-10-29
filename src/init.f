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

      subroutine init(nref,ncol,ranges,clabs,ilabs)
c
      implicit none
c
      include 'units.inc'
      include 'params.inc'
      include 'info.inc'
      include 'runs.inc'
      include 'cmpl.inc'
      include 'xml.inc'
c
      character*30 clabs(maxcol)
      integer nref(2), ncol, ilabs(maxcol), today(3), now(3), values(8)
      double precision ranges(2, maxcol)
c
c----+=================================================================
c
      call date_and_time(VALUES=values)
      today(1) = values(3)
      today(2) = values(2)
      today(3) = values(1)
      now(1:3) = values(5:7)
      if(today(3).lt.100) today(3)=today(3)+2000
      write(xml_timestamp,
     +  '(i4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2)')
     +  today(3),'-',today(2),'-',today(1),
     +  now(1),':',now(2),':',now(3)
c
      call ccp4h_init_lib(-1,-1)
c
      write(stdout,'(/)')
      write(stdout,'(3a)') ' #####################################'
     +  ,'#####################################'
     +  ,'#######################'
      write(stdout,'(a)')
     .  ' ## MRFANA - metrics for unmerged X-Ray reflection data'
      write(stdout,'(3a,/)') ' #####################################'
     +  ,'#####################################'
     +  ,'#######################'

      write(stdout,'(a,/)')
     +  ' Copyright (c) 2007, 2019 Global Phasing Ltd.'
      write(stdout,'(a)')
     +  ' Authors: Clemens Vonrhein, Claus Flensburg,'//
     +  ' Ian J. Tickle and'
      write(stdout,'(a,/)')
     +  '          Gerard Bricogne'

      write(stdout,'(3a)') ' #####################################'
     +  ,'#####################################'
     +  ,'#######################'

      write(stdout,'(/,a,i4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2,/)')
     .  ' Current date and time : ',today(3),'-',today(2),'-',today(1),
     +  now(1),':',now(2),':',now(3)
c
      lun = 10
      nref(1) = 0
      nref(2) = 0
      ncol = 0
      ranges = 0.0
      clabs = ' '
      ilabs = 0
c
      run_reslim(1,0) = -999999.
      run_reslim(2,0) =  999999.
      run_nref = 0
c
      cmpl_nshell = 0
c
c----+=================================================================
c
      return
      end
