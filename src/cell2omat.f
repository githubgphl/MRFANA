c     **********************************************************************
c
c      Copyright (c) 2017, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein and Gerard Bricogne
c
c     **********************************************************************

c     define orthogonalisation matrix

      subroutine cell2omat(cell,omat)
c
      implicit none
c
      double precision cell(6)
      double precision omat(3,3)
c
      double precision d2r,
     +  COSA,
     +  COSB,
     +  COSG,
     +  SINB,
     +  SING,
     +  COSAS,
     +  SINAS
      parameter (d2r = 0.01745329251994329576d0 )
c
c----+=================================================================
c
      COSA=cos(cell(4)*d2r)
      COSB=cos(cell(5)*d2r)
      COSG=cos(cell(6)*d2r)
      SINB=sin(cell(5)*d2r)
      SING=sin(cell(6)*d2r)
      COSAS = ((COSG*COSB)-COSA)/ (SINB*SING)
      SINAS = 1.0d0-(COSAS*COSAS)
      if(SINAS.lt.0.0d0)SINAS=0.0d0
      SINAS=sqrt(SINAS)

      omat = 0.0d0
c     XO along a  Zo along c*
      omat(1,1) =  cell(1)
      omat(1,2) =  cell(2)*COSG
      omat(1,3) =  cell(3)*COSB
      omat(2,2) =  cell(2)*SING
      omat(2,3) = -cell(3)*SINB*COSAS
      omat(3,3) =  cell(3)*SINB*SINAS
c
c----+=================================================================
c
      return
      end
