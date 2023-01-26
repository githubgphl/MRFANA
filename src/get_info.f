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

c     get all the information from reflection file - including the
c     number of reflections.

      subroutine get_info(filtyp,nref,ncol,ncoltot,ranges,clabs,ilabs)
c
      implicit none
c
      include 'args.inc'
      include 'units.inc'
      include 'params.inc'
      include 'lookup.inc'
      include 'info.inc'
      include 'xds.inc'
      include 'dset.inc'
      include 'cell.inc'
c
      character*(*) filtyp
      character*30 clabs(maxcol)
      integer nref(2), ncol, ncoltot, ilabs(maxcol)
      double precision ranges(2,maxcol)
      real ranges_t(2,maxcol)
c
      character*1024 line, clibd, syminfo
      character*10 version
      character*1 ctyps(maxcol)
      integer jcol,i,j,k,ido, idataset
      double precision t(maxcol)
c
      integer nsymp, nsym
      character*1 ltype
      character*10 pgnam
      real rsym_t(4,4,192)
      logical lprint
c
      integer c2i
      external c2i
c
c     for MTZ files
c
      character*20 ename(150)
      character*5  etype(150)
      character*4  extn(150)
      logical lskip
      integer icount, iprint, sortx(5)
c
      real datcell_t(6,maxdatasets)
      real datwave_t(maxdatasets)
c
c----+=================================================================
c
      syminfo = ' '
      clibd = ' '
      call getenv('SYMINFO',syminfo)
      if (syminfo.eq.' ') then
        call getenv('CLIBD',clibd)
        if (clibd.eq.' ') then
          call getenv('CCP4',clibd)
          if (clibd.eq.' ') then
            call error_exit('unable to find working CCP4 installation')
          endif
          clibd = clibd(1:len_trim(clibd)) // '/lib/data'
        endif
        syminfo = clibd(1:len_trim(clibd)) // '/syminfo.lib'
      end if
      lskip = .true.
      icount = 0
      call csetnv('SYMINFO',syminfo,ename,etype,extn,icount,lskip)
c
      nref(1) = 0
      nref(2) = 0
      ncol = 0

c     --------------------------------------------------------
c     XDS
c     --------------------------------------------------------
      if (filtyp(1:4).eq.'XDS_') then
        rewind(lun)
        xds_item_iset = 0
        jcol = 0
 10     read(lun,'(a)',end=20) line
c
c       reflection
c
 15     if (line(1:1).ne.'!') then
          nref(1) = nref(1) + 1
          if (ncol.gt.0) then
            read(line,*) (t(i),i=1,ncol)
            do i = 1, ncol
              if (nref(1).eq.1) then
                ranges(1,ilabs(i)) = t(i)
                ranges(2,ilabs(i)) = t(i)
              else
                ranges(1,ilabs(i)) = min(ranges(1,ilabs(i)),t(i))
                ranges(2,ilabs(i)) = max(ranges(2,ilabs(i)),t(i))
              end if
            end do
          end if
c
c         header
c
        else
c                            1234567890123456789012345678901234567890
          if (line(1:37).eq.'!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=')
     +      then
            read(line(38:),*) ncol
            ncoltot = ncol
c
c           add columns that will be generated on the fly
c
            ncoltot = ncol + 1
            clabs(ncoltot) = 'M/ISYM'
            ilabs(ncoltot) = ncoltot
            ncoltot = ncoltot + 1
            clabs(ncoltot) = 'BATCH'
            ilabs(ncoltot) = ncoltot
            ncoltot = ncoltot + 1
            clabs(ncoltot) = 'FLAG'
            ilabs(ncoltot) = ncoltot
          else if (line(1:6).eq.'!ITEM_') then
            jcol = jcol + 1
            i = index(line,'=')
            clabs(jcol) = line(7:(i-1))
            read(line((i+1):),*) ilabs(jcol)
            if (line(2:11).eq.'ITEM_ISET=') then
              read(line(12:),*) xds_item_iset
            end if
          else if ((filtyp.eq.'XDS_INTEGRATE').and.
     +        (line(1:17).eq.'!H,K,L,IOBS,SIGMA')) then
            k=0
            ido=1
            do while (ido.eq.1)
              i=2
              do while (line(i:i).eq.' ')
                i = i + 1
              end do
              j=index(line(i:),',')
              do while (j.gt.0)
                k = k + 1
                clabs(k) = line(i:(i-1+j-1))
                ilabs(k) = k
                i=i+j
                j=index(line(i:),',')
                if ((j.eq.0).and.(line(i:).ne.' ')) then
                  j = len_trim(line(i:)) + 1
                end if
              end do
              if (line(len_trim(line):len_trim(line)).eq.',') then
                ido=1
              else
                ido=0
              end if
              read(lun,'(a)',end=20) line
            end do
            goto 15
          else if (line(1:20).eq.'!SPACE_GROUP_NUMBER=') then
            read(line(21:),*) symm
            call msymlb(99,symm,spgrn,pgnam,nsymp,nsym,rsym_t)
            call epsln(nsym,nsymp,rsym_t,0)
            lprint = .false.
            call asuset(spgrn,symm,pgnam,nsym,rsym_t,nsymp,nsyml,lprint)
          else if (line(1:21).eq.'!UNIT_CELL_CONSTANTS=') then
            read(line(22:),*) xds_cell
            read(line(22:),*) cell
          else if (line(1:18).eq.'!X-RAY_WAVELENGTH=') then
            read(line(19:),*) xds_wave
          else if ((filtyp.eq.'XDS_XSCALE').and.
     +        (line(13:29).eq.'X-RAY_WAVELENGTH=')) then
            read(line(30:39),*) xds_wave
c123456789012345678901234567890123456789012345678901234567890
c!COMPRISES THE FOLLOWING SCALED INPUT FILES:
c! ISET=   1 X-RAY_WAVELENGTH=   1.03961 (<0 if unknown)
c! ISET=  13 X-RAY_WAVELENGTH=   1.03961 (<0 if unknown)
          else if ((filtyp.eq.'XDS_XSCALE').and.
     +        (line(16:32).eq.'X-RAY_WAVELENGTH=')) then
c123456789012345678901234567890123456789012345678901234567890
c!COMPRISES THE FOLLOWING SCALED INPUT FILES:
c! ISET=      1 INPUT_FILE=195706_1_E1/XDS_ASCII.HKL
c! ISET=      1 X-RAY_WAVELENGTH=   0.97941 (<0 if unknown)
            read(line(33:42),*) xds_wave
          end if
        end if
        goto 10
 20     rewind(lun)
c
c     --------------------------------------------------------
c     MTZ
c     --------------------------------------------------------
      else if  (filtyp(1:3).eq.'MTZ') then
        call lrinfo (lun,version,ncol,nref(1),ranges_t)

        l_have_ranges = .true.
        do i = 1, ncol
          if (ranges_t(1,i).gt.ranges_t(2,i)) then
            l_have_ranges = .false.
            exit
          end if
        end do
        if (l_have_ranges) ranges = ranges_t

        call lrclab (lun,clabs,ctyps,ncol)
        ncoltot = ncol
        do i = 1, ncol
          ilabs(i) = i
        end do
c
c       pgnam = point group name
c
        call lrsymi(lun,nsymp,ltype,symm,spgrn,pgnam)
        call lrsymm(lun,nsym,rsym_t)
        lprint = .false.
c
c       unmerged MTZ files from AIMLESS 0.1.29 (at least) store the Laue
c       Group in the SYMINF record of the MTZ file - and the CCP4
c       library routines for extracting the Point Group therefore return
c       the wrong item! Therefore we need to call MSYMLB again to get the
c       correct value directly from syminfo.lib.
c
        call msymlb(99,symm,spgrn,pgnam,nsymp,nsym,rsym_t)
        call epsln(nsym,nsymp,rsym_t,0)
        call asuset(spgrn,symm,pgnam,nsym,rsym_t,nsymp,nsyml,lprint)
        ndatasets = maxdatasets
        call lridx(lun,pname,xname,dname,isets,datcell_t,datwave_t,
     +    ndatasets)
        datcell = datcell_t
        do idataset = 1, ndatasets
          call cell2omat(datcell(1,idataset),
     +                   datcell_orth_mat(1,1,idataset))
          write(stdout,'(/,a,i3,a,6f10.5)')
     +      ' Dataset',idataset,' cell = ',
     +      datcell(1,idataset),
     +      datcell(2,idataset),
     +      datcell(3,idataset),
     +      datcell(4,idataset),
     +      datcell(5,idataset),
     +      datcell(6,idataset)
          write(stdout,'(a,3f12.5)') '   Orthogonalisation matrix = ',
     +      (datcell_orth_mat(1,i,idataset),i=1,3)
          write(stdout,'(30x,3f12.5)')
     +      (datcell_orth_mat(2,i,idataset),i=1,3)
          write(stdout,'(30x,3f12.5)')
     +      (datcell_orth_mat(3,i,idataset),i=1,3)

          call matinv3(datcell_orth_mat(1,1,idataset),
     +                 datcell_orth_mat_inv(1,1,idataset))
          write(stdout,'(a,3f12.5)') '   Inverse of above         = ',
     +      (datcell_orth_mat_inv(1,i,idataset),i=1,3)
          write(stdout,'(30x,3f12.5)')
     +      (datcell_orth_mat_inv(2,i,idataset),i=1,3)
          write(stdout,'(30x,3f12.5)')
     +      (datcell_orth_mat_inv(3,i,idataset),i=1,3)
        end do
c
        cell(1) = datcell(1,1)
        cell(2) = datcell(2,1)
        cell(3) = datcell(3,1)
        cell(4) = datcell(4,1)
        cell(5) = datcell(5,1)
        cell(6) = datcell(6,1)
        call lrsort(lun,sortx)
        if ((abs(sortx(1)).ne.1).or.
     +      (abs(sortx(2)).ne.2).or.
     +      (abs(sortx(3)).ne.3)    ) then
          write(stdout,'(/,a)')
     +      ' MTZ file is not marked as being sorted on H, K, L - '
     +      //'will possibly need to do internal sorting.'
        end if
      end if
c
      ncoltot = ncoltot + 1
      clabs(ncoltot) = 'RESOL'
      ilabs(ncoltot) = ncoltot
      ncoltot = ncoltot + 1
      clabs(ncoltot) = 'RESO3'
      ilabs(ncoltot) = ncoltot
      ncoltot = ncoltot + 1
      clabs(ncoltot) = 'RUN'
      ilabs(ncoltot) = ncoltot
      ncoltot = ncoltot + 1
      clabs(ncoltot) = 'OUT'
      ilabs(ncoltot) = ncoltot
      ncoltot = ncoltot + 1
      clabs(ncoltot) = 'RAND'
      ilabs(ncoltot) = ncoltot
      ncoltot = ncoltot + 1
      clabs(ncoltot) = 'ELLIPSOID'
      ilabs(ncoltot) = ncoltot
c
c     call some CCP4 setup routines
c
      call setrsl(sngl(cell(1)),sngl(cell(2)),sngl(cell(3)),
     +            sngl(cell(4)),sngl(cell(5)),sngl(cell(6)))
      write(stdout,'(/,a,6f10.5)')
     +  ' Overall cell = ',cell
      call cell2omat(cell, cell_orth_mat)
      write(stdout,'(a,3f12.5)') ' Orthogonalisation matrix = ',
     +  (cell_orth_mat(1,i),i=1,3)
      write(stdout,'(28x,3f12.5)') (cell_orth_mat(2,i),i=1,3)
      write(stdout,'(28x,3f12.5)') (cell_orth_mat(3,i),i=1,3)

      call matinv3(cell_orth_mat,cell_orth_mat_inv)
      write(stdout,'(a,3f12.5)') ' Inverse of above         = ',
     +  (cell_orth_mat_inv(1,i),i=1,3)
      write(stdout,'(28x,3f12.5)') (cell_orth_mat_inv(2,i),i=1,3)
      write(stdout,'(28x,3f12.5)') (cell_orth_mat_inv(3,i),i=1,3)

      iprint = 0
      call centric(nsym,rsym_t,iprint)
c
c     setup lookup table
c
      icol(i_h    ) = c2i(filtyp,'H'           ,ncoltot,clabs,ilabs)
      icol(i_k    ) = c2i(filtyp,'K'           ,ncoltot,clabs,ilabs)
      icol(i_l    ) = c2i(filtyp,'L'           ,ncoltot,clabs,ilabs)
      icol(i_i    ) = c2i(filtyp,'I'           ,ncoltot,clabs,ilabs)
      icol(i_sigi ) = c2i(filtyp,'SIGI'        ,ncoltot,clabs,ilabs)
      icol(i_xobs ) = c2i(filtyp,'XOBS'        ,ncoltot,clabs,ilabs)
      icol(i_yobs ) = c2i(filtyp,'YOBS'        ,ncoltot,clabs,ilabs)
      icol(i_zobs ) = c2i(filtyp,'ZOBS'        ,ncoltot,clabs,ilabs)
      icol(i_lp   ) = c2i(filtyp,'LP'          ,ncoltot,clabs,ilabs)
      icol(i_xcal ) = c2i(filtyp,'XCAL'        ,ncoltot,clabs,ilabs)
      icol(i_ycal ) = c2i(filtyp,'YCAL'        ,ncoltot,clabs,ilabs)
      icol(i_zcal ) = c2i(filtyp,'ZCAL'        ,ncoltot,clabs,ilabs)
      icol(i_misym) = c2i(filtyp,'M/ISYM'      ,ncoltot,clabs,ilabs)
      icol(i_batch) = c2i(filtyp,'BATCH'       ,ncoltot,clabs,ilabs)
      icol(i_resol) = c2i(filtyp,'RESOL'       ,ncoltot,clabs,ilabs)
      icol(i_reso3) = c2i(filtyp,'RESO3'       ,ncoltot,clabs,ilabs)
      icol(i_flag ) = c2i(filtyp,'FLAG'        ,ncoltot,clabs,ilabs)
      icol(i_run  ) = c2i(filtyp,'RUN'         ,ncoltot,clabs,ilabs)
      icol(i_out  ) = c2i(filtyp,'OUT'         ,ncoltot,clabs,ilabs)
      icol(i_f    ) = c2i(filtyp,'F'           ,ncoltot,clabs,ilabs)
      icol(i_sigf ) = c2i(filtyp,'SIGF'        ,ncoltot,clabs,ilabs)
      icol(i_rand ) = c2i(filtyp,'RAND'        ,ncoltot,clabs,ilabs)
      icol(i_frac ) = c2i(filtyp,'FRACTIONCALC',ncoltot,clabs,ilabs)
      icol(i_iset ) = c2i(filtyp,'ISET'        ,ncoltot,clabs,ilabs)
      icol(i_ell  ) = c2i(filtyp,'ELLIPSOID'   ,ncoltot,clabs,ilabs)
c
c----+=================================================================
c
      return
      end
