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

c     read command-line arguments as well as (option) arguments on
c     stdin:
c
c       RESOlution RUN <irun> <reslow> <reshigh>
c       NAME       RUN <irun> [...] DATAset <dname>

      subroutine get_args()
c
      implicit none
c
      include 'params.inc'
      include 'args.inc'
      include 'info.inc'
      include 'units.inc'
      include 'runs.inc'
      include 'reso.inc'
      include 'ice.inc'
      include 'xml.inc'
      include 'aniso.inc'
c
c----+=================================================================
c
      integer numarg, i, irun, istdin, j
      integer  COMMAND_ARGUMENT_COUNT
      double precision t, r1, r2, IceRingsUsedMargin
      double precision td
c
      integer, dimension(:), allocatable :: seed
      integer seed_size
c
      character*1024 arg, line, line0
      character*10240 long_line
      character*128 word, t1, t2, t3
c
c----+=================================================================
c
c     defaults
c
      iverb  = 0
      nshell = 0
      iset_nshell = 0
      nrefbin     = 1000
      L_BinningByEqualNumber = .false.
      N_BinningByEqualNumber = 0
      nrefbin_run =  200
      L_BinningByEqualNumber_run = .false.
      N_BinningByEqualNumber_run = 0
      lshell3 = .true.
      filnam = ' '
      xml_base = ' '
      reslow = 999999.9d0
      reshigh = 0.0d0
      reshigh_ell = 999999.9d0
      reslow_ell = 0.0d0
      current_reshigh = 0.0
      istdin = 0
      consider_ice = .true.
      nIce = 0
      iwide_format = 0
      check_part_min = 0.01d0
      check_part_max = 1.99d0

      l_skip_empty_runs = .false.

c     Rmerge Rmeas Rpim
c     denominator: 1 = summing over Ij
c                  2 =              <I>
      calc_type_r(1,1) = 2
      calc_type_r(1,2)  = 2
      calc_type_r(1,3)   = 2
c     <I> : 1 = unweighted
c           2 = weighted
      calc_type_r(2,1) = 2
      calc_type_r(2,2)  = 2
      calc_type_r(2,3)   = 2
c
      do irun = 1, maxruns
        run_reslim(1,irun) = 999999.0d0
        run_reslim(2,irun) = 0.01d0
        run_num(irun)      = 0
        run_i2n(irun)      = 0
        run_batlim(1,irun) = 9999999
        run_batlim(2,irun) = -999999
        run_nbatexcl(irun) = 0
        run_set(irun) = 0
      end do

c     to do with anisotropy
      l_ell = .false.
      do i = 1, 3
        do j = 1, 3
          EllFit_mat(i,j) = 0.0d0
        end do
        EllFit_res(i) = 0.0d0
      end do
c
c     default random number seed
c
      rand_seed = -1
c
      numarg = COMMAND_ARGUMENT_COUNT()
      i = 0
      do while (i.le.numarg)
        if (i.eq.0) then
          call GET_COMMAND_ARGUMENT ( i, pgmnam )
        else
          call GET_COMMAND_ARGUMENT ( i, arg )
c
          if (iverb.gt.0) then
            write(stdout,'(2a)') ' processing argument   ',
     +        arg(1:len_trim(arg))
          end if
c
c     -v : increase verbosity
c
          if (arg.eq.'-v') then
            iverb = iverb + 1
c
c     -ice/-noice/-icer <reslow> <reshigh>
c
          else if (arg.eq.'-ice') then
            consider_ice = .true.
          else if (arg.eq.'-noice') then
            consider_ice = .false.
          else if (arg.eq.'-icer') then
            consider_ice = .true.
            nIce = nIce + 1
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg )
            read(arg,*) IceRingsUsed(1,nIce)
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg )
            read(arg,*) IceRingsUsed(2,nIce)
            if (IceRingsUsed(1,nIce).lt.IceRingsUsed(2,nIce)) then
              IceRingsUsed(2,nIce)=IceRingsUsed(1,nIce)
              read(arg,*) IceRingsUsed(1,nIce)
            end if
c           safety margin
            IceRingsUsedMargin = 0.1*(IceRingsUsed(1,nIce)
     +        -IceRingsUsed(2,nIce))
c           report
            write(t1,'(f10.3)') IceRingsUsed(1,nIce)
            write(t2,'(f10.3)') IceRingsUsed(2,nIce)
            write(t3,'(f10.4)') IceRingsUsedMargin
            write(stdout,'(7a)') ' NOTE: taking ice-ring (',
     +        trim(adjustl(t1)),'-',
     +        trim(adjustl(t2)),
     +        ' A) into account with safety margin +- ',
     +        trim(adjustl(t3)),' A.'
            IceRingsUsed(1,nIce) = IceRingsUsed(1,nIce) +
     +        IceRingsUsedMargin
            IceRingsUsed(2,nIce) = IceRingsUsed(2,nIce) -
     +        IceRingsUsedMargin
c
c     -wide-format
c
          else if (arg.eq.'-wide-format') then
            iwide_format = 1
c
c     -skip_empty_runs
c
          else if (arg.eq.'-skip_empty_runs') then
            l_skip_empty_runs = .true.
c
c     -n <nshell> : number of bins (equal volume); a negative
c           value will do the 'standard' binning on 1/d**2
c
          else if (arg.eq.'-n') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg )
            read(arg,*) nshell
            if (nshell.lt.0) then
              nshell = -nshell
              lshell3 = .false.
            end if
            if (nshell.ge.maxshell) then
              write(stdout,'(a,i4)')
     +          ' increase MAXSHELL (and recompile) from currently'
     +          ,maxshell
              call error_exit('too many resolution shells requested')
            end if
            iset_nshell = 1
c
c     -nref <nrefbin> : number of measured reflections per bin
c
          else if (arg.eq.'-nref') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg )
            read(arg,*) nrefbin
            if (nrefbin.lt.0) then
              L_BinningByEqualNumber = .true.
              nrefbin=-nrefbin
              lshell3 = .false.
            end if
c
c     -nrefrun <nrefbin_run> : number of measured reflections per
c           bin in runs
c
          else if (arg.eq.'-nrefrun') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg )
            read(arg,*) nrefbin_run
            if (nrefbin_run.lt.0) then
              L_BinningByEqualNumber_run = .true.
              nrefbin_run=-nrefbin_run
            end if
c
c     -bins <r1>,<r2>,...,<rN> : explicit bin limits (reso >= rN)
c           with final bin going from rN to high resolution limit
c
          else if (arg.eq.'-bins') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, long_line )
            L_BinningByEqualNumber     = .true.
            L_BinningByEqualNumber_run = .true.
            N_BinningByEqualNumber     = 0
            N_BinningByEqualNumber_run = 0

            j = index(long_line,',')

            do while (j.ne.0)

              N_BinningByEqualNumber = N_BinningByEqualNumber + 1

              if (N_BinningByEqualNumber.gt.(maxshell-2)) then
                write(stdout,'(a)')
     +            ' increase MAXSHELL (and recompile)'//
     +            ' or specify fewer explicit bins'
                call error_exit('too many explicit bins requested')
              end if

              read(long_line(1:(j-1)),*)
     +          BinningByEqualNumber(N_BinningByEqualNumber)
              N_BinningByEqualNumber_run =
     +          N_BinningByEqualNumber_run + 1
              BinningByEqualNumber_run(N_BinningByEqualNumber_run) =
     +          BinningByEqualNumber(N_BinningByEqualNumber)

              if (iverb.gt.0) then
                if (N_BinningByEqualNumber.eq.1) then
                  write(stdout,'(/)')
                end if
                write(stdout,'(a,i3,1x,a,1x,f8.4,a)')
     +            ' setting explicit resolution bin ',
     +            N_BinningByEqualNumber,': R <=',
     +            BinningByEqualNumber(N_BinningByEqualNumber),' A'
              end if

              long_line = long_line((j+1):)
              j = index(long_line,',')
              if ((j.eq.0).and.(long_line.ne.' ')) then
                j = 1024
              end if

            end do

            if (iverb.gt.0) then
              if (N_BinningByEqualNumber.gt.1) then
                write(stdout,'(a,i3,1x,a,1x,f8.4,a,/)')
     +            ' setting explicit resolution bin ',
     +            (N_BinningByEqualNumber+1),': R > ',
     +            BinningByEqualNumber(N_BinningByEqualNumber),' A'
              end if
            end if
c
c     -xml <some.xml>
c
          else if (arg.eq.'-xml') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, xml_base )
c
c     -check_part <check_part_min> <check_part_max>
c
          else if (arg.eq.'-check_part') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg ) 
            read(arg,*) check_part_min
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg ) 
            read(arg,*) check_part_max
            if (check_part_min.gt.check_part_max) then
              td=check_part_min
              check_part_min=check_part_max
              check_part_max=td
            end if
c
c     -calc_type_r <1|2|3> <integer> <integer>
c
c           rmerge rmeas rpim
c
          else if (arg.eq.'-calc_type_r') then

            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg ) 
            read(arg,*) j

            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg ) 
            read(arg,*) calc_type_r(1,j)

            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg ) 
            read(arg,*) calc_type_r(2,j)
c
c     eigenvector matrix: eigenvalue/vector list
c
c     -ell <Res1> <EV1(1)> <EV1(2)> <EV1(3)> \
c          <Res2> <EV2(1)> <EV2(2)> <EV2(3)> \
c          <Res3> <EV3(1)> <EV3(2)> <EV3(3)> \
c                
c
          else if (arg.eq.'-ell') then
            do j = 1, 3
              i = i + 1
              call GET_COMMAND_ARGUMENT ( i, arg ) 
              read(arg,*) EllFit_res(j)
              EllFit_res_inv(j) = 1.0d0/EllFit_res(j)
              i = i + 1
              call GET_COMMAND_ARGUMENT ( i, arg )
              read(arg,*) EllFit_mat(1,j)
              i = i + 1
              call GET_COMMAND_ARGUMENT ( i, arg )
              read(arg,*) EllFit_mat(2,j)
              i = i + 1
              call GET_COMMAND_ARGUMENT ( i, arg )
              read(arg,*) EllFit_mat(3,j)
            end do
            call matinv3(EllFit_mat,EllFit_mat_inv)
            l_ell = .true.
            write(stdout,'(/,a,3F10.5)')
     +            '   Eigenvalue matrix         = ',
     +                                     EllFit_mat(1,1),
     +                                     EllFit_mat(1,2),
     +                                     EllFit_mat(1,3)
            write(stdout,'(31x,3F10.5)')   EllFit_mat(2,1),
     +                                     EllFit_mat(2,2),
     +                                     EllFit_mat(2,3)
            write(stdout,'(31x,3F10.5)')   EllFit_mat(3,1),
     +                                     EllFit_mat(3,2),
     +                                     EllFit_mat(3,3)
            write(stdout,'(/,a,3F10.5)')
     +            '   Ellipsoid rotation matrix = ',
     +                                     EllFit_mat_inv(1,1),
     +                                     EllFit_mat_inv(1,2),
     +                                     EllFit_mat_inv(1,3)
            write(stdout,'(31x,3F10.5)')   EllFit_mat_inv(2,1),
     +                                     EllFit_mat_inv(2,2),
     +                                     EllFit_mat_inv(2,3)
            write(stdout,'(31x,3F10.5)')   EllFit_mat_inv(3,1),
     +                                     EllFit_mat_inv(3,2),
     +                                     EllFit_mat_inv(3,3)
            write(stdout,'(/,a,3F10.5,/)')
     +            '   Eigenvalues               = ',
     +                                     EllFit_res
c
c     -r <reslow> <reshigh>
c
          else if (arg.eq.'-r') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg ) 
            read(arg,*) reslow
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg ) 
            read(arg,*) reshigh
            if (reshigh.gt.reslow) then
              t=reshigh
              reshigh=reslow
              reslow=t
            end if
            do irun = 1, maxruns
              run_reslim(1,irun) = reslow
              run_reslim(2,irun) = reshigh
            end do
c
c     -i : to trigger reading from stdin
c
          else if (arg.eq.'-i') then
            istdin = 1
c
c     -seed <rand_seed>
c
          else if (arg.eq.'-seed') then
            i = i + 1
            call GET_COMMAND_ARGUMENT ( i, arg )
            read(arg,*) rand_seed
c
c     -h
c
          else if (arg.eq.'-h') then
            write(stdout,'(/)')
            write(stdout,'(8(a,/))')
     +        ' USAGE: mrfana [-h] [-i] [-n <nshell>]'//
     +        '                            \',
     +        '               [-nref <nrefbin>]'//
     +        ' [-nrefrun <nrefbin_run>]         \',
     +        '               [-bins <r1>,<r2>,...,<rN>]'//
     +        ' [-r <reslow> <reshigh>] \',
     +        '               [-ell'//
     +        ' <Res1> <EV1(1)> <EV1(2)> <EV1(3)>'//
     +        '            \',
     +        '                     '//
     +        '<Res2> <EV2(1)> <EV2(2)> <EV2(3)>'//
     +        '            \',
     +        '                     '//
     +        '<Res3> <EV3(1)> <EV3(2)> <EV3(3)>]'//
     +        '           \',
     +        '               <MTZ file|XDS file>'

            write(stdout,'(/,2a)')
     +        '        -h                        : ',
     +        'show this help message'

            write(stdout,'(/,2a)')
     +        '        -i                        : ',
     +        'read additional cards from standard input'

            write(stdout,'(/,2a)')
     +        '        -n <nshell>               : ',
     +        'number of bins (equal volume); a negative'
            write(stdout,'(2a)')
     +        '                                    ',
     +        'value will do the ''standard'' binning on 1/d**2'

            write(stdout,'(/,2a)')
     +        '        -nref <nrefbin>           : ',
     +        'number of measured reflections per bin'
            write(stdout,'(2a,i6,a)')
     +        '                                    ',
     +        '(currently = ',nrefbin,'); a negative'
            write(stdout,'(2a)')
     +        '                                    ',
     +        'value will trigger binning by equal numbers'

            write(stdout,'(/,2a)')
     +        '        -nrefrun <nrefbin_run>    : ',
     +        'number of measured reflections per bin and run'
            write(stdout,'(2a,i6,a)')
     +        '                                    ',
     +        '(currently = ',nrefbin_run,'); a negative'
            write(stdout,'(2a)')
     +        '                                    ',
     +        'value will trigger binning by equal numbers'


            write(stdout,'(/,2a)')
     +        '        -bins <r1>,<r2>,...,<rN>  : ',
     +        'explicit bin limits (reso >= rN) with final bin'
            write(stdout,'(2a)')
     +        '                                    ',
     +        'going from rN to high resolution limit'

            write(stdout,'(//,2a)')
     +        '        -r <reslow> <reshigh>     : ',
     +        'resolution limits'


            write(stdout,'(//,2a)')
     +        '        -ell ...                  : ',
     +        'eigenvector matrix (eigenvalue/vector list)'

            write(stdout,'(/)')
c
            stop
c
          else if (arg(1:1).eq.'-') then
            call error_exit('unknown argument given = '//arg)
          else
            if (filnam.ne.' ') then
              call error_exit('more than one file name given')
            end if
            filnam = arg
          end if
        end if
        i = i + 1
      end do
c
c     checks
c
      if (filnam.eq.' ')
     +  call error_exit('no filename given on command-line')
c
      if (.not.L_BinningByEqualNumber) then
        if (nshell.gt.0) then
          if (lshell3) then
            write(stdout,'(/,a,i3,a)') ' binning will be in ',nshell
     +        ,' shells of equal volume'
          else
            write(stdout,'(/,a,i3,a)') ' binning will be in ',nshell
     +        ,' shells'
          end if
        end if
      end if

c     now read stdin
      numrun = 0
c
      if (istdin.eq.0) goto 100
c
      write(stdout,'(/,a,/)') ' expecting keywords from stdin:'
c
 10   read(stdin,'(a)',end=100) line
c
      line = trim(adjustl(line))
      if (line.eq.' ') goto 10
      line0 = line
c
c     comment
c
      if ((line(1:1).eq."!").or.(line(1:1).eq."#")) goto 10
c
      write(stdout,'(2a)') ' Data line--- ',line0(1:len_trim(line0))
c
      call getword(line,word)
c
c     RUN <irun> [SET <iset>] BATCH <ibatch1> TO <ibatch2> EXCLUDE <xbat1> <xbat2> ... <xbatN>
c
      if (word(1:3).eq.'RUN') then
        call rmword(line)
        call getword(line,word)
        read(word,*) irun
        numrun = numrun + 1
        if ((numrun.gt.maxruns).or.(irun.gt.maxruns)) then
          write(stdout,*) ' numrun, irun, MAXRUNS = ',
     +      numrun,irun,maxruns
          call error_exit('too many runs defined on RUN card')
        end if
c
c       numrun is the current counter/stack
c
c       run_num stores the specification (RUN <i>)
c
c       run_i2n goes back from <i> to stack - we need that in later
c       cards that also use RUN <i> syntax.
c
        run_num(numrun) = irun
        run_i2n(irun) = numrun
        run_batlim(1,numrun) = 9999999
        run_batlim(2,numrun) = -999999
        run_nbatexcl(numrun) = 0
        run_dname(numrun) = ' '
        run_set(numrun) = 0
        call rmword(line)
        call getword(line,word)
        if (word(1:3).eq.'SET') then
          call rmword(line)
          call getword(line,word)
          read(word,*) run_set(numrun)
          call rmword(line)
          call getword(line,word)
        end if
        if (word(1:4).eq.'BATC') then
          call rmword(line)
          call getword(line,word)
          read(word,*) run_batlim(1,numrun)
          call rmword(line)
          call getword(line,word)
          if (word(1:2).ne.'TO') then
            call error_exit('unknown sub-keyword on RUN card')
          end if
          call rmword(line)
          call getword(line,word)
          read(word,*) run_batlim(2,numrun)
          call rmword(line)
          call getword(line,word)
          if (word(1:4).eq.'EXCL') then
            call rmword(line)
            call getword(line,word)
            run_nbatexcl(numrun) = 0
            do while (word(1:1).ne.' ')
              run_nbatexcl(numrun) = run_nbatexcl(numrun) + 1
              if (run_nbatexcl(numrun).gt.maxrunexcl) then
                write(stdout,*) ' MAXRUNEXCL = ',maxrunexcl
                call error_exit('too many EXCLUDed batches on RUN card')
              end if
              read(word,*) run_batexcl(run_nbatexcl(numrun),numrun)
              call rmword(line)
              call getword(line,word)
            end do
          end if
        else
          call error_exit('unknown sub-keyword on RUN card')
        end if
c
c     RESOlution RUN <irun> <reslow> <reshigh>
c
      else if (word(1:4).eq.'RESO') then
        call rmword(line)
        call getword(line,word)
        if (word(1:3).ne.'RUN') then
          call error_exit('unknown sub-keyword on RESOLUTION card')
        end if
        call rmword(line)
        call getword(line,word)
        read(word,*) irun
        if (run_i2n(irun).eq.0) then
          call error_exit('undefined RUN used on RESOLUTION card')
        end if
        if (run_batlim(1,run_i2n(irun)).gt.
     +      run_batlim(2,run_i2n(irun))    ) then
          call error_exit('undefined RUN used on RESOLUTION card')
        end if
        call rmword(line)
        read(line,*) r1,r2
        run_reslim(1,run_i2n(irun)) =
     +    min(r1,run_reslim(1,run_i2n(irun)))
        run_reslim(2,run_i2n(irun)) =
     +    max(r2,run_reslim(2,run_i2n(irun)))
c
c     NAME RUN <irun> ... DATAset <dname>
c
      else if (word(1:4).eq.'NAME') then
        call rmword(line)
        call getword(line,word)
        if (word(1:3).ne.'RUN') then
          call error_exit('unknown sub-keyword on NAME card')
        end if
        call rmword(line)
        call getword(line,word)
        read(word,*) irun
        if (run_i2n(irun).eq.0) then
          call error_exit('undefined RUN used on NAME card')
        end if
        if (run_batlim(1,run_i2n(irun)).gt.
     +      run_batlim(2,run_i2n(irun))    ) then
          call error_exit('undefined RUN used on NAME card')
        end if
        call getword(line,word)
        do while ((word(1:4).ne.'DATA').and.(word.ne.' '))
          call rmword(line)
          call getword(line,word)
        end do
        if (word(1:4).ne.'DATA') then
          call error_exit('DATAset not found on NAME card')
        end if
        call rmword(line)
        call getword(line,run_dname(run_i2n(irun)))
c
      else if (word(1:3).eq.'END') then
        goto 100
c
c     anything else
c
      else
        call error_exit('unknown keyword')
      end if
c
      goto 10
c
 100  continue
c
cv
c     this part should only be done if all reflections are also part of
c     a RUN definition!
cv
      do irun = 1, numrun
        if (irun.eq.1) then
          r1 = run_reslim(1,irun)
          r2 = run_reslim(2,irun)
        else
          r1 = max(r1,run_reslim(1,irun))
          r2 = min(r2,run_reslim(2,irun))
        end if
      end do
      if (numrun.gt.0) then
        reslow  = min(reslow,r1)
        reshigh = max(reshigh,r2)
      end if
c
c     setting up random number generator
      if (rand_seed.lt.0) then
        call system_clock (rand_seed)
        write(stdout,'(/,a,i10)') ' Setting random seed to ',rand_seed
        write(stdout,'(a,/)')
     +    ' You can set this also with the -seed flag'
      end if
      seed_size = 1
      CALL RANDOM_SEED(SIZE = seed_size)
      allocate(seed(1:seed_size))
      seed = rand_seed
      CALL RANDOM_SEED(PUT = seed)
c
      if (xml_base.ne.' ') then
        i = len_trim(xml_base)
        if (xml_base((i-3):i).eq.'.xml') then
          xml_base((i-3):i) = '    '
        end if
      end if
         
c
      return
      end
