c     **********************************************************************
c
c      Copyright (c) 2015, 2019 Global Phasing Ltd.
c
c      This Source Code Form is subject to the terms of the Mozilla Public
c      License, v. 2.0. If a copy of the MPL was not distributed with this
c      file, You can obtain one at http://mozilla.org/MPL/2.0/.
c
c      Authors: Clemens Vonrhein, Claus Flensburg and Gerard Bricogne
c
c     **********************************************************************

c     standardised "table 1"

      subroutine write_table1(maxshell,mshell,iResoLimits,
     +  ResoLimits,Rval,NumIsigI,NumUnique,
     +  IsigImrg,Cmpltness,CmpltnessEll,Multiplicity,CC,SigAno,nrun)
c
      implicit none
c
      include 'units.inc'
      include 'xml.inc'
      include 'aniso.inc'
      include 'args.inc'
c
      integer maxshell, mshell, iResoLimits
      double precision ResoLimits(2,maxshell)
      double precision Rval(2,3,maxshell)
      integer NumIsigI(maxshell)
      integer NumUnique(2,maxshell)
      double precision IsigImrg(maxshell)
      double precision Cmpltness(2,maxshell)
      double precision CmpltnessEll(2,maxshell)
      double precision Multiplicity(2,maxshell)
      double precision CC(2,maxshell)
      double precision SigAno(maxshell)
      integer nrun
c
      integer iGeneral, iAnomalous
      parameter (iGeneral=1, iAnomalous=2)
c
      integer ierr, p(3), i, j
      character*30 t(3)
      character*12 CC123(3)
c
c----+=================================================================
c
      if (nrun.gt.1) then
        write(stdout,'(/,a,i3)')
     +    ' Number of RUNs (sweeps) contributing to this dataset = '
     +    ,nrun
      end if

      write(stdout,'(/,a77)')
     +  'Overall  InnerShell  OuterShell'
      write(stdout,'(a47,a30)')
     +  '---------------------------------------------',
     +  '------------------------------'
      write(stdout,'(a41,3f12.3)')
     +  '     Low resolution limit                ',
     +  ResoLimits(1,(mshell+1)),
     +  ResoLimits(1,1),
     +  ResoLimits(1,iResoLimits)
      write(stdout,'(a41,3f12.3)')
     +  '     High resolution limit               ',
     +  ResoLimits(2,(mshell+1)),
     +  ResoLimits(2,1),
     +  ResoLimits(2,iResoLimits)
      write(stdout,'(/)')

      write(stdout,'(a41,3f12.3)')
     +  '     Rmerge  (all I+ & I-)               ',
     +  Rval(iGeneral,1,(mshell+1)),
     +  Rval(iGeneral,1,1),
     +  Rval(iGeneral,1,iResoLimits)
      write(stdout,'(a41,3f12.3)')
     +  '     Rmerge  (within I+/I-)              ',
     +  Rval(iAnomalous,1,(mshell+1)),
     +  Rval(iAnomalous,1,1),
     +  Rval(iAnomalous,1,iResoLimits)

      write(stdout,'(a41,3f12.3)')
     +  '     Rmeas   (all I+ & I-)               ',
     +  Rval(iGeneral,2,(mshell+1)),
     +  Rval(iGeneral,2,1),
     +  Rval(iGeneral,2,iResoLimits)
      write(stdout,'(a41,3f12.3)')
     +  '     Rmeas   (within I+/I-)              ',
     +  Rval(iAnomalous,2,(mshell+1)),
     +  Rval(iAnomalous,2,1),
     +  Rval(iAnomalous,2,iResoLimits)

       write(stdout,'(a41,3f12.3)')
     +  '     Rpim    (all I+ & I-)               ',
     +  Rval(iGeneral,3,(mshell+1)),
     +  Rval(iGeneral,3,1),
     +  Rval(iGeneral,3,iResoLimits)
       write(stdout,'(a41,3f12.3)')
     +  '     Rpim    (within I+/I-)              ',
     +  Rval(iAnomalous,3,(mshell+1)),
     +  Rval(iAnomalous,3,1),
     +  Rval(iAnomalous,3,iResoLimits)

      write(stdout,'(a41,3i12)')  
     +  '     Total number of observations        ',
     +  NumIsigI((mshell+1)),
     +  NumIsigI(1),
     +  NumIsigI(iResoLimits)
      write(stdout,'(a41,3i12)') 
     +  '     Total number unique                 ',
     +  NumUnique(1,(mshell+1)),
     +  NumUnique(1,1),
     +  NumUnique(1,iResoLimits)
      write(stdout,'(a41,3f12.1)')
     +  '     Mean(I)/sd(I)                       ',
     +  IsigImrg((mshell+1)),
     +  IsigImrg(1),
     +  IsigImrg(iResoLimits)

      if (l_ell) then
        write(stdout,'(a41,3f12.1)')
     +    '     Completeness (spherical)            ',
     +    min(100.,100.*Cmpltness(iGeneral,(mshell+1))),
     +    min(100.,100.*Cmpltness(iGeneral,1)),
     +    min(100.,100.*Cmpltness(iGeneral,iResoLimits))
        write(stdout,'(a41,3f12.1)')
     +    '     Completeness (ellipsoidal)          ',
     +    min(100.,100.*CmpltnessEll(iGeneral,(mshell+1))),
     +    min(100.,100.*CmpltnessEll(iGeneral,1)),
     +    min(100.,100.*CmpltnessEll(iGeneral,iResoLimits))
      else
        write(stdout,'(a41,3f12.1)')
     +    '     Completeness                        ',
     +    min(100.,100.*Cmpltness(iGeneral,(mshell+1))),
     +    min(100.,100.*Cmpltness(iGeneral,1)),
     +    min(100.,100.*Cmpltness(iGeneral,iResoLimits))
      end if

      write(stdout,'(a41,3f12.1)')
     +  '     Multiplicity                        ',
     +  Multiplicity(iGeneral,(mshell+1)),
     +  Multiplicity(iGeneral,1),
     +  Multiplicity(iGeneral,iResoLimits)

      if (CC(iGeneral,(mshell+1)) .ge.1.01) then
        CC123(1) = '         NA '
      else
        write(CC123(1),'(f12.3)')
     +    max(-1.,min(1.,CC(iGeneral,(mshell+1))))
      endif
      if (CC(iGeneral,1)          .ge.1.01) then
        CC123(2) = '         NA '
      else
        write(CC123(2),'(f12.3)')
     +    max(-1.,min(1.,CC(iGeneral,1)))
      endif
      if (CC(iGeneral,iResoLimits).ge.1.01) then
        CC123(3) = '         NA '
      else
        write(CC123(3),'(f12.3)')
     +    max(-1.,min(1.,CC(iGeneral,iResoLimits)))
      end if
      write(stdout,'(a41,3a12)')
     +  '     CC(1/2)                             ',
     +  CC123
      write(stdout,'(/)')

      if (l_ell) then
        write(stdout,'(a41,3f12.1)')
     +    '     Anomalous completeness (spherical)  ',
     +    min(100.,100.*Cmpltness(iAnomalous,(mshell+1))),
     +    min(100.,100.*Cmpltness(iAnomalous,1)),
     +    min(100.,100.*Cmpltness(iAnomalous,iResoLimits))
        write(stdout,'(a41,3f12.1)')
     +    '     Anomalous completeness (ellipsoidal)',
     +    min(100.,100.*CmpltnessEll(iAnomalous,(mshell+1))),
     +    min(100.,100.*CmpltnessEll(iAnomalous,1)),
     +    min(100.,100.*CmpltnessEll(iAnomalous,iResoLimits))
      else
        write(stdout,'(a41,3f12.1)')
     +    '     Anomalous completeness              ',
     +    min(100.,100.*Cmpltness(iAnomalous,(mshell+1))),
     +    min(100.,100.*Cmpltness(iAnomalous,1)),
     +    min(100.,100.*Cmpltness(iAnomalous,iResoLimits))
      end if
      write(stdout,'(a41,3f12.1)')
     +  '     Anomalous multiplicity              ',
     +  Multiplicity(iAnomalous,(mshell+1)),
     +  Multiplicity(iAnomalous,1),
     +  Multiplicity(iAnomalous,iResoLimits)

      if (CC(iAnomalous,(mshell+1)) .ge.1.01) then
        CC123(1) = '         NA '
      else
        write(CC123(1),'(f12.3)')
     +    max(-1.,min(1.,CC(iAnomalous,(mshell+1))))
      endif
      if (CC(iAnomalous,1)          .ge.1.01) then
        CC123(2) = '         NA '
      else
        write(CC123(2),'(f12.3)')
     +    max(-1.,min(1.,CC(iAnomalous,1)))
      endif
      if (CC(iAnomalous,iResoLimits).ge.1.01) then
        CC123(3) = '         NA '
      else
        write(CC123(3),'(f12.3)')
     +    max(-1.,min(1.,CC(iAnomalous,iResoLimits)))
      end if
      write(stdout,'(a41,3a12)')
     +  '     CC(ano)                             ',
     +  CC123

      write(stdout,'(a41,3f12.3)')
     +  '     |DANO|/sd(DANO)                     ',
     +  SigAno((mshell+1)),
     +  SigAno(1),
     +  SigAno(iResoLimits)
      write(stdout,'(/)')

      if ((xml_base.ne.' ').and.(xml_file.ne.' ')) then
        open(71,file=xml_file,status='unknown',iostat=ierr)
        if (ierr.ne.0) then
        end if
        
        write(71,101)'  <AutoProc>'
        write(71,103)
     +    '    <spaceGroup>',
     +    TRIM(xml_sg),'</spaceGroup>'
        write(71,205)
     +    '    <wavelength>',
     +    xml_wave,'</wavelength>'
        write(71,203)
     +    '    <refinedCell_a>',
     +    xml_cell(1),'</refinedCell_a>'
        write(71,203)
     +    '    <refinedCell_b>',
     +    xml_cell(2),'</refinedCell_b>'
        write(71,203)
     +    '    <refinedCell_c>',
     +    xml_cell(3),'</refinedCell_c>'
        write(71,203)
     +    '    <refinedCell_alpha>',
     +    xml_cell(4),'</refinedCell_alpha>'
        write(71,203)
     +    '    <refinedCell_beta>',
     +    xml_cell(5),'</refinedCell_beta>'
        write(71,203)
     +    '    <refinedCell_gamma>',
     +    xml_cell(6),'</refinedCell_gamma>'
        write(71,101)
     +    '  </AutoProc>'
        write(71,101)
     +    '  <AutoProcScalingContainer>'
        write(71,101)
     +    '    <AutoProcScaling>'
        write(71,103)
     +    '      <recordTimeStamp>',
     +    xml_timestamp,'</recordTimeStamp>'
        if (l_ell) then
          write(71,101)
     +    '      <StaranisoEllipsoid>'
          write(71,101)
     +    '        <StaranisoEllipsoidRotationMatrix>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix11>',
     +                                 EllFit_mat_inv(1,1),
     +             '</StaranisoEllipsoidRotationMatrix11>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix12>',
     +                                 EllFit_mat_inv(1,2),
     +             '</StaranisoEllipsoidRotationMatrix12>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix13>',
     +                                 EllFit_mat_inv(1,3),
     +             '</StaranisoEllipsoidRotationMatrix13>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix21>',
     +                                 EllFit_mat_inv(2,1),
     +             '</StaranisoEllipsoidRotationMatrix21>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix22>',
     +                                 EllFit_mat_inv(2,2),
     +             '</StaranisoEllipsoidRotationMatrix22>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix23>',
     +                                 EllFit_mat_inv(2,3),
     +             '</StaranisoEllipsoidRotationMatrix23>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix31>',
     +                                 EllFit_mat_inv(3,1),
     +             '</StaranisoEllipsoidRotationMatrix31>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix32>',
     +                                 EllFit_mat_inv(3,2),
     +             '</StaranisoEllipsoidRotationMatrix32>'
          write(71,205)
     +    '          <StaranisoEllipsoidRotationMatrix33>',
     +                                 EllFit_mat_inv(3,3),
     +             '</StaranisoEllipsoidRotationMatrix33>'
          write(71,101)
     +    '        </StaranisoEllipsoidRotationMatrix>'
          write(71,101)
     +    '        <StaranisoEllipsoidEigenvalues>'
          write(71,205)
     +    '          <StaranisoEllipsoidEigenvalue1>',
     +                                 EllFit_res(1),
     +             '</StaranisoEllipsoidEigenvalue1>'
          write(71,205)
     +    '          <StaranisoEllipsoidEigenvalue2>',
     +                                 EllFit_res(2),
     +             '</StaranisoEllipsoidEigenvalue2>'
          write(71,205)
     +    '          <StaranisoEllipsoidEigenvalue3>',
     +                                 EllFit_res(3),
     +             '</StaranisoEllipsoidEigenvalue3>'
          write(71,101)
     +    '        </StaranisoEllipsoidEigenvalues>'
          write(71,101)
     +    '      </StaranisoEllipsoid>'
        end if
        write(71,101)'    </AutoProcScaling>'

        t(1) = 'overall'
        p(1) = mshell + 1
        t(2) = 'innerShell'
        p(2) = 1
        t(3) = 'outerShell'
        p(3) = iResoLimits
        do i = 1, 3
          j = p(i)
          write(71,101)'    <AutoProcScalingStatistics>'
          write(71,103)
     +      '      <scalingStatisticsType>',
     +      trim(t(i)),'</scalingStatisticsType>'

          write(71,203)
     +      '      <resolutionLimitLow>',
     +      ResoLimits(1,j),'</resolutionLimitLow>'
          write(71,203)
     +      '      <resolutionLimitHigh>',
     +      ResoLimits(2,j),'</resolutionLimitHigh>'
          write(71,203)
     +      '      <rMerge>',
     +      Rval(iGeneral,1,j),'</rMerge>'
          write(71,203)
     +      '      <rMeasWithinIPlusIMinus>',
     +      Rval(iAnomalous,2,j),'</rMeasWithinIPlusIMinus>'
          write(71,203)
     +      '      <rMeasAllIPlusIMinus>',
     +      Rval(iGeneral,2,j),'</rMeasAllIPlusIMinus>'
          write(71,203)
     +      '      <rPimWithinIPlusIMinus>',
     +      Rval(iAnomalous,3,j),'</rPimWithinIPlusIMinus>'
          write(71,203)
     +      '      <rPimAllIPlusIMinus>',
     +      Rval(iGeneral,3,j),'</rPimAllIPlusIMinus>'
          write(71,301)
     +      '      <nTotalObservations>',
     +      NumIsigI(j),'</nTotalObservations>'
          write(71,301)
     +      '      <nTotalUniqueObservations>',
     +      NumUnique(1,j),'</nTotalUniqueObservations>'
          write(71,201)
     +      '      <meanIOverSigI>',
     +      IsigImrg(j),'</meanIOverSigI>'

          if (l_ell) then
            write(71,201)
     +        '      <completenessSpherical>',
     +        min(100.,100.*Cmpltness(iGeneral,j)),
     +        '</completenessSpherical>'
            write(71,201)
     +        '      <completenessEllipsoidal>',
     +        min(100.,100.*CmpltnessEll(iGeneral,j)),
     +        '</completenessEllipsoidal>'
            write(71,201)
     +        '      <completeness>',
     +        min(100.,100.*CmpltnessEll(iGeneral,j)),
     +        '</completeness>'
          else
            write(71,201)
     +        '      <completeness>',
     +        min(100.,100.*Cmpltness(iGeneral,j)),
     +        '</completeness>'
          endif

          write(71,201)
     +      '      <multiplicity>',
     +      Multiplicity(iGeneral,j),'</multiplicity>'


          if (CC(iGeneral,j) .ge.1.01) then
            CC123(1) = '0.0'
          else
            write(CC123(1),'(f12.3)')
     +        max(-1.,min(1.,CC(iGeneral,j)))
          endif
          write(71,'(3a)')
     +      '      <ccHalf>',
     +      trim(adjustl(CC123(1))),'</ccHalf>'

          if (l_ell) then
            write(71,201)
     +        '      <anomalousCompletenessSpherical>',
     +        min(100.,100.*Cmpltness(iAnomalous,j)),
     +        '</anomalousCompletenessSpherical>'
            write(71,201)
     +        '      <anomalousCompletenessEllipsoidal>',
     +        min(100.,100.*CmpltnessEll(iAnomalous,j)),
     +        '</anomalousCompletenessEllipsoidal>'
            write(71,201)
     +        '      <anomalousCompleteness>',
     +        min(100.,100.*CmpltnessEll(iAnomalous,j)),
     +        '</anomalousCompleteness>'
          else
            write(71,201)
     +        '      <anomalousCompleteness>',
     +        min(100.,100.*Cmpltness(iAnomalous,j)),
     +        '</anomalousCompleteness>'
          end if
          write(71,201)
     +      '      <anomalousMultiplicity>',
     +      Multiplicity(iAnomalous,j),
     +      '</anomalousMultiplicity>'

          if (CC(iAnomalous,j) .ge.1.01) then
            CC123(1) = '0.0'
          else
            write(CC123(1),'(f12.3)')
     +        max(-1.,min(1.,CC(iAnomalous,j)))
          endif
          write(71,'(3a)')
     +      '      <ccAnomalous>',
     +      trim(adjustl(CC123(1))),'</ccAnomalous>'

          write(71,203)
     +      '      <DanoOverSigDano>',
     +      SigAno(j),'</DanoOverSigDano>'
          write(71,101)'    </AutoProcScalingStatistics>'
        end do
        write(71,101)'  </AutoProcScalingContainer>'
        
      end if
c
 101  format(1a)
 103  format(3a)
 201  format(a,f0.1,a)
 203  format(a,f0.3,a)
 205  format(a,f0.5,a)
 301  format(a,i0,a)
c
c----+=================================================================
c
      return
      end
