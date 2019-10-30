[Peck, A., Sunden, F., Andrews, L. D., Pande, V. S., & Herschlag,
D. (2016). Tungstate as a transition state analog for catalysis by
alkaline phosphatase. Journal of molecular biology, 428(13),
2758-2768.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6169531/)

To download the [diffraction data](https://data.sbgrid.org/dataset/456/) for [5C66](https://www.rcsb.org/structure/5C66):
```
  rsync -av rsync://data.sbgrid.org/10.15785/SBGRID/456 .
```
Processing with [autoPROC](https://www.globalphasing.com/autoproc/) ([Vonrhein et al, 2011](https://scripts.iucr.org/cgi-bin/paper?ba5166)) using
```
  process -I 456 -d autoPROC.01 | tee autoPROC.01.lis
```
then gives several [unmerged reflection files](https://www.globalphasing.com/autoproc/wiki/index.cgi?MTZforStaranisoServer) (in subdirectory
autoPROC.01):
```
  INTEGRATE.HKL
  XDS_ASCII.HKL
  aimless_alldata_unmerged.mtz
  aimless_unmerged.mtz
```
These can be analysed with MRFANA via e.g.
```
  mrfana -nref -1000 autoPROC.01/aimless_alldata_unmerged.mtz | tee aimless_alldata_unmerged.mrfana.log
  mrfana -n 20       autoPROC.01/aimless_unmerged.mtz         | tee aimless_unmerged.mrfana.log
  mrfana -n -20      autoPROC.01/XDS_ASCII.HKL                | tee XDS_ASCII.mrfana.log
```

Of course, the full autoPROC processing analysis contains additional information and plots - see [summary_inlined.html](http://htmlpreview.github.com/?https://github.com/githubgphl/MRFANA/blob/master/examples/5C66/autoPROC.01/summary_inlined.html) and e.g. [report_staraniso.pdf](autoPROC.01/report_staraniso.pdf). With the correct information about the anisotropy analysis (since the crystal diffracts anisotropically to 2.3A in two and 1.6A in the third direction), one can also run
```
  mrfana -nref -1000 autoPROC.01/aimless_alldata_unmerged.mtz \
    -ell 2.33 1.0 0.0 0.0 2.33 0.0 1.0 0.0 1.625 0.0 0.0 1.0 | tee aimless_alldata_unmerged.mrfana-aniso.log
```
