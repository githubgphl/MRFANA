# MRFANA

MRFANA is a Fortran program to calculate quality metrics for unmerged
reflection data from X-Ray diffraction experiments.


## Getting Started

You will need
* a Fortran compiler (e.g. Intel's ifort or GNU gfortran)
* a copy of the [CCP4](http://www.ccp4.ac.uk/) libraries compatible with your Fortran compiler

To compile your own copy of CCP4 libraries:
```
git clone https://github.com/githubgphl/ccp4io.git
cd ccp4io
dir=`pwd`

cd libccp4
FFLAGS="-O2 -fallow-argument-mismatch" ./configure --prefix=$dir
make
make install

cd ../mmdb
g++ -O2 -c *.cpp
ar r ../lib/libmmdb2.a *.o
ranlib ../lib/libmmdb2.a
```

Compile with e.g.
```
cd src
gfortran -o mrfana *.f *.f90 -L/where/ever/ccp4/lib -lccp4f -lccp4c   # maybe add -static to get static binaries
````
or
```
cd src
ifort -o mrfana *.f *.f90 -L/where/ever/ccp4/lib -lccp4f -lccp4c
```



If you have a standard CCP4 installation (that comes with shared
libraries), the following commands might also work on Linux (but it
depends on your compiler version):
```
cd src
gfortran -o mrfana *.f *.f90 -L$CLIB -lccp4f -Wl,-rpath,$CLIB
```
or
```
cd src
ifort -o mrfana *.f *.f90 -L$CLIB -lccp4f -Wl,-rpath,$CLIB 
```


## Running

For help see
```
./mrfana -h
```

MRFANA can read unmerged reflection data in the following formats:
* XDS: [INTEGRATE.HKL](http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#INTEGRATE.HKL), [XDS_ASCII.HKL](http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html#XDS_ASCII.HKL) or [XSCALE output](http://xds.mpimf-heidelberg.mpg.de/html_doc/xscale_parameters.html#OUTPUT_FILE=)
* MTZ: [unmerged multi-record files](http://www.ccp4.ac.uk/html/mtzformat.html) e.g. from [AIMLESS](http://www.ccp4.ac.uk/html/aimless.html),
[POINTLESS](http://www.ccp4.ac.uk/html/pointless.html) or [SCALA](http://www.ccp4.ac.uk/html/scala.html)

A reflection file can be given on the command-line, e.g. using

```
./mrfana aimless_unmerged.mtz
```

or

```
./mrfana XDS_ASCII.HKL
```

Further details can be found [on the autoPROC
wiki](https://www.globalphasing.com/autoproc/wiki/index.cgi?MrfanaRunning),
especially concerning the different binning methods.

## Authors

* **Clemens Vonrhein**
* **Claus Flensburg**
* **Ian J. Tickle**
* **Gerard Bricogne**

See also

* [autoPROC](https://www.globalphasing.com/autoproc/): [Vonrhein, C., Flensburg, C., Keller, P., Sharff, A., Smart, O., Paciorek, W., Womack, T. & Bricogne, G. (2011). Data processing and analysis with the autoPROC toolbox. Acta Cryst. D67, 293-302.](https://scripts.iucr.org/cgi-bin/paper?ba5166)
* [STARANISO](http://staraniso.globalphasing.org/): Tickle, I.J., Flensburg, C., Keller, P., Paciorek, W., Sharff, A., Vonrhein, C., Bricogne, G. (2018). STARANISO. Cambridge, United Kingdom: Global Phasing Ltd.

## License

This project is licensed under the Mozilla Public License, v. 2.0 -
see the [LICENSE.md](LICENSE.md) file for details. © 2007-2019 Global Phasing Ltd.

## Acknowledgments

* Phil Evans, Kay Diederichs and Manfred Weiss for discussions about
  the finer details of the formulae for the various quality metrics

* Wolfgang Kabsch ([XDS](http://xds.mpimf-heidelberg.mpg.de/)) and
  [CCP4](http://www.ccp4.ac.uk/) for great software that is used
  extensively within the [autoPROC
  toolbox](https://www.globalphasing.com/autoproc/)

* All members of the [Global Phasing
  Consortium](https://www.globalphasing.com/)


