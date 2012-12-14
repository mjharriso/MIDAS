#!/bin/csh -v


#module use /usr/local/Modules/x64/modulefiles
module load ifort
module load icc
module load netcdf/4.1.2


#:/net2/mjh/local/bin
setenv PATH /net2/mjh/local/python-2.7.2/bin:${PATH}

unlimit


\rm *.o *.mod
ln -s ../fms_build/*.{o,mod} .
f2py  --verbose  --debug --fcompiler=intelem  --f90flags="-fPIC -O2 -g -i4 -r8"  -c -m hinterp_mod *.o hinterp.f90 -L/usr/local/x64/netcdf-4.1.2/lib/shared -L/usr/local/x64/hdf5-1.8.6/lib/shared -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz 

mv hinterp_mod.so ../.
