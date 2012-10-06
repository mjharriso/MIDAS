#!/bin/csh -v

module load ifort
module load icc
module load netcdf/4.1.2

setenv PATH /net2/mjh/local/python-2.7.3/bin:/net2/mjh/local/bin:/usr/local/x64/netcdf-4.1.2/bin:/usr/local/x64/hdf5-1.8.6/bin:/usr/local/x64/intel/Compiler/11.1/073/bin/intel64:/home/gfdl/bin:/usr/local/OSoverlay/bin:/bin:/usr/bin:/usr/local/bin:/usr/local/explorer/bin


\rm *.o *.mod
ln -s ../fms_build/*.{o,mod} .
f2py   --fcompiler=intelem  --f90flags="-fPIC -O0 -g"  -c -m hinterp_mod *.o hinterp.f90 -L/usr/local/x64/netcdf-4.1.2/lib/shared -L/usr/local/x64/hdf5-1.8.6/lib/shared -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz  -L/net2/mjh/local/lib -lmpich -L/net2/mjh/local/lib -lmpl

mv hinterp_mod.so ../.
