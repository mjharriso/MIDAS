#!/bin/csh -fv


module use /usr/local/Modules/x64/modulefiles
module load intel_compilers
module load netcdf/4.1.2


#:/net2/mjh/local/bin
setenv PATH /net2/mjh/local/python-2.7.3/bin:/usr/local/x64/netcdf-4.1.2/bin:/usr/local/x64/hdf5-1.8.6/bin:/usr/local/x64/intel/Compiler/11.1/073/bin/intel64:/home/gfdl/bin:/usr/local/OSoverlay/bin:/bin:/usr/bin:/usr/local/bin:/usr/local/explorer/bin

unlimit


\rm *.o *.mod
ln -s ../fms_build/*.{o,mod} .
f2py  --verbose  --debug --fcompiler=intelem  --f90flags="-fPIC -O2 -g -i4 -r8"  -c -m hinterp_mod *.o hinterp.f90 -L/usr/local/x64/netcdf-4.1.2/lib/shared -L/usr/local/x64/hdf5-1.8.6/lib/shared -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz 

mv hinterp_mod.so ../.
