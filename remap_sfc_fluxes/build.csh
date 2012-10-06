#!/bin/csh -v

module purge
module load ifort
module load icc
module load netcdf/4.1.2

setenv PATH /net2/mjh/local/python-2.7.3/bin:/net2/mjh/local/bin:/usr/local/x64/netcdf-4.1.2/bin:/usr/local/x64/hdf5-1.8.6/bin:/usr/local/x64/intel/Compiler/11.1/073/bin/intel64:/home/gfdl/bin:/usr/local/OSoverlay/bin:/bin:/usr/bin:/usr/local/bin:/usr/local/explorer/bin


f2py --fcompiler=intelem --f90flags="-fPIC -O0 -g"  -c -m remap_sfc_fluxes remap_sfc_fluxes.f90

mv remap_sfc_fluxes.so ../.
