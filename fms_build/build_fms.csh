#!/bin/csh -v

set workdir = $cwd
set root = $cwd:h
set platform = intel
set mkmfTemplate = $root/fms/site/hpcs/$platform.mk
set npes = 1

set sharedir     = $root/fms/shared        
set includedir   = $sharedir/{include,mosaic}
set mppincludedir   = $sharedir/mpp/include
set MKMF         = $root/fms/bin/mkmf
#set cppDefs      = ("-Duse_netCDF -Duse_netCDF3 -Duse_libMPI" )
set cppDefs      = ("-Duse_netCDF -Duse_netCDF3" )
set LISTPATHS    = $root/fms/bin/list_paths

module purge
module load intel_compilers
module load netcdf/4.1.2


cd $root
\rm path_names*
$LISTPATHS $sharedir
mv -f path_names pt_orig
egrep -v "atmos_ocean_fluxes|coupler_types|coupler_util" pt_orig > path_names

cd $workdir
#\rm *.{o,mod}
$MKMF -m Makefile -a $root  -t $mkmfTemplate -p libfms.a -c "$cppDefs"  $root/path_names $includedir $mppincludedir  -L/usr/loca/lx64/netcdf-4.1.2/lib /usr/local/x64/netcdf-4.1.2/include

make NETCDF=3 libfms.a




if ( $status ) then
    unset echo
    echo ERROR: make failed.
    exit 1
else
    unset echo
    echo NOTE: make succeeded.
endif
