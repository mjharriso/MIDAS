#!/bin/csh -v

set workdir = $cwd
set root = $cwd:h:h
#set platform = linux
set platform = gfdl_hpcs

set maledir=$root/MOM6_ALE/src


set mkmfTemplate = $root/fms/site/$platform/mkmf.template
set npes = 1



set MKMF         = $root/fms/bin/mkmf
set cppDefs      = ("" )
set LISTPATHS    = $root/fms/bin/list_paths


if ($platform == 'gfdl_hpcs') then
  module load netcdf/4.2
  module load intel_compilers
endif




\rm *.{o,mod}
\rm path_names*
$LISTPATHS $maledir

$MKMF -m Makefile -a $workdir  -t $mkmfTemplate -p libale.a -c "$cppDefs"  path_names 

make DEBUG=0 libale.a


if ( $status ) then
    unset echo
    echo ERROR: make failed.
    exit 1
else
    unset echo
    echo NOTE: make succeeded.
endif

