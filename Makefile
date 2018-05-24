SHELL=bash

all: MIDAS

MIDAS:
	(conda create --name MIDAS)
	(git clone git@github.com:MJHarrison-GFDL/conda-recipes.git)
	(. activate MIDAS;cd conda-recipes/zlib;conda build .;conda install --use-local zlib)
	(. activate MIDAS;cd conda-recipes/hdf5;conda build .;conda install --use-local hdf5)
	(. activate MIDAS;cd conda-recipes/libnetcdf;conda build .;conda install --use-local libnetcdf)
	(. activate MIDAS;cd conda-recipes/libnetcdff;conda build .;conda install --use-local libnetcdff)
	(. activate MIDAS;pip install netCDF4)
	(. activate MIDAS;git clone git@github.com:mjharriso/MIDAS.git)
	(. activate MIDAS; cd MIDAS;git checkout dev/py36;. build.sh)





