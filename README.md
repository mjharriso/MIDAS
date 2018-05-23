[![Build Status](https://travis-ci.org/mjharriso/MIDAS.svg?branch=master)](https://travis-ci.org/mjharriso/MIDAS)

# DESCRIPTION


 MIDAS (MIDAS Midas is Data Analysis Software)
 is a Python package primarily for processing
 gridded data stored in CF-compliant NetCDF/HDF5 format
 (http://cfconventions.org).

 A handful of functions have been employed as class methods.
 For example: spatial interpolation between quadrilateral meshes using
 the FMS code developed at NOAA/GFDL; temporal interpolation between calendar
 dates (datetime); conservative re-mapping in the vertical dimension using
 MOM6/ALE developed at Princeton and NOAA/GFDL; spatial
 integration in one to three cartesian directions, e.g. 'X','XY' or 'XYZ';
 and temporal averaging (Datetime).


 MIDAS was first developed by Matthew Harrison 2011-2012 as an employee of NOAA in the
 GFDL Oceans and Sea Ice Processes Group.

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.



# CONDA INSTALL

> wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
> #(cd Downloads;./Miniconda3-latest-Linux-x86_64.sh)
> wget https://repo.anaconda.com/archive/Anaconda3-5.1.0-Linux-x86_64.sh
> (cd Downloads;./Anaconda3-5.1.0-Linux-x86_64.sh)
> source ~/.bashrc
> conda update conda

* Activate Conda and setup the default (root) environment

> . activate
> sudo apt-get instal libgfortran-6-dev
* Which version of gcc do I have? 
> gcc -v
...
gcc version 6.3.0 20170516 (Debian 6.3.0-18+deb9u1)


* For best results, build these libraries yourself - conda does not handle
* dependencies for linking c and c++ libraries to fortran APIs - consider yourself
* lucky if you can work with pre-compiled packages and associated libraries

> conda install conda-build
> git clone git@github.com:MJHarrison-GFDL/conda-recipes.git
> (cd conda-recipes/zlib;conda build .;conda install --use-local zlib)
> (cd conda-recipes/hdf5;conda build .;conda install --use-local hdf5)
> (cd conda-recipes/libnetcdf;conda build .;conda install --use-local libnetcdf)
>(cd conda-recipes/libnetcdff;conda build .;conda install --use-local libnetcdff)
> $nc-config --libs


* Optionally install with mpich2
* If you have root privleges


> (sudo apt-get install libmpich2-dev)

*Or else if you do not have root privleges

> (wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz;cd Downloads;tar xvf mpich-3.2.1.tar.gz;cd mpich-?3.2.1;./configure --enable-sha\
red --prefix=/home/$USER/anaconda3;make; make install)


* missing libcomm_err.so.3 at runtime?

> ln -s /home/$USER/anaconda3/pkgs/krb5-1.14.6-0/lib/libcom_err.so.3 /home/$USER/anaconda3/lib/.

# install the netCDF4 python API

> pip install netCDF4

* Setup a custom environment for MIDAS

> git clone git@github.com:mjharriso/MIDAS.git
> conda create --name MIDAS
> . deactivate
> . activate MIDAS
> (cd MIDAS;git checkout dev/py36;. build.sh)



* TROUBLESHOOTING

if you have a problem with libmkl missing:

> conda install nomkl numpy scipy scikit-learn numexpr
> conda remove mkl mkl-service



* EXAMPLES


> cd examples
> source activate MIDAS
> python contour_example.py # Fetches OpenDAP URL from NODC and plots with Matplotlib
> python hinterp_example.py # fms/horiz_interp does bi-linear interpolation from the original
>	                          # 1-deg to a 5-deg grid with masking
> python hist.py            # volume-weighted histogram of salinity in the Indian Ocean
> python subtile.py         # use subtiling algorithm to calculate un-weighted
> 		       		  # cell average bathymetry and roughness. Example along the Eastern US
> source deactivate


* HOW TO OBTAIN DOCUMENTATION

> ipython
>>> import midas.rectgrid as rectgrid
>>> rectgrid.[Tab]   # complete listing of methods
>>> rectgrid.quadmesh       # a generic rectangular grid description
>				    # Can be read from a file or provided as
>				    # 2-d lat/lon position arrays.
 >>> rectgrid.state?  # Description of state instance
>>> rectgrid.state.state.[Tab] # available methods for state objects
>>> rectgrid.state.volume_integral?  # integrate scalars over the domain
>				     # in 'X','Y','Z','XY' or 'XYZ'
>>> rectgrid.state?? # View the code
