[![Build Status](https://travis-ci.org/mjharriso/MIDAS.svg?branch=master)](https://travis-ci.org/mjharriso/MIDAS)

DESCRIPTION
===========

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



CONDA INSTALL
=============

############################
# INSTALL ANACONDA
############################

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
#(cd Downloads;./Miniconda3-latest-Linux-x86_64.sh)
wget https://repo.anaconda.com/archive/Anaconda3-5.1.0-Linux-x86_64.sh
(cd Downloads;./Anaconda3-5.1.0-Linux-x86_64.sh)
source ~/.bashrc
conda update conda

##############################
# Activate Conda and setup the default (root) environment
##############################
. activate
sudo apt-get instal libgfortran-6-dev
## Which version of gcc do I have? ##
gcc -v
Using built-in specs.
COLLECT_GCC=gcc
COLLECT_LTO_WRAPPER=/usr/lib/gcc/x86_64-linux-gnu/6/lto-wrapper
Target: x86_64-linux-gnu
Configured with: ../src/configure -v --with-pkgversion='Debian 6.3.0-18+deb9u1' --with-bugurl=file:///usr/share/doc/gcc-6/README.Bugs --enable-langua\
ges=c,ada,c++,java,go,d,fortran,objc,obj-c++ --prefix=/usr --program-suffix=-6 --program-prefix=x86_64-linux-gnu- --enable-shared --enable-linker-bui\
ld-id --libexecdir=/usr/lib --without-included-gettext --enable-threads=posix --libdir=/usr/lib --enable-nls --with-sysroot=/ --enable-clocale=gnu --\
enable-libstdcxx-debug --enable-libstdcxx-time=yes --with-default-libstdcxx-abi=new --enable-gnu-unique-object --disable-vtable-verify --enable-libmp\
x --enable-plugin --enable-default-pie --with-system-zlib --disable-browser-plugin --enable-java-awt=gtk --enable-gtk-cairo --with-java-home=/usr/lib\
/jvm/java-1.5.0-gcj-6-amd64/jre --enable-java-home --with-jvm-root-dir=/usr/lib/jvm/java-1.5.0-gcj-6-amd64 --with-jvm-jar-dir=/usr/lib/jvm-exports/ja\
va-1.5.0-gcj-6-amd64 --with-arch-directory=amd64 --with-ecj-jar=/usr/share/java/eclipse-ecj.jar --with-target-system-zlib --enable-objc-gc=auto --ena\
ble-multiarch --with-arch-32=i686 --with-abi=m64 --with-multilib-list=m32,m64,mx32 --enable-multilib --with-tune=generic --enable-checking=release --\
build=x86_64-linux-gnu --host=x86_64-linux-gnu --target=x86_64-linux-gnu
Thread model: posix
gcc version 6.3.0 20170516 (Debian 6.3.0-18+deb9u1)


############################
#For best results, build these libraries yourself - conda does not handle
#dependencies for linking c and c++ libraries to fortran APIs - consider yourself
#lucky if you can work with pre-compiled packages and associated libraries
############################
conda install conda-build
git clone git@github.com:MJHarrison-GFDL/conda-recipes.git
(cd conda-recipes/zlib;conda build .;conda install --use-local zlib)
(cd conda-recipes/hdf5;conda build .;conda install --use-local hdf5)
(cd conda-recipes/libnetcdf;conda build .;conda install --use-local libnetcdf)
(cd conda-recipes/libnetcdff;conda build .;conda install --use-local libnetcdff)
$nc-config --libs

############################
#Optionally install with mpich2
#If you have root privleges
############################

(sudo apt-get install libmpich2-dev)
############################
#Or else if you do not have root privleges
############################
(wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz;cd Downloads;tar xvf mpich-3.2.1.tar.gz;cd mpich-3.2.1;./configure --enable-sha\
red --prefix=/home/$USER/anaconda3;make; make install)


#missing libcomm_err.so.3 at runtime?
ln -s /home/$USER/anaconda3/pkgs/krb5-1.14.6-0/lib/libcom_err.so.3 /home/$USER/anaconda3/lib/.

#install the netCDF4 python API
pip install netCDF4


############################
# Setup a custom environment for MIDAS
############################


git clone git@github.com:mjharriso/MIDAS.git
conda create --name MIDAS
. deactivate
. activate MIDAS
(cd MIDAS;git checkout dev/py36;. build.sh)



TROUBLESHOOTING
===============

	if you have a problem with libmkl missing:

	conda install nomkl numpy scipy scikit-learn numexpr
	conda remove mkl mkl-service



EXAMPLES
========

	cd examples
	source activate MIDAS
	python contour_example.py # Fetches OpenDAP URL from NODC and plots with Matplotlib
	python hinterp_example.py # fms/horiz_interp does bi-linear interpolation from the original
	                          # 1-deg to a 5-deg grid with masking
	python hist.py            # volume-weighted histogram of salinity in the Indian Ocean
	python subtile.py         # use subtiling algorithm to calculate un-weighted
 		       		  # cell average bathymetry and roughness. Example along the Eastern US
	source deactivate


HOW TO OBTAIN DOCUMENTATION
===========================


	ipython
	>>> import midas.rectgrid as rectgrid
	>>> rectgrid.[Tab]   # complete listing of methods
	>>> rectgrid.quadmesh       # a generic rectangular grid description
				    # Can be read from a file or provided as
				    # 2-d lat/lon position arrays.
	>>> rectgrid.state?  # Description of state instance
	>>> rectgrid.state.state.[Tab] # available methods for state objects
	>>> rectgrid.state.volume_integral?  # integrate scalars over the domain
				     # in 'X','Y','Z','XY' or 'XYZ'
	>>> rectgrid.state?? # View the code
