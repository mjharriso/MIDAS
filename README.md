[![Build Status](https://travis-ci.org/mjharriso/MIDAS.svg?branch=master)](https://travis-ci.org/mjharriso/MIDAS)

# DESCRIPTION

 MIDAS (MIDAS Midas is Data Analysis Software)
 is a Python package primarily for processing
 gridded data stored in CF-compliant NetCDF/HDF5 format
 (http://cfconventions.org).

 * spatial interpolation between quadrilateral meshes
 * temporal interpolation between calendar dates (datetime)
 * conservative re-mapping in the vertical dimension using MOM6/ALE
 * spatial integration/averaging
 * temporal averaging (Datetime).

 MIDAS was first developed by Matthew Harrison as an employee of NOAA in the
 GFDL Oceans and Sea Ice Processes Group.

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

# CONDA INSTALL

1. Download and install miniconda

```
(wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh;./Miniconda3-latest-Linux-x86_64.sh)
```

or, alternatively full Anaconda

```
(wget https://repo.anaconda.com/archive/Anaconda3-5.1.0-Linux-x86_64.sh;./Anaconda3-5.1.0-Linux-x86_64.sh)
```

2. Update your existing shell to add conda to your path

```
source ~/.bashrc
```

3. Now update conda

```
conda update conda
```

4. Install libgfortran (if needed)

```
sudo apt-get install libgfortran-6-dev
```

5. Activate Conda and setup the default (root) environment

```
source activate
```

6. Build zlib/hdf5/libnetcdf/libnetcdff

For best results, build these libraries yourself - conda does not handle dependencies for linking c and c++ libraries to fortran APIs - consider yourself lucky if you can work with pre-compiled packages and associated libraries

```
conda install conda-build
git clone git@github.com:MJHarrison-GFDL/conda-recipes.git
(cd conda-recipes/zlib;conda build .;conda install --use-local zlib)
(cd conda-recipes/hdf5;conda build .;conda install --use-local hdf5)
(cd conda-recipes/libnetcdf;conda build .;conda install --use-local libnetcdf)
(cd conda-recipes/libnetcdff;conda build .;conda install --use-local libnetcdff)
```

7. Optionally install with mpich2 if you have root privileges

```
(sudo apt-get install libmpich2-dev)
```

Or else if you do not have root privileges

```
(wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz;cd Downloads;tar xvf mpich-3.2.1.tar.gz;cd mpich-?3.2.1;./configure --enable-sha\
red --prefix=$CONDA_PREFIX;make; make install)
```


8. install the netCDF4 python API

```
pip install netCDF4
```

9. Setup a custom environment for MIDAS

```
git clone git@github.com:mjharriso/MIDAS.git
conda create --name MIDAS
. deactivate
(. activate MIDAS;cd MIDAS;git checkout dev/py36;. build.sh)
```

**TROUBLESHOOTING**

if you have a problem with libmkl missing:

```
conda install nomkl numpy scipy scikit-learn numexpr
conda remove mkl mkl-service
```

missing libcomm_err.so.3 at runtime?

```
ln -s CONDA_PREFIX/pkgs/krb5-1.14.6-0/lib/libcom_err.so.3 $CONDA_PREFIX/lib/.
```


**EXAMPLES**

```
cd examples
source activate MIDAS
python contour_example.py # Fetches OpenDAP URL from NODC and plots with Matplotlib
python hinterp_example.py # fms/horiz_interp does bi-linear interpolation from the original
                 	  # 1-deg to a 5-deg grid with masking
python hist.py            # volume-weighted histogram of salinity in the Indian Ocean
python subtile.py         # use subtiling algorithm to calculate un-weighted
 	       		  # cell average bathymetry and roughness. Example along the Eastern US
source deactivate
```

**HOW TO OBTAIN DOCUMENTATION**

```
ipython
import midas.rectgrid as rectgrid
rectgrid.[Tab]   # complete listing of methods
rectgrid.quadmesh       # a generic rectangular grid description
				    # Can be read from a file or provided as
				    # 2-d lat/lon position arrays.
rectgrid.state?  # Description of state instance
rectgrid.state.state.[Tab] # available methods for state objects
rectgrid.state.volume_integral?  # integrate scalars over the domain
				     # in 'X','Y','Z','XY' or 'XYZ'
rectgrid.state?? # View the code
```
