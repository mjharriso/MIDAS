[![Build Status](https://travis-ci.org/mjharriso/MIDAS.svg?branch=master)](https://travis-ci.org/mjharriso/MIDAS)

# DESCRIPTION

 MIDAS (MIDAS Midas is Data Analysis Software)
 is a Python package primarily for processing
 gridded data stored in CF-compliant NetCDF/HDF5 format
 (http://cfconventions.org).

 * spatial interpolation between quadrilateral meshes
 * temporal interpolation between calendar dates (datetime)
 * conservative re-mapping in the vertical dimension (MOM6/ALE)
 * spatial integration/averaging with generalized horizontal/vertical coordinates
 * temporal averaging (Datetime).

 MIDAS was first developed by Matt Harrison as an employee of NOAA/GFDL
 and has been used for the generation of some of the realistic model configurations
 and post-run analysis scripts (https://github.com/NOAA-GFDL/MOM6-examples.git)


 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

# CONDA INSTALLATION

1. Download and install Anaconda

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
conda install conda-build
```

4. Install gfortran development libraries (if needed)

```
sudo apt-get install libgfortran-6-dev
```

optionally install mpich2 and associated libraries

```
(sudo apt-get install libmpich2-dev)
```

Or else if you do not have root privileges

```
(. activate;wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz;cd Downloads;tar xvf mpich-3.2.1.tar.gz;cd mpich-?3.2.1;./configure --enable-sha\
red --prefix=$CONDA_PREFIX;make; make install)
```


# MIDAS INSTALLATION

On linux platforms, simply type

```
make
```

This will take some time, since the Makefile will be downloading and compiling several large packages, like hdf5 and netcdf. Why are we compiling everything when there are pre-compiled binaries avaialble from the Anaconda cloud? Problems arise when linking c and c++ libraries to fortran APIs. There are often strict requirements for matching versions and glibc compatability which make using pre-compiled code unfeasible.

For platforms with multiple users, it is recommended that the compiled libraries and Python packages be made available through a local channel.

# MIDAS STEP-BY-STEP INSTALLATION

For best results, build the following libraries yourself - conda does not handle dependencies for linking c and c++ libraries to fortran APIs - consider yourself lucky if you can work with pre-compiled packages and associated libraries

```
(conda create --name MIDAS)
git clone git@github.com:MJHarrison-GFDL/conda-recipes.git
(. activate MIDAS; cd conda-recipes/zlib;conda build .;conda install --use-local zlib)
(. activate MIDAS; cd conda-recipes/hdf5;conda build .;conda install --use-local hdf5)
(. activate MIDAS; cd conda-recipes/libnetcdf;conda build .;conda install --use-local libnetcdf)
(. activate MIDAS; cd conda-recipes/libnetcdff;conda build .;conda install --use-local libnetcdff)
```

Install the netCDF4 python API

```
(. activate MIDAS;pip install netCDF4)
```

Install MIDAS

```
git clone git@github.com:mjharriso/MIDAS.git
(. activate MIDAS; cd MIDAS;git checkout dev/py36;. build.sh)
```

**TROUBLESHOOTING**

if you have a problem with libmkl missing:

```
(. activate MIDAS; conda install nomkl numpy scipy scikit-learn numexpr)
(. activate MIDAS; conda remove mkl mkl-service)
```

missing libcomm_err.so.3 at runtime?

```
(ln -s $CONDA_PREFIX/pkgs/krb5-1.14.6-0/lib/libcom_err.so.3 $CONDA_PREFIX/lib/.)
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

**MORE INFORMATION**

type

```
(. activate MIDAS;conda list)
```

this returns your current environment which should look something like this:

```
# Name                    Version                   Build  Channel
ca-certificates           2018.4.16                     0    conda-forge
curl                      7.60.0                        0    conda-forge
hdf5                      1.8.20                        0    local
krb5                      1.14.6                        0    conda-forge
libnetcdf                 4.4.1                         0    local
libnetcdff                4.4.4                         0    local
libssh2                   1.8.0                         2    conda-forge
openssl                   1.0.2o                        0    conda-forge
zlib                      1.2.11                        1    local

```