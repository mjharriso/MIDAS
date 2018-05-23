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

9. Install MIDAS in the root environment (libraries present)

```
git clone git@github.com:mjharriso/MIDAS.git
(cd MIDAS;git checkout dev/py36;. build.sh)
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

**MORE INFORMATION**

type

```
conda list
```

this returns your current environment which should look like this:

```
File Edit Options Buffers Tools Help
sqlalchemy                1.2.1            py36h14c3975_0
sqlite                    3.22.0               h1bed415_0
statsmodels               0.8.0            py36h8533d0b_0
sympy                     1.1.1            py36hc6d1c1c_0
tblib                     1.3.2            py36h34cf8b6_0
terminado                 0.8.1                    py36_1
testpath                  0.3.1            py36h8cadb63_0
tk                        8.6.7                hc745277_3
toolz                     0.9.0                    py36_0
tornado                   4.5.3                    py36_0
traitlets                 4.3.2            py36h674d592_0
typing                    3.6.2            py36h7da032a_0
unicodecsv                0.14.1           py36ha668878_0
unixodbc                  2.3.4                hc36303a_1
urllib3                   1.22             py36hbe7ace6_0
wcwidth                   0.1.7            py36hdf4376a_0
webencodings              0.5.1            py36h800622e_1
werkzeug                  0.14.1                   py36_0
wheel                     0.30.0           py36hfd4bba0_1
widgetsnbextension        3.1.0                    py36_0
wrapt                     1.10.11          py36h28b7045_0
xlrd                      1.1.0            py36h1db9f0c_1
xlsxwriter                1.0.2            py36h3de1aca_0
xlwt                      1.3.0            py36h7b00a1f_0
xz                        5.2.3                h55aa19d_2
yaml                      0.1.7                had09818_2
zict                      0.1.3            py36h3a3bf81_0
zlib                      1.2.11                        1    local

```