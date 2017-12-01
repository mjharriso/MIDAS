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


	Requirements:
	- bash
	- already installed mini-conda or anaconda (https://anaconda.org) and it is in your path.

	0. Configure Conda channels (~/.condarc):

	channels:
	- matthewharrison
	- ioos
	- scitools
	- conda-forge
	- defaults

	1. Create the enviromnent:

	conda env create -f MIDAS.yml

	2. Activate the enviroment to use MIDAS scripts

	source activate MIDAS

	3. Deactivate MIDAS environment when finished

	source deactivate

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
