DESCRIPTION
===========

 MIDAS (MIDAS Isolayer Data Analysis Software)
 is a Python package primarily for processing
 gridded data stored in CF-compliant NetCDF format
 (http://cfconventions.org). 

 A handful of functions have been employed as class methods.
 For example: spatial interpolation between quadrilateral meshes 
 (FMS) and temporal interpolation between calendar dates (Datetime); conservative
 re-mapping in the vertical dimension (MOM6/ALE); spatial
 integration in one to three cartesian directions, e.g. 'X','XY' or 'XYZ'; 
 and temporal averaging (Datetime).
 
 
 MIDAS (Modular Isosurface Data Analysis System) was first developed by 
 Matthew Harrison 2011-2012 as an employee of NOAA in the 
 GFDL Oceans and Climate Group.    
 

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.



CONDA INSTALL
=============



	For Anaconda users (Python only, no extensions) 
	
	conda install -c https://conda.binstar.org/matthewharrison midas
	
	
PYTHON INSTALL 
==============


	cd /home/$USER/$install_dir

	git clone https://github.com/mjharriso/MIDAS.git

	# Alternatively using ssh
	# git clone git@github.com:mjharriso/MIDAS.git

	cd MIDAS

	# This will install the package under the current directory
	# subsequently, if you want to use in a session, the PYTHONPATH
	# environment variable must be set prior to invoking python, or,
	# alternatively , scripts must contain the following:
	# import sys; sys.path.append('foo_dir/local/lib/python')
	# where foo_dir is the MIDAS directory
	
	setenv PYTHONPATH `pwd`/local/lib/python
	
	# simple install. Pure Python.

	make   

	# OR With F90 external modules using gFortran
	
	make -f Makefile_gfortran

	# OR On GFDL HPCS using Intel
	module load python
	module load netcdf/4.2
	module load intel_compilers
	make -f Makefile_GFDL
	

EXAMPLES
========

	cd examples
	python contour_example.py # Fetches OpenDAP URL from NODC and plots with Matplotlib
	python hinterp_example.py # fms/horiz_interp does bi-linear interpolation from the original
	                          # 1-deg to a 5-deg grid with masking
	python hist.py            # volume-weighted histogram of salinity in the Indian Ocean
	python subtile.py         # use subtiling algorithm to calculate un-weighted
 		       		  # cell average bathymetry and roughness. Example along the Eastern US 

	
HOW TO OBTAIN DOCUMENTATION
=====
	
	# STRONGLY RECOMMEND using the Ipython interpreter.
	# If you are having trouble, submit a help desk ticket
	# or contact your system administrator.
	
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
				     
	                                      
	                                      

