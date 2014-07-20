DESCRIPTION
===========

 MIDAS is a Python package primarily for processing
 gridded data stored in CF-compliant NetCDF format
 (http://cfconventions.org). Accompanying
 metadata information pertaining to variable attributes
 (units, standard name, ...) as well as the quadrilateral horizontal
 grid mesh on which the data reside and optional vertical cell boundaries
 between adjacent layers are stored (in a dictionary) for each field.

 To the extent possible, MIDAS follows CF conventions, however, the 
 encoding of grid cell and variable attributes is specific to this 
 application. The grid and field variable dictionaries are visible and malleable
 which is not without its disadvantages. It is entirely possible that
 this package will evolve to a more CF-like convention for class 
 methods. For example, instead of the current field.dict['units'] syntax
 a method-based approach, i.e. field.get_units(), would be desirable.
 

 A handful of functions have been employed which make use of 
 the stored information. For example: spatial interpolation between quad meshes 
 (FMS) and temporal interpolation between calendar dates (Datetime); spatial
 integration in one to three cartesian directions, e.g. 'X','XY' or 'XYZ'; 
 and temporal averaging (Datetime).
 
 
 MIDAS (Modular Isosurface Data Analysis System) was first developed by 
 Matthew Harrison 2011-2012 as an employee of NOAA in the 
 GFDL Oceans and Climate Group. The goal of this project is to produce a
 convenient and internally consistent environment in which to work with 
 gridded output from numerical models.    
 

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.



INSTALL 
=======



	For Anaconda users (Python only, no extensions) 
	
	conda install -c https://conda.binstar.org/matthewharrison midas
	
	
COMPLETE BUILD
==============

	cd /home/$USER/$install_dir

	git clone https://github.com/mjharriso/MIDAS.git

	# Alternatively using ssh
	# git clone git@github.com:mjharriso/MIDAS.git

	cd MIDAS

	# simple install. Pure Python.

	make   

	# With F90 external modules using gFortran
	#make -f Makefile_gfortran

	# On GFDL HPCS using Intel
	#module load python
	#module load netcdf/4.2
	#module load intel_compilers
	#make -f Makefile_GFDL
	

EXAMPLES
========

	cd examples
	python contour_example.py # Fetches OpenDAP URL from NODC and plots with Matplotlib
	python hinterp_example.py # fms/horiz_interp does bi-linear interpolation from the original
	                          # 1-deg to a 5-deg grid with masking
	python hist.py            # volume-weighted histogram of salinity in the Indian Ocean
	python subtile.py         # use subtiling algorithm to calculate un-weighted
 		       		  # cell average bathymetry and roughness. Example along the Eastern US 

	
USAGE
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
	>>> state.state.[Tab] # available methods for state objects
	>>> rectgrid.state.volume_integral?  # integrate scalars over the domain
				     # in 'X','Y','Z','XY','XZ',...
	                                      
	                                      

