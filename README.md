INSTALL
=======

        

	cd /home/foo/install_dir
	git clone git@github.com:mjharriso/MIDAS.git
	# or if SSH does not work, try HTTPS
	# git clone https://github.com/mjharriso/MIDAS.git
	cd MIDAS
	
	(cd fms;tar xvf fms_siena_201308.tar)
	
	
	# On GFDL HPCS (64-bit workstations of Analysis Cluster)
	# Using installed python,netcdf and intel modules
	# The NetCDF libraries are compiled with ifort so we 
	# are obliged to build the fms/shared code and the f2py interfaces
	# using the Intel compiler. With the current environment, the
	# C object files are not found in the libfms library, although
	# this does work with gfortran (no need to include the 'extra_objects'
	# flags.
	
	(cd fms_build;./build_fms.csh)
	
	module load python 
	module load netcdf/4.2
	module load intel_compilers
	
	(python setup.py config_fc --f90flags="-i4 -r8 -DPY_SOLO" --fcompiler=intelem build)

	#  sudo python setup.py install # for root users
	#  otherwise, install midas in your home directory
	#  and add this path to the environment variable PYTHONPATH. 
	#  This is the currently recommended method at GFDL.  
	#  I recommend that you do NOT load the "analysis_du_jour"
	#  module which can easily introduce conflicts. The recommendation
	#  is to load python/netcdf/intel_compilers as shown above prior
	#  to initiating your batch or interactive session.
	
	python setup.py install --home=/home/foo/local 
	setenv PYTHONPATH /home/foo/local/lib/python
	
	
	# On i686-linux platform with NetCDF library compiled using gFortran
	# NOTE: 'extra_objects' are not required but currently exist.
	
	(cd fms_build;cp build_fms.csh tmp;\
	sed -e 's/#set platform = linux/set platform = linux/' < tmp > tmp2;\
	sed -e 's/set platform = gfdl_hpcs/#set platform = gfdl_hpcs/' < tmp2 > build_fms.csh;\
	./build_fms.csh)
	
	python setup.py config_fc --fcompiler=gfortran --f90flags="-fcray-pointer -fdefault-real-8 -ffixed-line-length-132 -ffree-line-length-0 -DPY_SOLO" build

	#  sudo python setup.py install # for root users
	#  otherwise, install midas in your home directory
	#  and add this path to the environment variable PYTHONPATH. 
	#  This is the currently recommended method at GFDL.  
	
	python setup.py install --home=/home/foo/local 
	setenv PYTHONPATH /home/foo/local/lib/python



EXAMPLES
========

	cd examples
	python contour_example.py # Fetches OpenDAP URL from NODC and plots with Matplotlib
	python hinterp_example.py # fms/horiz_interp does bi-linear interpolation from the original
	                          # 1-deg to a 5-deg grid with masking
	python hist.py            # volume-weighted histogram of salinity in the Indian Ocean
	
USAGE
=====
	

	# STRONGLY RECOMMEND using the Ipython interpreter.
	# If you are having trouble, submit a help desk ticket
	# or contact your system administrator.
	
	>ipython
	>>> import midas
	>>> midas.[Tab]   # complete listing of methods 
	>>> midas.state?  # Description of state instance
	>>> midas.state.[Tab] # available methods for state objects
	>>> midas.state.volume_integrals?  # integrate scalars over the domain
	                                   # in 'X','Y','Z','XY','XZ',...
	                                   
	
	>>> midas.generic_rectgrid? # a generic rectangular grid description
				    # Can be read from a file or provided as
				    # 2-d lat/lon position arrays.
	>>> import midas.hinterp
	>>> midas.hinterp.hinterp_mod.hinterp? # docstring for Python interface
	                                       # to fms/shared/horiz_interp
	                                       # this module is not called directly in MIDAS,
	                                       # but via a class method.
	                                       

