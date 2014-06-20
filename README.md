INSTALL 
=======

	# (Tested on: Linux 3.8.0-30-generic i686 with GFortran and gcc version 4.7.3)
	
	# On i686-linux platform with NetCDF library compiled using gFortran
	# NOTE: 'extra_objects' are redundant in this environment.

        # From the shell
        
	cd /home/$USER/$install_dir
	
	git clone https://github.com/mjharriso/MIDAS.git
	# Alternatively using ssh
	# git clone git@github.com:mjharriso/MIDAS.git
	
	cd MIDAS
	
	(cd fms;tar xvf fms_tikal_201312.tar)
	
	# This compiles the FMS code using the proper F90/C compilers.
	
	(cd fms_build;cp build_fms.csh tmp;\
	sed -e 's/#set platform = linux/set platform = linux/' < tmp > tmp2;\
	sed -e 's/set platform = gfdl_hpcs/#set platform = gfdl_hpcs/' < tmp2 > build_fms.csh;\
	./build_fms.csh)
	
	(cd MOM6_ALE/build_ale;cp build_ale.csh tmp;\
	sed -e 's/#set platform = linux/set platform = linux/' < tmp > tmp2;\
	sed -e 's/set platform = gfdl_hpcs/#set platform = gfdl_hpcs/' < tmp2 > build_ale.csh;\
	./build_ale.csh)
	
	python setup.py config_fc --fcompiler=gfortran --f90flags="-fcray-pointer -fdefault-real-8 \
	-ffixed-line-length-132 -ffree-line-length-0 -DPY_SOLO" build

	#  sudo python setup.py install # for root users
	#  otherwise, install midas in your home directory
	#  and add this path to the environment variable PYTHONPATH. 
	# change --home=<path> to whatever you want.
	
	python setup.py install --home=/home/$USER/local 
	setenv PYTHONPATH /home/$USER/local/lib/python
	
	
	#  sudo python setup.py install # for root users
	


GFDL-HPCS INSTALL
=================

        # Log into a 64-bit workstation (preferred) or Analysis node 
        # From the shell
        
	cd /home/$USER/$install_dir
	
	git clone https://github.com/mjharriso/MIDAS.git
	# Alternatively using sshL
	# git clone git@github.com:mjharriso/MIDAS.git
	
	cd MIDAS
	
	(cd fms;tar xvf fms_tikal_201312.tar)
	
	#
	# NOTE:
	# =====
	# On GFDL HPCS (64-bit workstations and PPAN)
	# using installed python,netcdf and intel modules
	# The NetCDF libraries are compiled with ifort so we 
	# are obliged to build the fms/shared code and the f2py interfaces
	# using the Intel compiler. With the current environment, the
	# C object files are not found in the libfms library, although
	# this does work with gfortran (no need to include the 'extra_objects'
	# flags. in setup.py)
	#
	# This compiles the FMS code using the proper F90/C compilers.
	#
	#  OpenDAP currently disabled on Analysis Cluster. 
	#  Example scripts using OpenDap addresses
	#  will hang (unless timeouts are set)
	#
	#  With the exception of OpenDAP, the midas package is fully functional on both the
	#  GFDL 64-bit workstations and the Analysis cluster. 
	#
	#  Not yet tested on the PP cluster.
	#
	#  Recommend building on the workstations (where OpenDAP is available) if you have 
	#  enough hardware.
	#
	
	# If you have not already done so...
	
	module load python 
	module load netcdf/4.2
	module load intel_compilers
	
	
	# Build FMS
	
	(cd fms_build;./build_fms.csh)
	

	# Build ALE code
	
	(cd MOM6_ALE/build_ale;./build_ale.csh)
	
	# This builds the python interfaces to the underlying FMS and MOM6 codes
	
	(python setup.py config_fc --f90flags="-i4 -r8 -DPY_SOLO" --fcompiler=intelem build)

  	# This places the executables and modules where
  	# Python can find them. The environment variable is set here:
  	# Recommend installing your custom modules in the same location
  	
	python setup.py install --home=/home/$USER/local 
	setenv PYTHONPATH /home/$USER/local/lib/python
	

	
	
	



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
	
	ipython
	>>> import midas
	>>> midas.[Tab]   # complete listing of methods 
	>>> midas.state?  # Description of state instance
	>>> midas.state.[Tab] # available methods for state objects
	>>> midas.state.volume_integrals?  # integrate scalars over the domain
	                                   # in 'X','Y','Z','XY','XZ',...
	                                   
	
	>>> midas.rectgrid? # a generic rectangular grid description
				    # Can be read from a file or provided as
				    # 2-d lat/lon position arrays.
	>>> import midas.hinterp
	>>> midas.hinterp.hinterp_mod.hinterp? # docstring for Python interface
	                                       # to fms/shared/horiz_interp
	                                       # this module is not called directly in MIDAS,
	                                       # but via a class method.
	                                      
	                                      
UPDATING TO LATEST ON GitHub 
============================

	cd /home/foo/install_dir/MIDAS
	git status    
	git fetch
	# If you have updated your code, substitute a merge)
	git pull
	[python setup.py build from above]
	[python setup.py install from above]
	
