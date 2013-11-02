INSTALL
=======

        

	cd /home/foo/install_dir
	git clone git@github.com:mjharriso/MIDAS.git
	cd MIDAS
	
	(cd fms;tar xvf fms_siena_201308.tar)

	# On GFDL HPCS
	(cd fms_build;./build_fms.csh)
	(python setup.py config_fc --f90flags="-i4 -r8 -DPY_SOLO" --fcompiler=intelem build)
	# On i686-linux
	# (cd fms_build;cp build_fms.csh tmp;sed -e 's/#set platform = linux/set platform = linux/' < tmp > tmp2;sed -e 's/set platform = gfdl_hpcs/#set platform = gfdl_hpcs/' < tmp2 > build_fms.csh;./build_fms.csh)
	# python setup.py config_fc --fcompiler=gfortran --f90flags="-fcray-pointer -fdefault-real-8 -ffixed-line-length-132 -ffree-line-length-0 -DPY_SOLO" build


	#  sudo python setup.py install # for root users
	python setup.py install --home=/home/foo/local 
	setenv PYTHONPATH=/home/foo/local/lib/python

USAGE
=====


	# Using the Ipython interpreter
	>>> import midas
	>>> midas.[Tab]   # complete listing of methods 
	>>> midas.state?  # Description of state instance
	>>> midas.state.[Tab] # available methods for state objects
	>>> midas.state.volume_integrals?  # spatial integration of fields 
	
	>>> midas.generic_rectgrid? # a generic rectangular grid description
	>>> import midas.hinterp
	>>> midas.hinterp.hinterp_mod.hinterp? # docstring for F90 interface
	
	# See examples directory
	>>> cd /home/foo/install_dir/MIDAS/examples
	>>> python contour_example.py # get WOA09 slice and make a contour plot
	>>> python hinterp_example.py # Regrid to 5 degree grid
	>>> python hist.py # Create a volume-weighted histogram of salinity
	                   # data in the Pacific Ocean
	>>> python test_grid_overlay.py # Calculate un-weighted cell averages
					# maximum/minimum and rms deviations
					# from a best-fit planar surface using
					# Scipy solver.
	
	
