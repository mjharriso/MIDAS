INSTALL
=======
	# On GFDL HPCS
	(cd fms;tar xvf fms_siena_201308.tar)
	(cd fms_build;./build_fms.csh)
	(python setup.py config_fc --f90flags="-i4 -r8 -DPY_SOLO" --fcompiler=intelem build)
	# On i686-linux
	(cd fms;tar xvf fms_siena_201308.tar)
	(cf fms_build;cp build_fms.csh tmp;sed -e 's/#set platform = linux/set platform = linux/' < tmp > tmp;sed -e 's/set platform = gfdl_hpcs/#set platform = gfdl_hpcs/' < tmp > tmp;./build_fms.csh)
	python setup.py	config_fc --fcompiler=gnu --f90flags="-fcray-pointer -fdefault-real-8 -ffixed-line-length-132 -ffree-line-length-0 -DPY_SOLO" build
	python setup.py install

USAGE
=====


	# Using the Ipython interpreter
	>>> import midas
	>>> midas[Tab]   # complete listing of methods 
	>>> midas.generic_rectgrid? # docstring for this method
	>>> import midas.hinterp
	>>> midas.hinterp.hinterp_mod.hinterp? # docstring for F90 interface
	# See examples directory
