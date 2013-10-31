INSTALL
=======
	
	# On i686-linux
	(cd fms;tar xvf fms_siena_201308.tar)
	(cf fms_build;./build_fms.csh)
	python setup.py	build	
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