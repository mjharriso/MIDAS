"""
==============================

 MIDAS (Modular Iso-surface Data Analysis System)
 was first developed by Matthew Harrison
 starting in 2011 as an employee of NOAA in the 
 GFDL Oceans and Ice Sheet Processes Group. 

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

===============================
"""
from numpy.distutils.core import setup,Extension

doclines = __doc__.split("\n")


remap_sfc_fluxes = Extension(name = 'remap_sfc_fluxes',
                sources = ['remap_sfc_fluxes/remap_sfc_fluxes.f90'])

hinterp = Extension(name = 'fms_hinterp',
                include_dirs = ['fms_build'],
                library_dirs = ['fms_build'],
                libraries = ['fms','netcdf','netcdff'],
                extra_objects = ['fms_build/*.o'],
                sources = ['hinterp/hinterp.f90'])



vertmap_GOLD = Extension(name = 'vertmap_GOLD',
                sources = ['vertmap_GOLD/vertmap_GOLD.F90']
                )

vertmap_ALE = Extension(name = 'vertmap_ALE',
                include_dirs = ['MOM6_ALE/build_ale'],
                library_dirs = ['MOM6_ALE/build_ale'],
                libraries = ['ale'],                        
                sources = ['MOM6_ALE/pyale.f90'])


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "midas",
          version = '1.1',
          description = doclines[0],
          long_description = "\n".join(doclines[2:]),
          author = "Matthew Harrison",
          author_email = "matthew.harrison@noaa.gov",
          url = "none",
          license = 'CCL',
          platforms = ["any"],
          packages=['midas'],
          ext_modules = [hinterp,remap_sfc_fluxes,vertmap_GOLD,vertmap_ALE],
          )
    
