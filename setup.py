"""
==============================

 MIDAS specializes in handling data on a 2-dimensional
 lattice of grid cells lying in a 
 cartesian space or on a sphere.  The vertical
 coordinate accomodates generalized surfaces
 i.e. surfaces which are static dynamic in time.

 MIDAS (Modular Isosurface Data Analysis System)
 was first developed by Matthew Harrison
 in 2011-2012 as an employee of NOAA in the 
 GFDL Oceans and Climate Group. MIDAS is being made readily
 available through Github.  #NOAA#GFDL
 

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

===============================
"""

doclines = __doc__.split("\n")


from numpy.distutils.core import Extension

remap_sfc_fluxes = Extension(name = 'midas.remap_sfc_fluxes',
                sources = ['remap_sfc_fluxes/remap_sfc_fluxes.f90'])

hinterp = Extension(name = 'midas.hinterp',
                include_dirs = ['fms_build'],
                library_dirs = ['fms_build'],
                libraries = ['fms','netcdf','netcdff'],
                extra_objects = ['fms_build/*.o'],
                sources = ['hinterp/hinterp.f90'])



vertmap_GOLD = Extension(name = 'midas.vertmap_GOLD',
                sources = ['vertmap_GOLD/vertmap_GOLD.F90']
                )

vertmap_ALE = Extension(name = 'midas.vertmap_ALE',
                include_dirs = ['MOM6_ALE/build_ale'],
                library_dirs = ['MOM6_ALE/build_ale'],
                libraries = ['ale'],                        
                sources = ['MOM6_ALE/pyale.f90'])


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "midas",
          version = '1.0',
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
    
