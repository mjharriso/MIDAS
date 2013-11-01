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



from numpy.distutils.core import Extension

remap_sfc_fluxes = Extension(name = 'midas.remap_sfc_fluxes',
                sources = ['remap_sfc_fluxes/remap_sfc_fluxes.f90'])

hinterp = Extension(name = 'midas.hinterp',
                include_dirs = ['fms_build'],
                library_dirs = ['fms_build'],
                libraries = ['fms','netcdf','netcdff'],
                sources = ['hinterp/hinterp.f90'])



vertmap = Extension(name = 'midas.vertmap',
                sources = ['vertmap/midas_vertmap.F90']
                )


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(name = "midas",
          version = '1.0',
          description = "ModularIsosurfaceDataAnalysisSoftware",
          long_description = "specializes in handling gridded datasets for geophysical applications",
          author = "Matthew Harrison",
          author_email = "matthew.harrison@noaa.gov",
          url = "none",
          license = 'CCL',
          platforms = ["any"],
          packages=['midas'],
          ext_modules = [hinterp,remap_sfc_fluxes,vertmap],
          )
    
