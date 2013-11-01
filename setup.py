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
    
