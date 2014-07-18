"""
 MIDAS (Modular Isosurface Data Analysis System) was first developed by 
 Matthew Harrison 2011-2012 as an employee of NOAA in the 
 GFDL Oceans and Climate Group. The goal of this project is to produce a
 class for handling finite volume representations of numerical model output
 produced by GOLD and later MOM6. The resulting class instances can be spatially
 and temporally processed and saved for additional analysis and visualization.


 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

===============================

 These are REQUIRED packages. If any of these are not present, you
 need to install into your current version of Python. (* == Optional but recommended )

 URLS:
 http://code.google.com/p/netcdf4-python/
 *http://www.scipy.org/
"""

import rectgrid
import rectgrid_gen
import profiles
import wright_eos

if __name__ == "__main__":
    import doctest
    doctest.testmod()
