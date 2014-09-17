"""
 MIDAS (Modular Isolayer Data Analysis System)

 matthew.harrison@noaa.gov (2011-)

 MIDAS is primarily designed to read self-documenting
 files per the CF convention (http://cfconventions.org/)
 describing scalar and vector fields on quadrilateral grid
 meshes in generalized vertical coordinates.

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


import profiles 
import wright_eos
import utils
import rectgrid
import rectgrid_gen
import rectgrid_utils



# end optional packages

if __name__ == "__main__":
    import doctest
    doctest.testmod()
