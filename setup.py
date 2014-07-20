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

doclines = __doc__.split("\n")


if __name__ == '__main__':
    from distutils.core import setup
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
          )
    
