"""
==============================

 This Python module contains 
 classes for handling data on a 2-dimensional
 lattice of grid cells lying in a 
 cartesian space or on a sphere.  The vertical
 coordinate accomodates Generalized or Lagrangian surfaces
 i.e. surfaces which are dynamic, s(x,y,k,t)
 or static in time, s(x,y,k) or s(x,y).

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

 These are REQUIRED packages. If any of these are not present, you
 need to install into your current version of Python. (* == Optional but recommended )

 URLS:
 http://code.google.com/p/netcdf4-python/
 *http://www.scipy.org/
 http://matplotlib.org/
 *https://pypi.python.org/pypi/seawater/
"""


import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
from  datetime import *
import pickle
import matplotlib.pyplot as plt

# End REQUIRED packages

# MIDAS package

from midas_utils import *
from midas_grid_gen import *
from wright_eos import *
from midas_profiles import *

# Optional packages

try:
    import scipy as sp
    HAVE_SCIPY=True
except:
    HAVE_SCIPY=False
    pass

try:
    import matplotlib.animation as animation
    HAVE_MATPLOTLIB_ANIMATION=True
except:
    HAVE_MATPLOTLIB_ANIMATION=False
    pass

try:
    from seawater import gibbs
except:
    pass

# end optional packages


# Constants 

epsln = 1.e-20
Omega=7.295e-5

# Turn on for more verbose output

DEBUG = 0
   
class rectgrid(object):

  """A rectgrid object describes a horizontal lattice0 on
  a sphere or plane. When instantiating from a file, a rectgrid 
  class is constructed by reading the lat-lon list associated 
  with (var) from file (path) which are placed
  at the centers of the cells at positions x_T and y_T
  Additionally, a grid can be constructed using a supergrid 
  object (see midas_grid_gen.py).

  Additional lattice points are needed to define the 
  cell perimeters. The primary lattice of (T) cell centers, 
  and the perimeter lattice of (Q) points are located at coordinates 
  (y_T,x_T) and (y_T_bounds,x_T_bounds) respectively.  

  There are (jm,im) H cell points and (jm+1,im+1) perimeter
  locations. 
                X_T_bounds
                   V    
           +-------+-------+-------+
           :       :       :       :
           :       :       :       :
           :       :       :       :           
           +-------+-------Q-------:
           :       :       :       :
           :       :   T   :       :< Y_T
           :       :       :       :           
           +-------Q-------+-------+< Y_T_bounds
                       ^   ^
                      X_T

  Cell metrics defining the cell widths and areas are defined by
  dxh,dyh and Ah respectively. Cell orientation is defined by angle_dx.

  """
  
  def __init__(self,path=None,cyclic=False,tripolar_n=False,var=None,simple_grid=False,supergrid=None,lon=None,lat=None,is_latlon=True,is_cartesian=False,grid_type='generic',):

    # """
    
    # >>> grid=rectgrid('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',var='t_an',cyclic=True)
    # >>> print grid.lonh[0],grid.lonh[-1]
    # 0.5 359.5
    # >>> print grid.lath[0],grid.lath[-1]
    # -89.5 89.5
    # """

    self.is_latlon=is_latlon
    self.is_cartesian=is_cartesian
    self.have_metrics = False
    self.simple_grid=simple_grid

    if np.logical_and(is_cartesian,is_latlon):
        print 'Select either is_latlon or is_cartesian, not both'
        return None

    self.yDir=1

    self.cyclic_x = cyclic
    self.tripolar_n = tripolar_n
    
    if supergrid is not None:
        var_dict = {}
        x=supergrid.x; y=supergrid.y
        grid_x=supergrid.grid_x; grid_y=supergrid.grid_y
        self.lonh=grid_x[1::2]
        self.lath=grid_y[1::2]
        self.lonq=grid_x[0::2]
        self.latq=grid_y[0::2]
        self.x_T=x[1::2,1::2]
        self.y_T=y[1::2,1::2]
        self.x_T_bounds=x[0::2,0::2]
        self.y_T_bounds=y[0::2,0::2]

        self.jm = self.x_T.shape[0]
        self.im = self.x_T.shape[1]        

        self.wet = np.ones((self.jm,self.im))
        self.cyclic_x=supergrid.dict['cyclic_x']

        if supergrid.have_metrics:
            dx=supergrid.dx
            dy=supergrid.dy
            dxh = dx[:,::2]+dx[:,1::2]
            self.dxh=0.5*(dxh[:-1:2,:]+dxh[2::2,:])
            dyh = dy[::2,:]+dy[1::2,:]
            self.dyh=0.5*(dyh[:,:-1:2]+dyh[:,2::2])
            self.Ah = self.dxh*self.dyh
            self.angle_dx=supergrid.angle_dx[1::2,1::2]
            self.have_metrics=True
        else:
            self.have_metrics=False


        return

    if path is not None:
        if type(path) == 'netCDF4.Dataset':
            f=path
        else:
            f=nc.Dataset(path)


            
    if grid_type == 'gold_geometry':
        self.x_T = f.variables['geolon'][:]
        self.x_T_bounds = f.variables['geolonb'][:]
        xtb0=2.0*self.x_T_bounds[:,0]-self.x_T_bounds[:,1]
        xtb0=xtb0[:,np.newaxis]
        self.x_T_bounds = np.hstack((xtb0,self.x_T_bounds))
        xtb0=self.x_T_bounds[0,:]
        self.x_T_bounds = np.vstack((xtb0,self.x_T_bounds))    
        self.y_T = f.variables['geolat'][:]
        self.y_T_bounds = f.variables['geolatb'][:]
        ytb0=2.0*self.y_T_bounds[0,:]-self.y_T_bounds[1,:]
        self.y_T_bounds = np.vstack((ytb0,self.y_T_bounds))
        ytb0=self.y_T_bounds[:,0]
        ytb0=ytb0[:,np.newaxis]
        self.y_T_bounds = np.hstack((ytb0,self.y_T_bounds))
        self.lonh = f.variables['lonh'][:] 
        self.lath = f.variables['lath'][:]
        self.lonq = f.variables['lonq'][:]
        self.lonq = np.hstack((self.lonq[0]-(self.lonq[1]-self.lonq[0]),self.lonq))
        self.latq = f.variables['latq'][:]
        self.latq = np.hstack((self.latq[0]-(self.latq[1]-self.latq[0]),self.latq))
        self.D    = f.variables['D'][:]
        self.f    = f.variables['f'][:]

        try:
            self.dxh    = f.variables['dxh'][:]
        except:
            self.dxh    = f.variables['dxT'][:]

        try:
            self.dyh    = f.variables['dyh'][:]
        except:
            self.dyh    = f.variables['dyT'][:]

            
        self.Ah    = f.variables['Ah'][:]

        self.wet    = f.variables['wet'][:]

        self.im = np.shape(self.lonh)[0]
        self.jm = np.shape(self.lath)[0]

        self.have_metrics=True

        return

    if grid_type == 'mom4_gridspec':
        self.x_T = f.variables['geolon_t'][:]
        self.x_T_bounds = f.variables['geolon_e'][:]
        xtb0=2.0*self.x_T_bounds[:,0]-self.x_T_bounds[:,1]
        xtb0=xtb0[:,np.newaxis]
        self.x_T_bounds = np.hstack((xtb0,self.x_T_bounds))
        xtb0=self.x_T_bounds[0,:]
        self.x_T_bounds = np.vstack((xtb0,self.x_T_bounds))    
        self.y_T = f.variables['geolat_t'][:]
        self.y_T_bounds = f.variables['geolat_n'][:]
        ytb0=2.0*self.y_T_bounds[0,:]-self.y_T_bounds[1,:]
        self.y_T_bounds = np.vstack((ytb0,self.y_T_bounds))
        ytb0=self.y_T_bounds[:,0]
        ytb0=ytb0[:,np.newaxis]
        self.y_T_bounds = np.hstack((ytb0,self.y_T_bounds))
        self.lonh = f.variables['gridlon_t'][:] 
        self.lath = f.variables['gridlat_t'][:]
        self.lonq = f.variables['gridlon_c'][:]
        self.lonq = np.hstack((self.lonq[0]-(self.lonq[1]-self.lonq[0]),self.lonq))
        self.latq = f.variables['gridlat_c'][:]
        self.latq = np.hstack((self.latq[0]-(self.latq[1]-self.latq[0]),self.latq))

        self.D    = f.variables['ht'][:]
        self.f    = 2.0*Omega*np.sin(self.y_T*np.pi/180.)

        self.dxh    = f.variables['dxt'][:]
        self.dyh    = f.variables['dyt'][:]

        self.Ah    = self.dxh*self.dyh

        self.wet    = f.variables['wet'][:]


        self.im = np.shape(self.lonh)[0]
        self.jm = np.shape(self.lath)[0]

        self.have_metrics=True

        return

    f = None
    if lon is not None and lat is not None:
        self.lonh=sq(lon[0,:])
        self.lath=sq(lat[:,0])
        self.x_T=lon.copy()
        self.y_T=lat.copy()
    if lon is None and lat is None:
        f=nc.Dataset(path)
  
    
    if var is None and f is not None:
      print """ Need to specify a variable from which to create a
                dummy grid since a valid grid option was not
                specified """
      raise
    else:

      var_dict = {}
      var_dict['X']=None
      var_dict['Y']=None
      var_dict['Z']=None
      var_dict['T']=None

      var_dict['type'] = 'T'

      if f is not None:
          for n in range(0,f.variables[var].ndim):
              dimnam = f.variables[var].dimensions[n]
              dim = f.variables[dimnam]
              cart = get_axis_cart(dim,dimnam)
              if cart is not None:
                  var_dict[cart]=dimnam

          if var_dict['X'] is not None and lon is None:
              lon_axis = f.variables[var_dict['X']]
              dir=get_axis_direction(lon_axis)
              self.lonh = sq(f.variables[var_dict['X']][:])
              if dir == -1:
                  self.lonh=self.lonh[::-1]
              
          elif lon is not None:
              self.lonh = sq(lon[0,:])
          else:
              print "Longitude axis not detected "
              raise
          
          if var_dict['Y'] is not None and lat is None:
              lat_axis = f.variables[var_dict['Y']]
              dir=get_axis_direction(lat_axis)
              self.lath = sq(f.variables[var_dict['Y']][:])
              if dir == -1:
                  self.yDir=-1
                  self.lath=self.lath[::-1]
          elif lat is not None:
              self.lath=sq(lat[:,0])
          else:
              print "Latitude axis not detected "
              raise
      


      self.lonq = 0.5*(self.lonh + np.roll(self.lonh,-1))

      if np.isscalar(self.lonq):
          self.lonq=np.array([self.lonq])

      if np.size(self.lonq) > 2:
          self.lonq[-1] = 2.0*self.lonq[-2] -self.lonq[-3]

      if np.size(self.lonq) > 1:
          lon0 = self.lonq[0]-2.0*(self.lonq[0]-self.lonh[0])
      else:
          lon0 = self.lonq[0]
      
      self.lonq=np.hstack(([lon0],self.lonq))
      self.latq = 0.5*(self.lath + np.roll(self.lath,-1))

      if np.isscalar(self.latq):
          self.latq=np.array([self.latq])

      if np.size(self.latq) > 2:
          self.latq[-1] = 2.0*self.latq[-2]-self.latq[-3]

      if np.size(self.latq) > 1:
          lat_end=self.latq[-1]+2.0*(self.lath[-1]-self.latq[-1])
          self.latq=np.concatenate((self.latq,[lat_end]))      


      try:
          self.im = len(self.lonh)
      except:
          self.lonh = np.array([self.lonh])
          self.im = 1

      try:
          self.jm = len(self.lath)
      except:
          self.lath = np.array([self.lath])          
          self.jm = 1
      
      if simple_grid is True:
          self.simple_grid = True
          self.cyclic_x = cyclic          
          return


      if np.logical_and(self.im > 1, self.jm > 1):
          self.x_T,self.y_T = np.meshgrid(self.lonh,self.lath)

          xtb=0.5*(self.x_T + np.roll(self.x_T,shift=1,axis=1))
          xtb0=2.0*xtb[:,-1]-xtb[:,-2]
          xtb0=xtb0[:,np.newaxis]
          xtb=np.hstack((xtb,xtb0))
          xtb0=2.0*xtb[:,1]-xtb[:,2]
          xtb[:,0]=xtb0
          self.x_T_bounds=xtb.copy()
          xtb0=self.x_T_bounds[-1,:]
          self.x_T_bounds = np.vstack((self.x_T_bounds,xtb0))    
          
          ytb=0.5*(self.y_T + np.roll(self.y_T,shift=1,axis=0))
          ytb0=2.0*ytb[-1,:]-ytb[-2,:]
          ytb0=ytb0[np.newaxis,:]
          ytb=np.vstack((ytb,ytb0))
          ytb0=2.0*ytb[1,:]-ytb[2,:]
          ytb[0,:]=ytb0
          self.y_T_bounds=ytb.copy()
          ytb0=self.y_T_bounds[:,-1]
          ytb0=ytb0[:,np.newaxis]
          self.y_T_bounds = np.hstack((self.y_T_bounds,ytb0))          
          
          dx = (self.x_T_bounds - np.roll(self.x_T_bounds,axis=1,shift=1))
          dx=dx[1:,1:]
          dx=dx*np.pi/180.
          self.dxh = dx*R_earth*np.cos(self.y_T*np.pi/180.)
          dy = (self.y_T_bounds - np.roll(self.y_T_bounds,axis=0,shift=1))
          dy = dy[1:,1:]
          dy = dy*np.pi/180.
          self.dyh = dy*R_earth
          self.Ah=self.dxh*self.dyh
          self.have_metrics=True
      else:
          self.x_T=self.lonh; self.y_T=self.lath
          

      self.cyclic_x = cyclic

      if self.cyclic_x:
        self.xmod_len = self.x_T_bounds[0,-1] - self.x_T_bounds[0,0] 

      self.wet = np.ones((self.jm,self.im))


      
  def find_geo_bounds(self,x=None,y=None):
      """
      Returns the bounds of the grid. Currently this is of limited use 
      for generalized horizontal coordinates (based on lonh/lath).
      

      >>> from midas import *
      >>> x=np.arange(0.5,360.5,1.0);y=np.arange(-89.5,90.5,1.0)
      >>> X,Y=np.meshgrid(x,y)
      >>> grid=rectgrid(lon=X,lat=Y,cyclic=True)
      >>> xs,xe,ys,ye=grid.find_geo_bounds(x=(20.,50.),y=(-10.,10.))
      >>> print xs,xe,grid.lonh[xs],grid.lonh[xe]
      20 50 20.5 50.5
      >>> print ys,ye,grid.lath[ys],grid.lath[ye]
      79 99 -10.5 9.5
      """


      xs=0; xe=self.im
      ys=0; ye=self.jm

      if x is not None:
          [xs,xe]=find_axis_bounds(self.lonq,x=x,modulo_360=self.cyclic_x)

      if y is not None:
          [ys,ye]=find_axis_bounds(self.latq,x=y)

      return xs,xe,ys,ye
          
  def geo_region(self,y=None,x=None,name=None):

     """
     Returns a \"section\" dictionary for sampling a contiguous region
     based on geographical boundaries. For a lat-lon region returned, 
     the mesh of grid cells is bounded by the first grid tracer point 
     reached to the right of the left domain boundary and the last grid cell
     to the left of the right boundary.
     
     Currently this is of limited use for generalized horizontal 
     coordinates (based on lonh/lath). [ Instead I typically use 
     the indexed_region method to slice general grids]. 
      

     >>> from midas import *
     >>> x=np.arange(0.5,360.5,1.0);y=np.arange(-89.5,90.5,1.0)
     >>> X,Y=np.meshgrid(x,y)
     >>> grid=rectgrid(lon=X,lat=Y,cyclic=True)
     >>> section=grid.geo_region(x=(-30.,20.),y=(-10.,10.))
     >>> print section['xax_data'][0],section['xax_data'][-1]
     330.5 379.5
     >>> print section['yax_data'][0],section['yax_data'][-1]
     -10.5 8.5

     """

  
     section={}

     xs,xe,ys,ye = self.find_geo_bounds(x=x,y=y)

     if xe == xs:
         xe=xe+2


     if ye == ys:
         ye=ye+2

     section['y']=np.arange(ys,ye)
     section['yax_data']= self.lath[section['y']]


     if xe>=xs:
         section['x']=np.arange(xs,xe)
         section['x_read']=section['x']
         section['xax_data']=self.lonh[section['x']]
     else:
         section['x_read']=[np.arange(xs,self.im)]
         section['x_read'].append(np.arange(0,xe))
         xind = np.hstack((section['x_read'][0],section['x_read'][1]))
         section['x'] = xind
         section['xax_data']=self.lonh[xind]

     lonh=section['xax_data'].copy()

     if not np.isscalar(lonh):
         lon0=lonh[0]
     else:
         lon0=lonh

     if self.cyclic_x:
         lonh[lonh<lon0]=lonh[lonh<lon0]+360.

     section['xax_data']=lonh

     section['lon0']=lon0


     section['name'] = name
     section['parent_grid'] = self

     return section

  def indexed_region(self,i=None,j=None,name=None):

    """
    Returns a \"section\" dictionary for sampling a contiguous region. 
    based on index coordinates.
      

    >>> from midas import *
    >>> x=np.arange(0.5,360.5,1.0);y=np.arange(-89.5,90.5,1.0)
    >>> X,Y=np.meshgrid(x,y)
    >>> grid=rectgrid(lon=X,lat=Y,cyclic=True)
    >>> section=grid.indexed_region(i=(20,20))
    >>> print section['xax_data'][0],section['xax_data'][-1]
    20.5 20.5
    >>> print section['yax_data'][0],section['yax_data'][-1]
    -89.5 89.5
    """
  
    section={}

    if j is not None:
      ys=j[0];ye=j[1]
      section['y']=np.arange(ys,ye+1)
      section['yax_data']= self.lath[section['y']]
    else:
      section['y']=None
      section['yax_data']= self.lath
      
    if i is not None:
        xs=i[0];xe=i[1]
        if xe>=xs:
            section['x']=np.arange(xs,xe+1)
        else:
            section['x']=np.hstack((np.arange(xs,self.im),np.arange(0,xe)))
        section['xax_data']= self.lonh[section['x']]
    else:
      section['x']=None

    section['name'] = name
    section['parent_grid'] = self

    return section
      
  def extract(self,geo_region=None):
    """
    Returns new grid object using a \"section\" dictionary 
    created using the \"geo_region\" or \"indexed_region\"
    method.
      

    >>> from midas import *
    >>> import hashlib
    >>> x=np.arange(0.5,360.5,1.0);y=np.arange(-89.5,90.5,1.0)
    >>> X,Y=np.meshgrid(x,y)
    >>> grid=rectgrid(lon=X,lat=Y,cyclic=True)
    >>> section=grid.geo_region(x=(-30.,20.),y=(-10.,10.))
    >>> new_grid = grid.extract(section)
    >>> hash=hashlib.md5(new_grid.x_T)
    >>> hash.update(new_grid.y_T)
    >>> print hash.hexdigest()
    465f6aa54e3a672bb4a1d8cc02b7e96c
    """

    if geo_region is None:
        grid=copy.copy(self)
        return grid
    else:
      grid = copy.copy(self)

      x_section = geo_region['x']
      y_section = geo_region['y']      

      if x_section is not None:
          xb_section = np.hstack((x_section,x_section[-1]+1))
      if y_section is not None:
          yb_section = np.hstack((y_section,y_section[-1]+1))        

      grid.lath = np.take(self.lath,y_section,axis=0)
      grid.latq = np.take(self.latq,yb_section,axis=0)

      grid.lonh=np.take(self.lonh,x_section,axis=0)
      grid.lonq=np.take(self.lonq,xb_section,axis=0)

      lon0=grid.lonq[0]

      grid.lonh[grid.lonh<lon0]=grid.lonh[grid.lonh<lon0]+360.
      grid.lonq[grid.lonq<lon0]=grid.lonq[grid.lonq<lon0]+360.

      grid.x_T = np.take(np.take(self.x_T,y_section,axis=0),x_section,axis=1)
#      grid.x_T[grid.x_T<lon0]=grid.x_T[grid.x_T<lon0]+360.
      grid.x_T_bounds = np.take(np.take(self.x_T_bounds,yb_section,axis=0),xb_section,axis=1)
#      grid.x_T_bounds[grid.x_T_bounds<lon0]=grid.x_T_bounds[grid.x_T_bounds<lon0]+360.

      grid.y_T = np.take(np.take(self.y_T,y_section,axis=0),x_section,axis=1)
      grid.y_T_bounds = np.take(np.take(self.y_T_bounds,yb_section,axis=0),xb_section,axis=1)

      if hasattr(grid,'D'):
          grid.D = np.take(np.take(self.D,y_section,axis=0),x_section,axis=1)

      if hasattr(grid,'f'):
          grid.f = np.take(np.take(self.f,y_section,axis=0),x_section,axis=1)

      if hasattr(grid,'wet'):
          grid.wet = np.take(np.take(self.wet,y_section,axis=0),x_section,axis=1)

          
      if grid.have_metrics:
          grid.dxh = np.take(np.take(self.dxh,y_section,axis=0),x_section,axis=1)
          grid.dyh = np.take(np.take(self.dyh,y_section,axis=0),x_section,axis=1)
          grid.Ah = np.take(np.take(self.Ah,y_section,axis=0),x_section,axis=1)

          
      if hasattr(grid,'mask'):
          grid.mask = np.take(np.take(self.mask,y_section,axis=0),x_section,axis=1)

      grid.im = np.shape(grid.lonh)[0]
      grid.jm = np.shape(grid.lath)[0]

      return grid

  def add_mask(self,field,path=None):

    """
       Add a 2-D mask to the grid. The mask can have any values, e.g.
       a different number for each ocean basin. This can be used to
       defne other masked arrays using mask_where, for instance.

    """

    if path is not None:
      f=nc.Dataset(path)
      
    if field in f.variables:
      self.mask = f.variables[field][:]
    else:
      print ' Field ',field,' is not present in file ',path
      return None
    

    return None


class state(object):
  """Returns a model state of (fields). The default is 
  to extract the entire data domain at all time levels from (path). Use
  (geo_region) to store a section of the horizontal grid. Use
  (time_indices) or (date_bounds) to extract along the record dimension.
  (z_indices) can be used to read a contiguous number of vertical layers or
  levels.

  stagger:
     11 - centered at grid tracer (T) points, (lath,lonh)
     21 - centered on north face of tracer cell
     12 - centered on east face of tracer cell
     22 - centered on north-east corner of tracer cell
     01 - centered on south face of tracer cell
     10 - centered on west face of tracer cell
     00 - centered on south-west corner of tracer cell          
     
     
  NOTE:The layer (interfaces) variable must be stored in order to calculate
  accurate finite volume integrals.

  """
  
  def __init__(self,path=None,grid=None,geo_region=None,time_indices=None,date_bounds=None,z_indices=None,fields=None,default_calendar=None,MFpath=None,interfaces=None,path_interfaces=None,MFpath_interfaces=None,stagger=None,verbose=True):
    """
#    >>> from midas import *
#    >>> import hashlib
#    >>> grid=rectgrid('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',var='t_an',cyclic=True)
#    >>> IO = grid.geo_region(x=(30,120.),y=(-30.,25.))
#    >>> S=state('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',grid=grid,geo_region=IO,fields=['t_an'])
#    >>> hash=hashlib.md5(S.t_an)
#    >>> print hash.hexdigest()
  
    """

    if path is not None:
      f=nc.Dataset(path)
      self.path = path
      self.is_MFpath = False
    elif MFpath is not None:
      f=nc.MFDataset(MFpath)
      self.path = MFpath
      self.is_MFpath = True      
    else:
      f=None

    self.rootgrp = f


    self.variables = {}
    self.var_dict = {}    

    z_indices_in = z_indices

    if fields is None:
      if grid is not None:
        new_grid = grid.extract(geo_region)
        self.grid = new_grid
      return

    if stagger is  None:
        stagger = {}
        for v in fields:
            stagger[v]='00'
    else:
        var_stagger = stagger.copy()
        if len(var_stagger) != len(fields):
            print """ Need to provide stagger for each field """
            return
        stagger = {}
        n=0
        for v in fields:
            stagger[v]=var_stagger[n]
            n=n+1

    if grid is None:
        if path is not None:
            grid = rectgrid(path,var=fields[0])
        else:
            grid = rectgrid(self.rootgrp,var=fields[0])

        if geo_region is not None:
            new_grid=grid.extract(geo_region)
            self.grid = new_grid
        else:
            self.grid = grid
    else:
        new_grid = grid.extract(geo_region)
        self.grid = new_grid
        
    for v in fields:
       try: 
         self.variables[v] = self.rootgrp.variables[v]  # netCDF4 variable object
       except KeyError:
         print """ Variable named %(v)s does not exist in file named %(path)s
         Aborting ... """%{'v':v,'path':self.path}
         raise
          
       var_dict = {}


       var_dict['X'] = None; var_dict['Y'] = None
       var_dict['Z'] = None; var_dict['T'] = None

       var_dict['path']=self.path

       var_dict['stagger']=stagger[v]
       
       var_dict['is_MFpath']=self.is_MFpath

       for n in range(0,self.variables[v].ndim):

# For each dimension , determine its Cartesian attribute
# This is a critical step! Things may appear to work successfully
# but fail at a later stage if the X/Y/Z/T orientation of the
# data stored is not associated.


         dimnam = self.variables[v].dimensions[n]
         dim = self.rootgrp.variables[dimnam]
         cart = get_axis_cart(dim,dimnam)

         if cart is not None:
           var_dict[cart]=self.variables[v].dimensions[n]
           
         if cart == 'Z':
           var_dict['Zdir'] = get_axis_direction(dim)
           try:
               zunits = getattr(dim,'units')
           except:
               zunits='none'
           

           if interfaces is not None:
               var_dict['Ztype'] = 'Generalized'
           else:
               var_dict['Ztype'] = 'Fixed'
               
           try:
             var_dict['Zb'] = getattr(self.rootgrp.variables[var_dict['Z']],'bounds')
           except:
               try:
                   var_dict['Zb'] = getattr(self.rootgrp.variables[var_dict['Z']],'edges')
               except:
                   var_dict['Zb'] = None


             
         if cart == 'T':
           nt = len(self.rootgrp.variables[var_dict['T']][:])
           t_indices = np.arange(0,nt)
           var_dict['tax_data'] = self.rootgrp.variables[var_dict['T']][:]
           try:
               var_dict['tunits'] = self.rootgrp.variables[var_dict['T']].units
           except:
               var_dict['tunits'] = 'none'
               
           try:
             var_dict['calendar'] = string.lower(self.rootgrp.variables[var_dict['T']].calendar)
           except:
             if default_calendar is not None:
               var_dict['calendar']=default_calendar
             else:
               var_dict['calendar']=None

           var_dict['tunits']=var_dict['tunits'].replace('0000','0001')

           
           if var_dict['calendar']=='365_days':
               var_dict['calendar']='365_day'

           if var_dict['tunits'].count('month') > 0:
               mon_to_day=365.0/12.0
               var_dict['tunits'] = 'days'+var_dict['tunits'][6:]
               var_dict['tax_data']=var_dict['tax_data'][:]*mon_to_day
               
           if var_dict['calendar'] is not None:
             var_dict['dates'] = num2date(var_dict['tax_data'],var_dict['tunits'],var_dict['calendar'])

           self.calendar = var_dict['calendar']
           
       if var_dict['T'] is not None:
         if time_indices is not None and date_bounds is not None:
           print """
              Invalid options. select either time_indices or date_bounds
              not both. """
           raise()
         try:
           Tb = getattr(self.rootgrp.variables[var_dict['T']],'bounds')
         except:
           Tb = None
           
         if Tb is not None:
           if self.rootgrp.variables[Tb].ndim == 2:
             var_dict['tbax_data'] = self.rootgrp.variables[Tb][:,0]
             
             try:
               tb_last = self.rootgrp.variables[Tb][-1,1]
             except:
               tb_last = self.rootgrp.variables[Tb][0,1]

             var_dict['tbax_data'] = np.hstack((var_dict['tbax_data'],[tb_last]))
           else:
             var_dict['tbax_data'] = self.rootgrp.variables[Tb][:]

         else:
           tdat=var_dict['tax_data']
           if len(tdat) > 1:
             tint=np.hstack((1.5*tdat[0]-0.5*tdat[1],0.5*(tdat[0:-1]+tdat[1:])))
             tint=np.hstack((tint,tint[-1]+tdat[-1]-tdat[-2]))
             var_dict['tbax_data']=tint
           else:
             var_dict['tbax_data']=None

         if var_dict['tbax_data'] is not None and var_dict['calendar'] is not None:
           
           var_dict['date_bounds'] = num2date(var_dict['tbax_data'],var_dict['tunits'],var_dict['calendar'])
           
         if time_indices is not None:
           t_indices = time_indices
           tb_indices = np.hstack((time_indices,time_indices[-1]+1))
           var_dict['tax_data']=var_dict['tax_data'][t_indices]
           if var_dict['tbax_data'] is not None:
             var_dict['tbax_data']=var_dict['tbax_data'][tb_indices]
           if var_dict['calendar'] is not None:
             var_dict['dates']=var_dict['dates'][t_indices]
             if var_dict['tbax_data'] is not None:
               var_dict['date_bounds']=var_dict['date_bounds'][tb_indices]
         elif date_bounds is not None:
           if var_dict['calendar'] is not None:
             ts,te = find_date_bounds(var_dict['date_bounds'][:],date_bounds[0],date_bounds[1])
             t_indices=np.arange(ts,te)
             tb_indices=np.arange(ts,te+1)             
             var_dict['tax_data']=var_dict['tax_data'][t_indices]
             var_dict['tbax_data']=var_dict['tbax_data'][tb_indices]                        
             var_dict['dates']=var_dict['dates'][t_indices]
             var_dict['date_bounds']=var_dict['date_bounds'][tb_indices]             
           else:
             print """
             Calendar is inactive for %(field)s """%{'field':v}
             raise()
       else:
         t_indices = None
         


       if var_dict['Z'] is not None:
           try:
               var_dict['zunits'] =   self.rootgrp.variables[var_dict['Z']].units
           except:
               var_dict['zunits'] = 'none'
               
           if z_indices is not None:
               tmp = np.array([z_indices[-1]+1])
               z_interfaces = np.hstack((z_indices,tmp))
               nz = len(z_indices)
           else:
               nz = len(self.rootgrp.variables[var_dict['Z']][:])
               z_indices=np.arange(0,nz)
               if len(z_indices) == 0:
                   z_indices = 0
               z_interfaces = np.arange(0,nz+1)

           var_dict['zax_data']= self.rootgrp.variables[var_dict['Z']][z_indices]

           zdat=f.variables[var_dict['Z']][z_indices]
           if len(zdat) > 1:
               zint=np.hstack((1.5*zdat[0]-0.5*zdat[1],0.5*(zdat[0:-1]+zdat[1:])))
               zint=np.hstack((zint,zint[-1]+zdat[-1]-zdat[-2]))
               var_dict['zbax_data']=zint
           else:
               var_dict['zbax_data']=None
           
           if var_dict['Zb'] is not None:
                 
               var_dict['zbax_data'] = self.rootgrp.variables[var_dict['Zb']][:]

               if var_dict['Ztype'] is 'Fixed':
                   if self.rootgrp.variables[var_dict['Zb']].ndim == 2:
                       zb_last = self.rootgrp.variables[var_dict['Zb']][-1,1]
                       var_dict['zbax_data']=var_dict['zbax_data'][:,0]
                       var_dict['zbax_data'] = np.hstack((var_dict['zbax_data'],[zb_last]))                     
               else:
                   zdat=var_dict['zax_data']
                   if len(zdat) > 1:
                       zint=np.hstack((1.5*zdat[0]-0.5*zdat[1],0.5*(zdat[0:-1]+zdat[1:])))
                       zint=np.hstack((zint,zint[-1]+zdat[-1]-zdat[-2]))
                       var_dict['zbax_data']=zint
                   else:
                       var_dict['zbax_data']=None

           if var_dict['Z'] == 'zi':   ### Using the "zi"dimension name to detect if this is an interface variable. Only works for GOLD output
               var_dict['interface_variable']=True
           else:
               var_dict['interface_variable']=False
           
       else:
         z_indices = None
         z_interfaces = None
         var_dict['interface_variable']=False

       var_dict['z_indices']=z_indices
       var_dict['zb_indices']=z_interfaces

       x_indices = np.arange(0,self.grid.im)
       x_indices_read = np.arange(0,self.grid.im)
       y_indices = np.arange(0,self.grid.jm)
       
       if var_dict['Y'] is not None and var_dict['X'] is not None:
           try:
               var_dict['xunits'] =   self.rootgrp.variables[var_dict['X']].units
           except:
               var_dict['xunits'] =  'none'
           try:
               var_dict['yunits'] =   self.rootgrp.variables[var_dict['Y']].units
           except:
               var_dict['yunits'] =  'none'
           if geo_region is not None:
               var_dict['yax_data'] = geo_region['yax_data'][:]
               y_indices=geo_region['y'][:]
               x_indices=geo_region['x'][:]
               x_indices_read = geo_region['x_read'][:]
               var_dict['xax_data'] = geo_region['xax_data'][:]
               var_dict['yax_data'] = geo_region['yax_data'][:]
               if grid.yDir == -1:
                   var_dict['yax_data']=var_dict['yax_data'][::-1]
           else:
               var_dict['yax_data'] = grid.lath
               var_dict['xax_data'] = grid.lonh
               if grid.yDir == -1:
                   var_dict['yax_data']=var_dict['yax_data'][::-1]

       if var_dict['X'] is None:
         x_indices = None


       if var_dict['Y'] is None:
         y_indices = None         


       slice_read = [];shape_read = []; shape_out=[]; slice_out = []


       for s in [t_indices,z_indices,y_indices]:
           if s is not None:
               if len(s) == 1:  # This seems necessary due to a NC4 bug
                   slice_read.append(s[0])
               else:
                   slice_read.append(s)
               slice_out.append(np.arange(0,s.shape[0]))             
               shape_read.append(s.shape[0])
               shape_out.append(s.shape[0])           
           else:
               slice_out.append([0])
               shape_out.append(1)
               shape_read.append(1)

       if geo_region is not None or x_indices is not None:
           slice_read.append(x_indices)
           slice_out.append(np.arange(0,x_indices.shape[0]))
           shape_read.append(x_indices.shape[0])
           shape_out.append(x_indices.shape[0])
       else:
           slice_out.append([0])
           shape_out.append(1)
           shape_read.append(1)

       slice_int_read = [];shape_int_read = []; shape_int_out=[]; slice_int_out = []
       for s in [t_indices,z_interfaces,y_indices]:
         if s is not None:
             if len(s) == 1:  # This seems necessary due to a NC4 bug
                 slice_int_read.append(s[0])
             else:
                 slice_int_read.append(s)
             slice_int_out.append(np.arange(0,s.shape[0]))                      
             shape_int_read.append(s.shape[0])
             shape_int_out.append(s.shape[0])           
         else:
             slice_int_out.append([0])
             shape_int_out.append(1)
             shape_int_read.append(1)

       if geo_region is not None or x_indices is not None:
           slice_int_read.append(x_indices)
           slice_int_out.append(np.arange(0,x_indices.shape[0]))
           shape_int_read.append(x_indices.shape[0])
           shape_int_out.append(x_indices.shape[0])
       else:
         slice_int_out.append([0])
         shape_int_out.append(1)
         shape_int_read.append(1)

       if var_dict['interface_variable'] == False:
         self.slice_read = slice_read
         self.slice_out = slice_out
         self.shape_read = shape_read
         self.shape_out = shape_out

         self.slice_int_read = slice_int_read
         self.slice_int_out = slice_int_out
         self.shape_int_read = shape_int_read
         self.shape_int_out = shape_int_out       


       data_read = np.array(self.rootgrp.variables[v][slice_read])
       data_read = np.reshape(data_read,(shape_read))

       vars(self)[v] = data_read

       if grid is not None:
           if self.grid.yDir == -1:
               vars(self)[v] = vars(self)[v][:,:,::-1,:]
           
       if DEBUG == 1 or verbose == True:
         print " Successfully extracted data named %(nam)s from %(fil)s "%{'nam':v,'fil':self.path}
         print " Resulting shape = ",vars(self)[v].shape
         print " Dictionary keys = ", var_dict.keys()

                  


       var_dict['masked']=False
       var_dict['_FillValue'] = None
       var_dict['missing_value'] = None

         
       for i in f.variables[v].ncattrs():
         if i == 'units':
           var_dict['units']=f.variables[v].units
         if i == '_FillValue':
           var_dict['_FillValue'] = f.variables[v]._FillValue
         if i == 'missing_value':
           var_dict['missing_value'] = f.variables[v].missing_value

       if var_dict['_FillValue'] is not None or var_dict['missing_value']  is not None:
         if var_dict['_FillValue'] is not None:
             vars(self)[v] = np.ma.masked_where(np.abs(vars(self)[v] - var_dict['_FillValue']) < np.abs(var_dict['_FillValue']*.01),vars(self)[v])
             var_dict['missing_value']=var_dict['_FillValue']
         else:
             vars(self)[v] = np.ma.masked_where(np.abs(vars(self)[v] - var_dict['missing_value']) < np.abs(var_dict['missing_value']*0.01),vars(self)[v])
             var_dict['_FillValue']=var_dict['missing_value']

         var_dict['masked']=True
         
       var_dict['z_interfaces'] = None


       if interfaces is not None:

         interfaces_exist=False
         
         try:
           vars(self)[interfaces]
           interfaces_exist=True
         except:
           pass

         if not interfaces_exist:

             if type(interfaces) == str:
                 if path_interfaces is not None:
                     f_interfaces = nc.Dataset(path_interfaces)
                     data_int_read = np.reshape(np.array(f_interfaces.variables[interfaces][slice_int_read]),(shape_int_read))
                 elif MFpath_interfaces is not None:
                     f_interfaces = nc.MFDataset(MFpath_interfaces)
                     data_int_read = np.reshape(np.array(f_interfaces.variables[interfaces][slice_int_read]),(shape_int_read))
                     path_interfaces = MFpath_interfaces
                 else:
                     data_int_read = np.reshape(np.array(self.rootgrp.variables[interfaces][slice_int_read]),(shape_int_read))
                     path_interfaces = self.path
           
                 vars(self)[interfaces] = data_int_read

                 if DEBUG == 1 or verbose == True:
                     print " Successfully extracted interface data named %(nam)s from %(fil)s "%{'nam':interfaces,'fil':path_interfaces}
                     print " Resulting shape = ",vars(self)[interfaces].shape
                     
             else:
                 zint=np.take(np.take(np.take(np.take(interfaces,slice_int_read[3],axis=3),slice_int_read[2],axis=2),slice_int_read[1],axis=1),slice_int_read[0],axis=0)
                 data_int_read = np.reshape(zint,(shape_int_read))

                 interfaces_name='e'
                 vars(self)[interfaces_name] = data_int_read

             self.interfaces=interfaces
           



       else:
           self.interfaces = None
           
       if var_dict['Z'] is not None and var_dict['Ztype'] is not 'Generalized':
         if var_dict['zbax_data'] is not None:
           tmp = np.reshape(var_dict['zbax_data'][z_interfaces],(nz+1,1,1))
           if geo_region is not None:
             ny = len(geo_region['y'])
             nx = len(geo_region['x'])
           else:
             ny = len(self.rootgrp.variables[var_dict['Y']][:])
             nx = len(self.rootgrp.variables[var_dict['X']][:])
              
           var_dict['z_interfaces']  = var_dict['Zdir']*np.tile(tmp,(1,ny,nx))
           tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))             
           var_dict['z'] = tmp[0:-1,:,:]
           var_dict['dz'] = var_dict['z_interfaces'] - np.roll(var_dict['z_interfaces'],axis=0,shift=-1)
           var_dict['dz'] = var_dict['dz'][0:-1,:,:]
             
       if var_dict['Z'] is not None and var_dict['Ztype'] is 'Generalized' and interfaces is not None:
           if type(interfaces) == str:
               var_dict['z_interfaces']  = vars(self)[interfaces]
           else:
               var_dict['z_interfaces']  = interfaces
               
           tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=1,shift=-1))             
           var_dict['z'] = tmp[:,0:-1,:,:]
           var_dict['dz'] = var_dict['z_interfaces'] - np.roll(var_dict['z_interfaces'],axis=1,shift=-1)
           var_dict['dz'] = var_dict['dz'][:,0:-1,:,:]
         
       try:
         time_avg_info = getattr(self.rootgrp.variables[v],'time_avg_info')
       except:
         time_avg_info = ""
         
       if 'average_DT' in time_avg_info:
           try:
               dt = self.rootgrp.variables['average_DT'][t_indices]
               units = self.rootgrp.variables['average_DT'].units
           except:
               dt = np.ones((len(t_indices)))
               units = None

           if units == 'days':
               dt = dt*86400. # convert to seconds
           elif units == 'hours':
               dt = dt*3600. # convert to seconds
       else:
         if var_dict['T'] is not None:           
           dt = np.ones((len(t_indices)))
         else:
           dt = 0.0



       var_dict['dt'] = dt

       var_dict['rootgrp'] = self.rootgrp
          
       self.var_dict[v] = var_dict
       self.geo_region = geo_region

       del t_indices,y_indices,x_indices
           

    f.close()
    
  def add_field(self,field,path=None,MFpath=None,use_interfaces=False,verbose=True):
    """
    Add a field to the existing state (e.g. more tracers).
    either from the current root file or from an alternate path.
    """

    if path is None and MFpath is None:
      f = self.rootgrp
      path = self.path
    elif path is not None:
      f=nc.Dataset(path)
    elif MFpath is not None:
      f=nc.MFDataset(MFpath)      
      


    if field in f.variables:
      nam = string.join(['self',field],sep='.')
      var_dict = {}
      var_dict['rootgrp'] = f
      self.variables[field] = f.variables[field]  # netCDF4 variable object
      var_dict['X'] = None; var_dict['Y'] = None
      var_dict['Z'] = None; var_dict['T'] = None
    else:
      print ' Field ',field,' is not present in file ',path

    if MFpath is None:
      var_dict['path'] = path
    else:
      var_dict['path'] = MFpath

    t_indices = None
    
    for n in range(0,self.variables[field].ndim):
      dimnam = self.variables[field].dimensions[n]
      dim    = f.variables[dimnam]
      cart   = get_axis_cart(dim,dimnam)
      
      if cart is not None:
        var_dict[cart]=self.variables[field].dimensions[n]

      try:
        var_dict['Zb'] = getattr(f.variables[var_dict['Z']],'bounds')
      except:
        var_dict['Zb'] = None


      if cart == 'T':
        t_indices = self.slice_read[0]
        nt = len(f.variables[var_dict['T']][t_indices])        
        var_dict['tax_data'] = f.variables[var_dict['T']][t_indices]
        var_dict['tunits'] = f.variables[var_dict['T']].units
        var_dict['calendar'] = self.calendar

        if var_dict['tunits'].count('month') > 0:
            mon_to_day=365.0/12.0
            var_dict['tunits'] = 'days'+var_dict['tunits'][6:]
            var_dict['tax_data']=var_dict['tax_data'][:]*mon_to_day
        
        if var_dict['calendar'] is not None:
          var_dict['dates'] = num2date(var_dict['tax_data'],var_dict['tunits'],var_dict['calendar'])
          
        var_dict['t_indices'] = t_indices


    
    if var_dict['T'] is not None:
        try:
            Tb = getattr(f.variables[var_dict['T']],'bounds')
        except:
            Tb = None

        if Tb is not None:
            if f.variables[Tb].ndim == 2:
                var_dict['tbax_data'] = f.variables[Tb][:,0]
                
                try:
                    tb_last = f.variables[Tb][-1,1]
                except:
                    tb_last = f.variables[Tb][0,1]

                var_dict['tbax_data'] = np.hstack((var_dict['tbax_data'],[tb_last]))
            else:
                var_dict['tbax_data'] = f.variables[Tb][:]

            if t_indices is not None:
                tb_indices = np.hstack((t_indices,t_indices[-1]+1))
                if var_dict['tbax_data'] is not None:
                    var_dict['tbax_data']=var_dict['tbax_data'][tb_indices]
                if var_dict['tbax_data'] is not None and var_dict['calendar'] is not None:
                    var_dict['date_bounds']=num2date(var_dict['tbax_data'],var_dict['tunits'],var_dict['calendar'])
                
        else:
            tdat=var_dict['tax_data']
            if len(tdat) > 1:
                tint=np.hstack((1.5*tdat[0]-0.5*tdat[1],0.5*(tdat[0:-1]+tdat[1:])))
                tint=np.hstack((tint,tint[-1]+tdat[-1]-tdat[-2]))
                var_dict['tbax_data']=tint
            else:
                var_dict['tbax_data']=None

        if var_dict['tbax_data'] is not None and var_dict['calendar'] is not None:
            var_dict['date_bounds'] = num2date(var_dict['tbax_data'],var_dict['tunits'],var_dict['calendar'])
           
      
    if var_dict['Z'] is not None:
        z_indices = self.slice_read[1]
        zdat=f.variables[var_dict['Z']][z_indices]        
        nz=len(z_indices)
        z_interfaces = self.slice_int_read[1]
        var_dict['zax_data']= f.variables[var_dict['Z']][z_indices]
        if var_dict['Zb'] is not None:
            var_dict['zbax_data'] = f.variables[var_dict['Zb']][z_interfaces]
        if len(zdat) > 1:
            zint=np.hstack((1.5*zdat[0]-0.5*zdat[1],0.5*(zdat[0:-1]+zdat[1:])))
            zint=np.hstack((zint,zint[-1]+zdat[-1]-zdat[-2]))
            var_dict['zbax_data']=zint[z_interfaces]
        else:
            var_dict['zbax_data']=None

        if var_dict['Z'] == 'zi':
            var_dict['interface_variable']=True
        else:
            var_dict['interface_variable']=False

        if self.interfaces is not None:
            var_dict['Ztype'] = 'Generalized'
        else:
            var_dict['Ztype'] = 'Fixed'
            
    else:
        z_indices = None
        z_interfaces = None
        var_dict['interface_variable']=False
      
    var_dict['z_indices']=z_indices
    var_dict['zb_indices']=z_interfaces

    
    geo_region = self.geo_region
    
          
    if var_dict['interface_variable']:
      data_read = np.reshape(np.array(f.variables[field][self.slice_int_read]),(self.shape_int_read))
    else:
      data_read = np.reshape(np.array(f.variables[field][self.slice_read]),(self.shape_read))

    vars(self)[field] = data_read
    
    if var_dict['Y'] is not None and var_dict['X'] is not None:
      var_dict['xunits'] =   var_dict['rootgrp'].variables[var_dict['X']].units
      var_dict['yunits'] =   var_dict['rootgrp'].variables[var_dict['Y']].units 
      if geo_region is not None:
        var_dict['yax_data'] = geo_region['yax_data'][:]
        var_dict['xax_data'] = geo_region['xax_data'][:]
      else:
        ny = len(var_dict['rootgrp'].variables[var_dict['Y']][:])
        y_indices = np.arange(0,ny)
        var_dict['yax_data'] = var_dict['rootgrp'].variables[var_dict['Y']][y_indices]
        nx = len(var_dict['rootgrp'].variables[var_dict['X']][:])
        x_indices = np.arange(0,nx)
        var_dict['xax_data'] = var_dict['rootgrp'].variables[var_dict['X']][x_indices]

         
    if DEBUG == 1 or verbose == True:
      print " Successfully extracted data named %(nam)s from %(fil)s "%{'nam':field,'fil':var_dict['path']}
      print " Resulting shape = ",vars(self)[field].shape                 


    var_dict['masked']=False
    var_dict['_FillValue'] = None
    var_dict['missing_value'] = None

    
    for i in f.variables[field].ncattrs():
      if i == 'units':
        var_dict['units']=f.variables[field].units      
      if i == '_FillValue':
        var_dict['_FillValue'] = f.variables[field]._FillValue
      if i == 'missing_value':
        var_dict['missing_value'] = f.variables[field].missing_value

    if var_dict['_FillValue'] is not None or var_dict['missing_value']  is not None:
      if var_dict['_FillValue'] is not None:
        vars(self)[field] = np.ma.masked_where(np.abs(vars(self)[field] - var_dict['_FillValue']) < np.abs(var_dict['_FillValue']*.01),vars(self)[field])
      else:
        vars(self)[field] = np.ma.masked_where(np.abs(vars(self)[field] - var_dict['missing_value']) < np.abs(var_dict['missing_value']*.01),vars(self)[field])             

      var_dict['masked']=True

    
    var_dict['z_interfaces'] = None
       
    if var_dict['Z'] is not None:
      if var_dict['zbax_data'] is not None and use_interfaces==False:
        tmp = np.reshape(var_dict['zbax_data'][z_interfaces],(nz+1,1,1))
        if geo_region is not None:
          ny = len(geo_region['y'])
          nx = len(geo_region['x'])          
        else:
          ny = len(f.variables[var_dict['Y']][:])
          nx = len(f.variables[var_dict['X']][:])

        var_dict['Zdir']=get_axis_direction(f.variables[var_dict['Z']])
        var_dict['z_interfaces']  = var_dict['Zdir']*np.tile(tmp,(1,ny,nx))
        tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))             
        var_dict['z'] = tmp[0:-1,:,:]
        var_dict['dz'] = var_dict['z_interfaces'] - np.roll(var_dict['z_interfaces'],axis=0,shift=-1)
        var_dict['dz'] = var_dict['dz'][0:-1,:,:]
      if var_dict['Z'] is not None and use_interfaces and self.interfaces is not None:           
        var_dict['z_interfaces']  = vars(self)[self.interfaces]
        tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=1,shift=-1))             
        var_dict['z'] = tmp[:,0:-1,:,:]
        var_dict['dz'] = var_dict['z_interfaces'] - np.roll(var_dict['z_interfaces'],axis=1,shift=-1)
        var_dict['dz'] = var_dict['dz'][:,0:-1,:,:]
      if use_interfaces and self.interfaces is None:
          print """ use_interfaces=True and self.interfaces is None in call to add_field"""
          return None
               
    try:
      time_avg_info = getattr(f.variables[field],'time_avg_info')
    except:
      time_avg_info = ""
         
    if 'average_DT' in time_avg_info:
      dt = f.variables['average_DT'][t_indices]
      units = var_dict['rootgrp'].variables['average_DT'].units
      if units == 'days':
        dt = dt*86400. # convert to seconds
      elif units == 'hours':
        dt = dt*3600. # convert to seconds
    else:
      if var_dict['T'] is not None:           
        dt = np.ones((nt))
      else:
        dt = 0.0

    var_dict['dt'] = dt      


    self.var_dict[field] = var_dict

    return None


    
    
  def create_field(self,expression,name,var_dict=None):
    """
    Define a new field and add to the existing state.
    The result of (expression) evaluate to (self.name)
    with associated (var_dict).
    """

    cmd = string.join(['self.',name,'=',expression],sep='')
    exec(cmd)

    self.variables[name]=name

    if var_dict is not None:
      self.var_dict[name]=dict.copy(var_dict)
    else:
      self.var_dict[name] = {}


    self.var_dict[name]['history']=expression
                  
    return None


  def add_field_from_array(self,var,name,var_dict=None,history=None):
    """
    Define a new field from an existing array and add to the current state.
    """

    self.variables[name]=name

    vars(self)[name]=var
    
    if var_dict is not None:
      self.var_dict[name]=dict.copy(var_dict)
    else:
      self.var_dict[name] = {}


      
    self.var_dict[name]['history']=history
                  
    return None
    

  def del_field(self,field):
    """
    Delete (field) 
    """
    
    f = self.rootgrp
    cmd = string.join(['del(self.',field,')'],sep='')
    exec(cmd)
    
    self.variables.pop(field)
    self.var_dict.pop(field)

    return None

  def rename_field(self,field,field_new):
    """
    Rename (field) 
    """
    
    f = self.rootgrp
    cmd = string.join(['self.',field_new,'=vars(self)[\'',field,'\'].copy()'],sep='')    

    exec(cmd)
    
    self.variables.pop(field)
    self.variables[field_new]=field_new
    self.var_dict[field_new] = dict.copy(self.var_dict[field])
    self.var_dict.pop(field)

    return None

  def mask_where(self,field=None,condition=None):
    """
    Mask (field) where (condition)
    """

    
    if condition is None or field is None:
      return None
    else:
      cmd = 'result=self.'+condition
      exec(cmd)
      shape_result = result.shape
      shape = vars(self)[field].shape

      if shape_result[-1] != shape[3] or shape_result[-2] != shape[2]:
        print """
        x-y Shape mismatch in self.mask_where"""
        return None

      if len(shape_result) ==2:
        result=np.reshape(result,(1,1,result.shape[0],result.shape[1]))
        result=np.tile(result,(shape[0],shape[1],1,1))
      elif len(shape_result) == 3:
        result=np.reshape(result,(1,1,result.shape[0],result.shape[1]))
        result=np.tile(result,(shape[0],1,1,1))
      elif len(shape_result) == 4:
        if shape_result[1] == 1:
          result=np.tile(result,(1,shape[1],1,1))                
        elif shape_result[0] == 1 and shape_result[1] == shape[1]:
          result=np.tile(result,(shape[0],1,1,1))
        elif shape_result[0] == 1 and shape_result[1] == 1:
          result=np.tile(result,(shape[0],shape[1],1,1))
        else:
          print """
          Error expanding mask """
          return None
        
      vars(self)[field] = np.ma.masked_where(result,vars(self)[field])
      self.var_dict[field]['masked']=True

  def unmask(self,field=None):
    """
    Remove mask
    """

    if field is None:
      return None

    vars(self)[field].soften_mask() 
    vars(self)[field].mask = False

    
  def add_interface_bounds(self,field=None):
    """
    Add interfaces at cell boundaries (ny+1,nx+1).

    ***NOTE: Assumes cyclic x at this point.
    ***TODO
    
    """
    
    if field is None:
      return None

#    if self.var_dict[field]['Ztype'] is not 'Fixed':
    e=self.var_dict[field]['z_interfaces']
    eb = 0.5*(e+np.roll(e,shift=-1,axis=3))
    self.var_dict[field]['z_interfaces_ew']=np.concatenate((np.take(eb,[-1],axis=3),eb),axis=3)
    eb = 0.5*(e+np.roll(e,shift=-1,axis=2))
    self.var_dict[field]['z_interfaces_ns']=np.concatenate((np.take(eb,[0],axis=2),eb),axis=2)
#    else:
#        e=self.var_dict[field]['z_interfaces']
#        eb = 0.5*(e+np.roll(e,shift=-1,axis=2))
#        self.var_dict[field]['z_interfaces_ew']=np.concatenate((np.take(eb,[-1],axis=2),eb),axis=2)
#        eb = 0.5*(e+np.roll(e,shift=-1,axis=1))
#        self.var_dict[field]['z_interfaces_ns']=np.concatenate((np.take(eb,[0],axis=1),eb),axis=1)            
    
  def fill_interior(self,field=None,smooth=False,num_pass=10):
    """
    Fill interior above the topography .

    >>> from midas import *
    >>> from midas import vertmap
    >>> import hashlib
    >>> x=np.linspace(0.,np.pi,100)
    >>> X,Y=np.meshgrid(x,x)
    >>> sgrid=supergrid(xdat=X,ydat=Y,config='cartesian')
    >>> grid=rectgrid(supergrid=sgrid,cyclic=True)
    >>> X=grid.x_T;Y=grid.y_T
    >>> grid.wet=np.ones(X.shape)
    >>> grid.D=np.ones(X.shape)
    >>> S=state(grid=grid)
    >>> a=np.cos(X)*np.sin(Y)
    >>> a=a[np.newaxis,np.newaxis,:]
    >>> S.add_field_from_array(a,'a')
    >>> S.a[0,0,::5,::5]=-1.e20
    >>> S.a=np.ma.masked_where(S.a==-1.e20,S.a)
    >>> S.var_dict['a']['masked']=True
    >>> print np.ma.sum(S.a)
    -20.1663931523
    >>> S.fill_interior('a')
    >>> print np.ma.sum(S.a)
    -18.8982071632
    """

    from midas import vertmap

    FVal_=-1.e34

  
    if field is None:
      print """
       Please specify a field name
      """
      return None

    if self.var_dict[field]['masked'] is False:
      print """
       Input field needs to be masked in  order to use fill
      """
      return None

    
    if self.var_dict[field]['_FillValue'] is not None:
        FillValue = self.var_dict[field]['_FillValue']          
    elif self.var_dict[field]['missing_value'] is not None:
        FillValue = self.var_dict[field]['missing_value']
    else:
        
        FillValue = FVal_
        self.var_dict[field]['_FillValue'] = FVal_
        self.var_dict[field]['missing_value'] = FVal_                


    val = vars(self)[field].copy()
    val.set_fill_value(FillValue)

    val_prev = np.zeros([val.shape[2],val.shape[3]])
    
    wet = self.grid.wet


    for i in np.arange(0,val.shape[0]):
      if val.shape[1] == 1:
        tmp = np.squeeze(np.take(np.take(val,[i],axis=0),[0],axis=1))
        mask_in = np.ma.getmask(tmp)
        fill = np.zeros([tmp.shape[0],tmp.shape[1]])
        fill[np.logical_and(wet == 1.0,mask_in) ]=1.0        
        mask_out=np.zeros([tmp.shape[0],tmp.shape[1]])
        mask_out[wet == 0.0]=1
        good = np.zeros([tmp.shape[0],tmp.shape[1]])
        good[~mask_in]=1
        v_filled = np.zeros([tmp.shape[1],tmp.shape[0]])
        v_filled=vertmap.midas_vertmap.fill_miss_2d(tmp.T,good.T,fill.T,cyclic_x=self.grid.cyclic_x,tripolar_n=self.grid.tripolar_n,smooth=smooth,num_pass=num_pass)
        v_filled=v_filled.T
        v_filled[mask_out==1]=FillValue
        val[i,0,:]=v_filled[:]
      else:
        for j in np.arange(0,val.shape[1]):
            zbot = self.var_dict[field]['z_interfaces'][j+1,:]
            ztop = self.var_dict[field]['z_interfaces'][j,:]        
            tmp = np.squeeze(np.ma.take(np.ma.take(val,[i],axis=0),[j],axis=1))
            mask_in = np.ma.getmask(tmp)

            # fill has a value of one over points
            # which are in the interior at the current depth
            # and zero otherwise
        
            fill = np.zeros([tmp.shape[0],tmp.shape[1]])
            fill[np.logical_and(np.logical_and(wet == 1.0,-ztop < self.grid.D),mask_in) ]=1

            mask_out=np.logical_or(wet == 0.0,-ztop > self.grid.D)
            
            good = np.zeros([tmp.shape[0],tmp.shape[1]])
            good[~mask_in]=1


            v_filled = np.zeros([tmp.shape[1],tmp.shape[0]])

            if j>0:
                # initialize with nearest fill or value at previous level
                v_filled=vertmap.midas_vertmap.fill_miss_2d(tmp.T,good.T,fill.T,val_prev.T,cyclic_x=self.grid.cyclic_x,tripolar_n=self.grid.tripolar_n,smooth=smooth,num_pass=num_pass)
            else:
                v_filled=vertmap.midas_vertmap.fill_miss_2d(tmp.T,good.T,fill.T,cyclic_x=self.grid.cyclic_x,tripolar_n=self.grid.tripolar_n,smooth=smooth,num_pass=num_pass)

        
            v_filled=v_filled.T
            val_prev = v_filled.copy()


            v_filled[mask_out==1]=FillValue                


            val[i,j,:]=v_filled.copy()


    val=np.ma.masked_where(np.abs(val-FillValue)<1.e-4*np.abs(FillValue),val)

    vars(self)[field]=val

    return None

    
  def volume_integral(self,field=None,axis=None,normalize=True):
    """

    Calculate a finite-volume-weighted integral of (field)
    , axis can be along one or any combination of XYZ axes,
    i.e. 'X','Y','Z','XY','XZ','YZ','XYZ'.
    
    """


    if self.var_dict[field]['interface_variable']:
        is_interface_var = True
    else:
        is_interface_var = False
    
    
    sout=vars(self)[field]
    
    var_dict = dict.copy(self.var_dict[field]) # inherit variable dictionary from parent 

    if self.var_dict[field]['Z'] is not None and not is_interface_var :
      if self.var_dict[field]['Ztype'] is 'Fixed':

        dz=self.var_dict[field]['dz'][:].reshape(1,sout.shape[1],sout.shape[2],sout.shape[3])
        dz = np.tile(dz,(sout.shape[0],1,1,1))


        if self.var_dict[field]['masked']:
          dz_masked = np.ma.masked_where(sout.mask,dz)
        else:
          dz_masked = dz

      else:
        dz=self.var_dict[field]['dz'][:]


        if self.var_dict[field]['masked']:
          dz_masked = np.ma.masked_where(sout.mask,dz)
        else:
          dz_masked = dz
          
    else:
      dz = np.ones((sout.shape))      
      dz_masked = dz

    dy=self.grid.dyh
    dy = np.tile(dy,(sout.shape[0],sout.shape[1],1,1))


    dx=self.grid.dxh
    dx = np.tile(dx,(sout.shape[0],sout.shape[1],1,1))        
        
    if self.var_dict[field]['masked']:
      dy_masked = np.ma.masked_where(sout.mask,dy)
      dx_masked = np.ma.masked_where(sout.mask,dx)
    else:
      dy_masked = dy
      dx_masked = dx

    if axis.upper() == 'Z':

      if var_dict['Z'] is None:
        return None

      if normalize:
          result = np.sum(sout*dz,axis=1)/np.sum(dz_masked,axis=1)
          result=np.reshape(result,(result.shape[0],1,result.shape[1],result.shape[2]))
          name = field+'_zav'
      else:
          result = np.sum(sout*dz,axis=1)
          result=np.reshape(result,(result.shape[0],1,result.shape[1],result.shape[2]))
          name = field+'_zint'
          
      vars(self)[name]=result

      
      var_dict['rootgrp']=None
            
      var_dict['z_indices']=[0]

      if var_dict['zax_data'] is not None:
          var_dict['zax_data']=[np.mean(var_dict['zax_data'])]
          var_dict['zbax_data']=[var_dict['zbax_data'][0],var_dict['zbax_data'][-1]]
      else:
          var_dict['zax_data']=[0.0]
          var_dict['zbax_data']=[0.0,1.0]                    

      var_dict['Z']=None



      if self.var_dict[field]['Ztype'] is 'Fixed':
          zlow=np.take(var_dict['z_interfaces'],[-1],axis=0)                
          zup=np.take(var_dict['z_interfaces'],[0],axis=0)
          var_dict['z_interfaces']=np.concatenate((zup,zlow),axis=0)
          result=var_dict['z_interfaces'][0,:]-var_dict['z_interfaces'][-1,:]
          result=np.reshape(result,(1,1,result.shape[0],result.shape[1]))          
      else:
          zlow=np.take(var_dict['z_interfaces'],[-1],axis=1)                
          zup=np.take(var_dict['z_interfaces'],[0],axis=1)
          var_dict['z_interfaces']=np.concatenate((zup,zlow),axis=1)
          result=var_dict['z_interfaces'][:,0,:]-var_dict['z_interfaces'][:,-1,:]
          result=np.reshape(result,(result.shape[0],1,result.shape[1],result.shape[2]))          

      var_dict['dz']=result
      var_dict['z']=np.mean(var_dict['z_interfaces'],axis=1)
      var_dict['z_interfaces']=np.reshape(var_dict['z_interfaces'],(result.shape[0],2,result.shape[2],result.shape[3]))

      if var_dict['Ztype'] is 'Fixed':
        var_dict['z']=np.take(var_dict['z'],[0],axis=0)
        var_dict['dz']=np.take(var_dict['dz'],[0],axis=0)
        var_dict['z_interfaces']=np.take(var_dict['z_interfaces'],[0],axis=1)                
      

      self.var_dict[name]=var_dict
      self.variables[name]=name
      
        
    if axis.upper() == 'Y':

      if normalize:
          result = np.sum(sout*dz_masked*dy_masked*dx_masked,axis=2)/np.sum(dz_masked*dy_masked*dx_masked,axis=2)
          result=np.reshape(result,(sout.shape[0],sout.shape[1],1,sout.shape[3]))
          name = field+'_yav'
      else:
          result = np.sum(sout*dz_masked*dy_masked*dx_masked,axis=2)
          result=np.reshape(result,(sout.shape[0],sout.shape[1],1,sout.shape[3]))
          name = field+'_yint'
          
      vars(self)[name]=result



      if var_dict['Z'] is not None:
        var_dict['dz'] = np.sum(dz_masked*dz_masked*dy_masked*dx_masked,axis=2)/np.sum(dz_masked*dy_masked*dx_masked,axis=2)      
        var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],sout.shape[1],1,sout.shape[3]))      
        if var_dict['Ztype'] is 'Fixed':
            z0 = np.take(var_dict['z_interfaces'],[0],axis=0)
        else:
            z0 = np.take(var_dict['z_interfaces'],[0],axis=1)
        dy0 = np.take(dy_masked,[0],axis=1)
        dx0 = np.take(dx_masked,[0],axis=1)
        result = np.sum(z0*dy0*dx0,axis=2)/np.sum(dy0*dx0,axis=2)
        result = np.reshape(result,(sout.shape[0],1,1,sout.shape[3]))
        var_dict['z_interfaces']=-np.cumsum(var_dict['dz'],axis=1)
        var_dict['z_interfaces']=np.concatenate((result,var_dict['z_interfaces']),axis=1)


        tmp=0.5*(var_dict['z_interfaces'] + np.roll(var_dict['z_interfaces'],axis=1,shift=-1))
        tmp = tmp[:,0:-1,:,:]

        if var_dict['Ztype'] is 'Fixed':
          var_dict['z']=np.take(tmp,[0],axis=0)
          var_dict['dz']=np.take(var_dict['dz'],[0],axis=0)
          var_dict['z_interfaces']=np.take(var_dict['z_interfaces'],[0],axis=0)                
        else:
          var_dict['z']=tmp
        
      var_dict['rootgrp']=None
      var_dict['Y']=None
      var_dict['yax_data']=[np.mean(var_dict['yax_data'])]
      var_dict['rootgrp']=None
      

      self.var_dict[name]=var_dict
      self.variables[name]=name

        
    if axis.upper() == 'X':

      if normalize:
          result = np.sum(sout*dz_masked*dy_masked*dx_masked,axis=3)/np.sum(dz_masked*dy_masked*dx_masked,axis=3)
          result=np.reshape(result,(sout.shape[0],sout.shape[1],sout.shape[2],1))
          name = field+'_xav'
      else:
          result = np.sum(sout*dz_masked*dy_masked*dx_masked,axis=3)
          result=np.reshape(result,(sout.shape[0],sout.shape[1],sout.shape[2],1))
          name = field+'_xint'
          
      vars(self)[name]=result


      if var_dict['Z'] is not None:
        
        var_dict['dz'] = np.sum(dz_masked*dz_masked*dy_masked*dx_masked,axis=3)/np.sum(dz_masked*dy_masked*dx_masked,axis=3)      
        var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],sout.shape[1],sout.shape[2],1))
        if var_dict['Ztype'] is 'Fixed':
            z0 = np.take(var_dict['z_interfaces'],[0],axis=0)
        else:
            z0 = np.take(var_dict['z_interfaces'],[0],axis=1)
        dy0 = np.take(dy_masked,[0],axis=1)
        dx0 = np.take(dx_masked,[0],axis=1)
        result = np.sum(z0*dy0*dx0,axis=3)/np.sum(dy0*dx0,axis=3)
        result = np.reshape(result,(sout.shape[0],1,sout.shape[2],1))
        var_dict['z_interfaces']=-np.cumsum(var_dict['dz'],axis=1)
        var_dict['z_interfaces']=np.concatenate((result,var_dict['z_interfaces']),axis=1)
        tmp=0.5*(var_dict['z_interfaces'] + np.roll(var_dict['z_interfaces'],axis=1,shift=-1))
        tmp = tmp[:,0:-1,:,:]

        if var_dict['Ztype'] is 'Fixed':
          var_dict['z']=np.take(tmp,[0],axis=0)
          var_dict['dz']=np.take(var_dict['dz'],[0],axis=0)
          var_dict['z_interfaces']=np.take(var_dict['z_interfaces'],[0],axis=0)                
        else:
          var_dict['z']=tmp
        


      var_dict['rootgrp']=None
      var_dict['X']=None
      var_dict['xax_data']=[np.mean(var_dict['xax_data'])]
      var_dict['rootgrp']=None
      

      self.var_dict[name]=var_dict
      self.variables[name]=name
      

    if axis.upper() == 'XY' or axis.upper() == 'YX':

      if normalize:
          result = np.sum(np.sum(sout*dz_masked*dy_masked*dx_masked,axis=3),axis=2)/np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=2)
          result=np.reshape(result,(sout.shape[0],sout.shape[1],1,1))
          name = field+'_xyav'
      else:
          result = np.sum(np.sum(sout*dz_masked*dy_masked*dx_masked,axis=3),axis=2)
          result=np.reshape(result,(sout.shape[0],sout.shape[1],1,1))
          name = field+'_xyint'

          
      vars(self)[name]=result

      if var_dict['Z'] is not None:      
        var_dict['dz'] = np.sum(np.sum(dz_masked*dz_masked*dy_masked*dx_masked,axis=3),axis=2)/np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=2)
        var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],sout.shape[1],1,1))
        if not is_interface_var :
            result = np.sum(np.sum(dz_masked*dz_masked*dy_masked*dx_masked,axis=3),axis=2)/np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=2)
            result = np.reshape(result,(sout.shape[0],sout.shape[1],1,1))
            var_dict['dz']=result.copy()
            zi = -np.cumsum(result,axis=1)
            if var_dict['Ztype'] is 'Fixed':
                z0 = [0]
                z0 = np.reshape(z0,(1,1,1,1))
                zi=np.take(zi,[0],axis=0)
                zi=np.concatenate((z0,zi),axis=1)
                var_dict['z_interfaces']=zi.copy()
                tmp=0.5*(zi + np.roll(zi,axis=0,shift=-1))
                tmp = tmp[0:-1,:,:]
                var_dict['z']=tmp                
            else:
                z0 = [0]
                z0 = np.reshape(z0,(1,1,1,1))                
                zi=np.concatenate((z0,zi),axis=1)            
                var_dict['z_interfaces']=zi.copy()
                tmp=0.5*(zi + np.roll(zi,axis=1,shift=-1))
                tmp = tmp[:,0:-1,:,:]                
                var_dict['z']=tmp

        else:
            var_dict['z_interfaces']=None
            


      var_dict['X']=None
      var_dict['Y']=None      
      var_dict['xax_data']=[np.mean(var_dict['xax_data'])]
      var_dict['yax_data']=[np.mean(var_dict['yax_data'])]      
      var_dict['rootgrp']=None
      

      self.var_dict[name]=var_dict
      self.variables[name]=name
      

      
    if axis.upper() == 'XZ' or axis.upper() == 'ZX':

      if var_dict['Z'] is None or var_dict['X'] is None:
        return None

      if normalize:
          result = np.sum(np.sum(sout*dz_masked*dy_masked*dx_masked,axis=3),axis=1)/np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=1)
          result=np.reshape(result,(result.shape[0],1,result.shape[1],1))
          name = field+'_xzav'
      else:
          result = np.sum(np.sum(sout*dz_masked*dy_masked*dx_masked,axis=3),axis=1)
          result=np.reshape(result,(result.shape[0],1,result.shape[1],1))
          name = field+'_xzint'
          
      vars(self)[name]=result
      
      var_dict['dz'] = np.sum(np.sum(dz_masked*dz_masked*dy_masked*dx_masked,axis=3),axis=1)/np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=1)
      var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],1,sout.shape[2],1))            

      z0 = np.take(var_dict['z_interfaces'],[0],axis=1)
      dy0 = np.take(dy_masked,[0],axis=1)
      dx0 = np.take(dx_masked,[0],axis=1)      
      result = np.sum(z0*dy0*dx0,axis=3)/np.sum(dy0*dx0,axis=3)
      result = result[:,:,:,np.newaxis]

      var_dict['z_interfaces']=-np.reshape(np.sum(var_dict['dz'],axis=1),(sout.shape[0],1,sout.shape[2],1))
    
      var_dict['z_interfaces']=np.concatenate((result,var_dict['z_interfaces']),axis=1)

      var_dict['z_indices']=[0]
      var_dict['zax_data']=[np.mean(var_dict['zax_data'])]
      var_dict['zbax_data']=[var_dict['zbax_data'][0],var_dict['zbax_data'][-1]]
      
      if var_dict['Ztype'] == 'Fixed':
        var_dict['z_interfaces']=np.take(var_dict['z_interfaces'],[0],axis=0)

      var_dict['Z']=None
      var_dict['X']=None
      var_dict['xax_data']=[np.mean(var_dict['xax_data'])]
      var_dict['rootgrp']=None
      

      self.var_dict[name]=var_dict
      self.variables[name]=name
      

    if axis.upper() == 'YZ' or axis.upper() == 'ZY':

      if var_dict['Z'] is None or var_dict['Y'] is None:
        return None

      if normalize:
          result = np.sum(np.sum(sout*dx_masked*dy_masked*dz_masked,axis=2),axis=1)/np.sum(np.sum(dx_masked*dy_masked*dz_masked,axis=2),axis=1)
          result=np.reshape(result,(result.shape[0],1,1,result.shape[1]))
          name = field+'_yzav'
      else:
          result = np.sum(np.sum(sout*dx_masked*dy_masked*dz_masked,axis=2),axis=1)
          result=np.reshape(result,(result.shape[0],1,1,result.shape[1]))
          name = field+'_yzint'
          
      vars(self)[name]=result
      
      var_dict['dz'] = np.sum(np.sum(dz_masked*dz_masked*dy_masked*dx_masked,axis=2),axis=1)/np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=2),axis=1)
      var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],1,1,sout.shape[3]))                  

      z0 = np.take(var_dict['z_interfaces'],[0],axis=1)
      dy0 = np.take(dy_masked,[0],axis=1)
      dx0 = np.take(dx_masked,[0],axis=1)      
      result = np.sum(z0*dy0*dx0,axis=2)/np.sum(dy0*dx0,axis=2)
      result = np.reshape(result,(sout.shape[0],1,1,sout.shape[3]))
      var_dict['z_interfaces']=-np.reshape(np.sum(var_dict['dz'],axis=1),(sout.shape[0],1,1,sout.shape[3]))
    
      var_dict['z_interfaces']=np.concatenate((result,var_dict['z_interfaces']),axis=1)

      var_dict['z_indices']=[0]
      var_dict['zax_data']=[np.mean(var_dict['zax_data'])]
      var_dict['zbax_data']=[var_dict['zbax_data'][0],var_dict['zbax_data'][-1]]
      
      if var_dict['Ztype'] == 'Fixed':
        var_dict['z_interfaces']=np.take(var_dict['z_interfaces'],[0],axis=0)

      var_dict['Z']=None
      var_dict['Y']=None
      var_dict['yax_data']=[np.mean(var_dict['yax_data'])]      
      var_dict['rootgrp']=None
      

      self.var_dict[name]=var_dict
      self.variables[name]=name

    if axis.upper() == 'XYZ':

      if var_dict['Z'] is None or var_dict['Y'] is None or var_dict['X'] is None:
        return None      

      if normalize:
          result = np.sum(np.sum(np.sum(sout*dx_masked*dy_masked*dz_masked,axis=3),axis=2),axis=1)/np.sum(np.sum(np.sum(dx_masked*dy_masked*dz_masked,axis=3),axis=2),axis=1)
          result=np.reshape(result,(result.shape[0],1,1,1))
          name = field+'_xyzav'
      else:
          result = np.sum(np.sum(np.sum(sout*dx_masked*dy_masked*dz_masked,axis=3),axis=2),axis=1)
          result=np.reshape(result,(result.shape[0],1,1,1))
          name = field+'_xyzint'
          
      vars(self)[name]=result
      
        
      var_dict['Z']=None
      var_dict['X']=None
      var_dict['xax_data']=[np.mean(var_dict['xax_data'])]
      var_dict['Y']=None
      var_dict['yax_data']=[np.mean(var_dict['yax_data'])]      
      var_dict['rootgrp']=None
      

      self.var_dict[name]=var_dict
      self.variables[name]=name

  def time_avg(self,field=None,vol_weight=True,target=None):
    """

    Calculate a finite-volume-weighted average of (field)
    along the time dimension.
    
    """

    missing=-1.e34
    
    if self.var_dict[field]['T'] is None:
      return None
      

    cmd = string.join(['sout=self.',field],sep='')
    exec(cmd)

    var_dict = dict.copy(self.var_dict[field]) # inherit variable dictionary from parent 


    if target is not None:
        nt = len(target['tax_data'])
    else:
        nt = 1
        
    if vol_weight == True:
        if self.var_dict[field]['Z'] is not None:
            if self.var_dict[field]['Ztype'] is 'Fixed':
                
                w=self.var_dict[field]['dz'][:]
                w = np.tile(w,(sout.shape[0],1,1,1))
        
                if self.var_dict[field]['masked']:
                    w_masked = np.ma.masked_where(sout.mask,w)
                else:
                    w_masked = w
            else:
                w=self.var_dict[field]['dz'][:]
          
                if self.var_dict[field]['masked']:
                    w_masked = np.ma.masked_where(sout.mask,w)
                else:
                    w_masked = w
        else:
            w = np.ones((sout.shape))

            if self.var_dict[field]['masked']:
                w_masked = np.ma.masked_where(sout.mask,w)
    else:
        w = np.ones((sout.shape))      
        w_masked = w

        if self.var_dict[field]['masked']:
            w_masked = np.ma.masked_where(sout.mask,w)        
        
            
    dt=self.var_dict[field]['dt'][:]

    shape=sout.shape
    w_masked = w_masked*np.tile(np.reshape(dt,(shape[0],1,1,1)),(1,shape[1],shape[2],shape[3]))
    mask_in = np.zeros((w_masked.shape[0],1,w_masked.shape[2],w_masked.shape[3]))
    tmp = w_masked[:,0,:]
    tmp = tmp[:,np.newaxis,:]
    mask_in[tmp.mask==False]=1.0
    result = np.ma.zeros((np.hstack((nt,sout.shape[1:]))))
    result[:]=missing
    result = np.ma.zeros((np.hstack((nt,sout.shape[1:]))))
    ntimes = np.ma.zeros((np.hstack((1,1,sout.shape[2:]))))
    ntimes[:]=missing
    nzp=shape[1]+1

    if var_dict['Z'] is not 'Fixed':
        interfaces= np.ma.zeros((np.hstack((nt,nzp,shape[2:]))))
    else:
        interfaces = np.ma.zeros((np.hstack((nzp,shape[2:]))))

    interfaces[:]=missing
    avg_time = 0.


    for j in np.arange(0,nt):

      if target is not None:
          ts,te = find_date_bounds(self.var_dict[field]['dates'],target['date_bounds'][j],target['date_bounds'][j+1])
          print 'j,ts,te= ',j,ts,te
          
          if ts != -1 and te != -1:
              if self.var_dict[field]['masked']:              
                  result[j,:,:,:]=np.ma.sum(sout[ts:te,:,:,:]*w_masked[ts:te,:,:,:],axis=0)/np.ma.sum(w_masked[ts:te,:,:,:],axis=0)
              else:
                  result[j,:,:,:]=np.sum(sout[ts:te,:,:,:]*w_masked[ts:te,:,:,:],axis=0)/np.sum(w_masked[ts:te,:,:,:],axis=0)

              if var_dict['z_interfaces'] is not None:
                  if var_dict['Ztype'] is not 'Fixed':
                      if self.var_dict[field]['masked']:
                          interfaces[j,:,:,:]=np.ma.sum(var_dict['z_interfaces'][ts:te,:,:,:]*dt[ts:te],axis=0)/np.ma.sum(dt[ts:te])                          
                      else:
                          interfaces[j,:,:,:]=np.sum(var_dict['z_interfaces'][ts:te,:,:,:]*dt[ts:te],axis=0)/np.sum(dt[ts:te])
      else:
          if self.var_dict[field]['masked']:          
              result[j,:,:,:]=np.ma.sum(sout*w_masked,axis=0)/np.ma.sum(w_masked,axis=0)
          else:
              result[j,:,:,:]=np.sum(sout*w_masked,axis=0)/np.sum(w_masked,axis=0)              
          
          result = np.ma.masked_where(result == missing, result)

          if var_dict['Z'] is not None:
              if var_dict['Ztype'] is 'Fixed':
                  if var_dict['z_interfaces'] is not None:        
                      interfaces=var_dict['z_interfaces'][:]
                      dz=np.ma.zeros(np.hstack(sout.shape[1:]))
                      for k in np.arange(0,sout.shape[1]):
                          dz[k,:,:]=interfaces[k,:,:]-interfaces[k+1,:,:]        
              else:
                  dz = np.ma.zeros((np.hstack((nt,sout.shape[1:]))))

              if var_dict['z_interfaces'] is not None:
                  if var_dict['Ztype'] is not 'Fixed':        
                      for k in np.arange(0,sout.shape[1]):
                          dz[:,k,:,:]=interfaces[:,k,:,:]-interfaces[:,k+1,:,:]
                  else:
                      for k in np.arange(0,sout.shape[1]):
                          dz[k,:,:]=interfaces[k,:,:]-interfaces[k+1,:,:]

    if var_dict['Z'] is not None:
        if var_dict['Ztype'] is not 'Fixed':        
            var_dict['z_interfaces']=interfaces
            var_dict['dz']=np.squeeze(dz)
            tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))
            var_dict['z'] = tmp[:,0:-1,:,:]
            var_dict['z'] = np.reshape(var_dict['z'],(nt,shape[1],shape[2],shape[3]))
            var_dict['dz'] = np.reshape(var_dict['dz'],(nt,shape[1],shape[2],shape[3]))            

        else:
            var_dict['z_interfaces']=sq(interfaces)
            dz=np.ma.zeros(np.hstack(sout.shape[1:]))
            for k in np.arange(0,sout.shape[1]):
                dz[k,:,:]=interfaces[k,:,:]-interfaces[k+1,:,:]                    
            var_dict['dz']=dz
            tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))
            var_dict['z'] = tmp[0:-1,:,:]
            var_dict['z'] = np.reshape(var_dict['z'],(shape[1],shape[2],shape[3]))
            var_dict['dz'] = np.reshape(var_dict['dz'],(shape[1],shape[2],shape[3]))

    if target is not None:
        var_dict['tbax_data'] = target['tbax_data']
        var_dict['tax_data'] = target['tax_data']        
        date_bounds = target['date_bounds']
        var_dict['dates'] = target['dates']
        var_dict['T'] = 'time'
        var_dict['dt']=target['dt']

    else:
        tbax_data = [var_dict['tbax_data'][0],var_dict['tbax_data'][-1]]
        var_dict['tbax_data']=tbax_data
        date_bounds = [var_dict['date_bounds'][0],var_dict['date_bounds'][-1]]
        avg_time = np.sum(dt*var_dict['tax_data'])/np.sum(dt)
        var_dict['dates'] = num2date(avg_time,var_dict['tunits'],var_dict['calendar'])
        var_dict['T'] = 'time'
        var_dict['dt']=[np.sum(dt)]      
        var_dict['tax_data'] = avg_time    


    var_dict['rootgrp']=None

    
    

    name = field+'_tav'
    vars(self)[name]=result
    self.var_dict[name]=var_dict.copy()
    self.variables[name]=name

    ntimes=np.sum(mask_in,axis=0)
    ntimes=ntimes[np.newaxis,:]

    var_dict['T']=None
    var_dict['Z']=None    
    
    name = field+'_nsamp'
    vars(self)[name]=ntimes
    self.var_dict[name]=var_dict.copy()
    self.variables[name]=name    

  def monthly_avg(self,field=None,vol_weight=True, DEBUG=False):
    """

    Calculate a finite-volume-weighted average of (field)
    along the time dimension for each month.
    
    """

    if self.var_dict[field]['T'] is None:
      return None
      

    cmd = string.join(['sout=self.',field],sep='')
    exec(cmd)

    var_dict = dict.copy(self.var_dict[field]) # inherit variable dictionary from parent 

    if vol_weight == True:
        #
        # weights are equal to the thickness
        #
        if self.var_dict[field]['Z'] is not None:
            if self.var_dict[field]['Ztype'] is 'Fixed':

                dz=self.var_dict[field]['dz'][:]
                w = np.tile(dz,(sout.shape[0],1,1,1))
        
                if self.var_dict[field]['masked']:
                    w_masked = np.ma.masked_where(sout.mask,w)
                else:
                    w_masked = w
            else:
                w=self.var_dict[field]['dz'][:]
          
                if self.var_dict[field]['masked']:
                    w_masked = np.ma.masked_where(sout.mask,w)
                else:
                    w_masked = w

        else:
            # weights are unity
            w = np.ones((sout.shape))

            if self.var_dict[field]['masked']:
                w_masked = np.ma.masked_where(sout.mask,w)
    else:
        w = np.ones((sout.shape))      
        w_masked = w

        if self.var_dict[field]['masked']:
            w_masked = np.ma.masked_where(sout.mask,w)
            
    dt=self.var_dict[field]['dt'][:]


    shape=sout.shape
    result = np.ma.zeros((np.hstack((12,sout.shape[1:]))))
    nzp=shape[1]+1
    interfaces = np.ma.zeros((np.hstack((12,nzp,shape[2:]))))
    nsamp=np.zeros((result.shape[0],1,result.shape[2],result.shape[3]))
    months=get_months(self.var_dict[field]['dates'])

    weights=np.zeros((12,shape[1],shape[2],shape[3]))

    for i in  np.arange(0,sout.shape[0]):
      tmp=sout[i,:]
      result_ptr=result[months[i]-1,:]
      w_ptr=w_masked[i,:]
      result_ptr[tmp.mask==False]=tmp[tmp.mask==False]*w_ptr[tmp.mask==False] + result_ptr[tmp.mask==False]
      wgt_ptr=weights[months[i]-1,:]
      wgt_ptr[tmp.mask==False]=wgt_ptr[tmp.mask==False]+w_ptr[tmp.mask==False]

      if var_dict['Z'] is not None:
          if var_dict['Ztype'] is not 'Fixed' and var_dict['z_interfaces'] is not None:
              int_ptr=interfaces[months[i]-1,:]
              zi=var_dict['z_interfaces'][i,:]
              int_ptr[tmp.mask==False]=zi[tmp.mask==False]+int_ptr[tmp.mask==False]

      tmp = w_masked[i,:]
      tmp[tmp.mask==False]=1.0
      tmp=np.ma.filled(tmp,0.)
      nsamp[months[i]-1,:]=nsamp[months[i]-1,:] + tmp


    weights = np.ma.masked_where(weights==0.,weights)

    if DEBUG:
        print 'weights=',weights
        print 'result 001=',result
        
    result = result / weights

    if DEBUG:
        print 'result 002=',result
        
    if self.var_dict[field]['masked']:
        result = np.ma.masked_where(nsamp==0.,result)

    elif np.min(nsamp == 0.):
        result = np.ma.masked_where(nsap==0.,result)
    
    if var_dict['Z'] is not None:
        if var_dict['Ztype'] is not 'Fixed' and var_dict['z_interfaces'] is not None:
            interfaces = interfaces/nsamp
            interfaces = np.ma.masked_where(nsamp==0.,interfaces)
      



    if var_dict['Z'] is not None:                
        if var_dict['Ztype'] is not 'Fixed' and var_dict['z_interfaces'] is not None:
            dz = np.ma.zeros((np.hstack((12,sout.shape[1:]))))
            
        for i in np.arange(0,12):
            for k in np.arange(0,sout.shape[1]):
                dz[i,k,:,:]=interfaces[i,k,:,:]-interfaces[i,k+1,:,:]

        if var_dict['Ztype'] is not 'Fixed':
            var_dict['z_interfaces']=interfaces
            var_dict['dz']=dz
            
        if self.var_dict[field]['Ztype'] is 'Fixed':
            tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))
            var_dict['z'] = tmp[0:-1,:,:]
        else:
            tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=1,shift=-1))
            var_dict['z'] = tmp[:,0:-1,:,:]          

    mod_yr=0001
    dates=make_monthly_axis(mod_yr)
    date0=datetime(mod_yr,1,1)

    var_dict['tax_data']=[]
    for i in np.arange(0,12):
        tdelta=dates[i]-date0
        var_dict['tax_data'].append(tdelta.days)
        
    var_dict['dates'] = dates
    var_dict['T'] = 'time'

    var_dict['rootgrp']=None

    
    name = field+'_monthly'
    vars(self)[name]=result
    self.var_dict[name]=var_dict.copy()
    self.variables[name]=name      

    name = field+'_monthly_nsamp'
    vars(self)[name]=nsamp
    var_dict['Z']=None
    self.var_dict[name]=var_dict
    self.variables[name]=name      
    

  def time_interp(self,field=None,target=None,vol_weight=True,name=None):
    """

    Interpolate (field) from its time axis to the (target) axis
    using linear interpolation.  If (vol_weight) is True, the
    cell volume is interpolated as well.

    (target) is a midas dictionary
    
    """

    if name is None:
        name = field+'_tinterp'
    
    if self.var_dict[field]['T'] is None:
      return None
      
    if target is None:
        return None
    
    cmd = string.join(['sout=self.',field],sep='')
    exec(cmd)

    var_dict = dict.copy(self.var_dict[field]) # inherit variable dictionary from parent 

    do_interfaces=False
    do_vol=False

  
    if vol_weight == True:
        if self.var_dict[field]['Z'] is not None:
            do_vol=True
            if self.var_dict[field]['Ztype'] is 'Fixed':

                dz=self.var_dict[field]['dz'][:]
                w = np.tile(dz,(sout.shape[0],1,1,1))
        
                if self.var_dict[field]['masked']:
                    w_masked = np.ma.masked_where(sout.mask,w)
                else:
                    w_masked = w
            else:
                do_interfaces=True
                w=self.var_dict[field]['dz'][:]
          
                if self.var_dict[field]['masked']:
                    w_masked = np.ma.masked_where(sout.mask,w)
                else:
                    w_masked = w
        else:
            w = np.ones((sout.shape))      
            w_masked = w
    else:
        w = np.ones((sout.shape))      
        w_masked = w

        
    shape=sout.shape
    nt = len(target['dates'])

    result = np.ma.zeros((np.hstack((nt,sout.shape[1:]))))
    nzp=shape[1]+1

    if do_interfaces:
        result2 = np.ma.zeros((np.hstack((nt,nzp,shape[2:]))))
        

    t1,t2,w1,w2 = time_interp_weights(self.var_dict[field]['dates'],target['dates'])

    for i in  np.arange(0,nt):

        if w1[i]>1.0:
            result[i,:]=-999.
        else:
            result[i,:,:,:]=sout[t1[i],:,:,:]*w_masked[t1[i],:,:,:]*w1[i] + sout[t2[i],:,:,:]*w_masked[t2[i],:,:,:]*w2[i]

        if do_interfaces:
            if w1[i]>1.0:
                result2[i,:]=-999.
            else:
                result2[i,:,:,:]=var_dict['z_interfaces'][t1[i],:,:,:]*w1[i] + var_dict['z_interfaces'][t2[i],:,:,:]*w2[i]
                    


    if do_interfaces:
        dz = np.ma.zeros((np.hstack((nt,sout.shape[1:]))))        
        for i in np.arange(0,nt):
            for k in np.arange(0,sout.shape[1]):
                dz[i,k,:,:]=result2[i,k,:,:]-result2[i,k+1,:,:]

        var_dict['z_interfaces']=result2
        var_dict['dz']=dz
        tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=1,shift=-1))
        var_dict['z'] = tmp[:,0:-1,:,:]          

    var_dict['tax_data']=target['tax_data']
    var_dict['tunits']=target['tunits']
    var_dict['dates'] = target['dates']
    var_dict['T'] = target['T']

    self.var_dict[name]=var_dict
    self.variables[name]=name      
    vars(self)[name]=np.ma.masked_where(result==-999.,result)
    
  def monthly_anom(self,field=None):
    """

    Calculate anomaly of (field) with 
    respect to a monthly climatology for the corresponding month
    
    """

    if self.var_dict[field]['T'] is None:
      return None
      

    cmd = string.join(['sout=self.',field],sep='')
    exec(cmd)

    var_dict = dict.copy(self.var_dict[field]) # inherit variable dictionary from parent 

    clim_name = field+'_monthly'
    if clim_name not in self.var_dict.keys():
        self.monthly_avg(field)

    cmd = string.join(['climout=self.',clim_name],sep='')
    exec(cmd)
    
    shape=sout.shape
    result = np.ma.zeros((sout.shape))
    months=get_months(self.var_dict[field]['dates'])
    
    for i in  np.arange(0,sout.shape[0]):
      result[i,:,:,:]=sout[i,:,:,:] - climout[months[i]-1,:,:,:]


    var_dict['rootgrp']=None

    
    name = field+'_monthly_anom'
    vars(self)[name]=result
    self.var_dict[name]=var_dict
    self.variables[name]=name      

      
  def time_smooth(self,field,window_len=None,window=None):
    import scipy.signal as signal
    """
    smooth the data in time using a window with requested size and shape.

    input:
    window_len: the dimension of the smoothing window; should be an odd integer
    window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    flat window will produce a moving average smoothing.
    
  
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    signal.lfilter

  """

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
      print """
      Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
      """
      return None


    s=vars(self)[field]
  
    if window == 'flat': 
      w=np.ones(window_len,'d')
    else:
      w=eval('np.'+window+'(window_len)')

    w = np.reshape(w/w.sum(),(window_len,1,1,1))

    W=np.tile(w,(1,s.shape[1],s.shape[2],s.shape[3]))

    result=signal.fftconvolve(s,W,'same')


    result=np.ma.masked_array(result)
    result[0:window_len/2,:].mask=True
    result[-window_len/2:,:].mask=True

    result[0:window_len/2,:]=np.nan
    result[-window_len/2:,:]=np.nan
    
    cmd = 'self.%(f)s_sm_%(l)s = result'%{'f':field,'l':window_len}
    exec(cmd)
    
    name='%(f)s_sm_%(l)s'%{'f':field,'l':window_len}
    self.variables[name]=name

    self.var_dict[name]=dict.copy(self.var_dict[field])

    return None

  def remap_to_isopycnals(self,temp_name='temp',salt_name='salt',R=None,p_ref=2.e7,wet=None,nkml=None,nkbl=None,hml=None,fit_target=False):
    """

    Remap z or z-like data from geopotential surfaces to a potential
    density coordinate with respect to reference pressure (p_ref) [Pa]
    The default reference depth is about 2km in the ocean.

    (temp) and (salt) are used to calculate potential density on the
    input data\'s native vertical grid.  T/S values are extrapolated from
    the input data along level surfaces above the level of topography prior
    to calculating density. Optional (nkml,nkbl,hml) where a finite thickness
    surface layer is preferred (i.e. to initialize a model with a bulk mixed
    layer). Specifying (fit_target) adjusts remapped layer T/S to closely
    match the target density for non-vanishing layers.
    
    
    """

    from midas import vertmap 

    if self.var_dict[temp_name]['Ztype'] is not 'Fixed':
      print """
        Need to provide z-space temperature and salinity in call to remap_Z_to_layers
        """
      return None

    print '...Filling interior points in z-space'
    
    self.fill_interior(temp_name)
    self.fill_interior(salt_name)

    
    self.create_field('wright_eos(vars(self)[\''+temp_name+'\'],vars(self)[\''+salt_name+'\'],'+str(p_ref)+')','sigma',var_dict=self.var_dict[temp_name])

    
    Rb=np.zeros(len(R)+1)
    Rb[0]=0.0
    Rb[1:-1]=0.5*(R[0:-1]+R[1:])
    Rb[-1]=2.0*Rb[-2]

    nlevs=self.sigma.count(axis=1)

    if nlevs.ndim > 2:
        nlevs=nlevs[0,:]
    nlevs=np.squeeze(nlevs.T)
    depth=self.grid.D.T

    zax=self.var_dict[temp_name]['zax_data']
    zbax=self.var_dict[temp_name]['zbax_data']
    
    nt=self.sigma.shape[0];nz=self.sigma.shape[1];ny=self.sigma.shape[2];nx=self.sigma.shape[3]

    if wet is not None:
      wet_=wet.T
    else:
      wet_=np.ones([nx,ny])


    zax=zax.astype('float64')
    Rb=Rb.astype('float64')
    depth=depth.astype('float64')

    land_fill = self.var_dict[temp_name]['missing_value']

    nlay=len(R)
    temp=np.zeros((nx,ny,nlay,nt))
    salt=np.zeros((nx,ny,nlay,nt))
    zi=np.zeros((nx,ny,nlay+1,nt))        
    R=R.astype('float64')
    
    for n in np.arange(nt):
        rho=self.sigma[n,:].T
        zi_n = zi[:,n]
        print '...Finding interface positions '
        zi_n=vertmap.midas_vertmap.find_interfaces(rho,zax,Rb,depth,nlevs,nkml,nkbl,hml)
        zi_n[:,:,1:][zi_n[:,:,1:]>-hml]=-hml
        ptemp=vars(self)[temp_name][n,:].T
        salinity=vars(self)[salt_name][n,:].T
        print '...Remapping temperature '
        temp_n=vertmap.midas_vertmap.tracer_z_init(ptemp,-zbax,zi_n,nkml,nkbl,land_fill,wet_,len(R),nlevs)
        temp_n=temp_n.astype('float64')
        print '...Remapping salinity '
        salt_n=vertmap.midas_vertmap.tracer_z_init(salinity,-zbax,zi_n,nkml,nkbl,land_fill,wet_,len(R),nlevs)
        salt_n=salt_n.astype('float64')        
        h_n=zi_n-np.roll(zi_n,axis=2,shift=-1)
        h_n=h_n.astype('float64')        
        if fit_target:
            print '...Adjusting temp/salt to fit target densities '
            vertmap.midas_vertmap.determine_temperature(temp_n,salt_n,R,p_ref,10,land_fill,h_n,nkml+nkbl+1)
        zi[:,:,:,n]=zi_n
        temp[:,:,:,n]=temp_n
        salt[:,:,:,n]=salt_n
        
    h=zi-np.roll(zi,axis=2,shift=-1)
    h=h[:,:,0:-1,:]
    
#    if fields is not None:
#        for fld in fields:
#            expr = fld+'=vars(self)[\''+fld+'\'].T'
#            exec(expr)
#            expr=fld+'_remap=vertmap.midas_vertmap.tracer_z_init('+fld+',-zbax,zi,nkml,nkbl,land_fill,wet_,len(R),nlevs)'
#            exec(expr)
#            expr=fld+'_remap='+fld+'_remap.astype(\'float64\')'
#            exec(expr)

            

    

    temp=np.ma.masked_where(np.abs(temp-land_fill)<1.e-5,temp)
    salt=np.ma.masked_where(np.abs(salt-land_fill)<1.e-5,salt)

    temp=temp.T
    salt=salt.T
    h=h.T
    zi=zi.T

    shape_out=temp.shape
    
    
#    if fields is not None:
#        for fld in fields:
#            expr=fld+'_remap=np.ma.masked_where(np.abs('+fld+'_remap-land_fill)<1.e-5,'+fld+'_remap)'
#            exec(expr)
#            expr=fld+'_remap='+fld+'_remap.T'
#            exec(expr)
#            expr=fld+'_remap=np.reshape('+fld+'_remap,(1,shape_out[0],shape_out[1],shape_out[2]))'
#            exec(expr)
            
#    temp=np.reshape_out(temp,(1,shape_out[0],shape_out[1],shape_out[2]))
#    salt=np.reshape_out(salt,(1,shape_out[0],shape_out[1],shape_out[2]))
#    h=np.reshape_out(h,(1,shape_out[0],shape_out[1],shape_out[2]))
#    zint=np.reshape_out(zi,(1,shape_out[0]+1,shape_out[1],shape_out[2]))            

        
    zout = 0.5*(zi + np.roll(zi,axis=1,shift=-1))
    zout = zout[:,0:-1,:]

    var_dict=dict.copy(self.var_dict[temp_name])
    var_dict['z_indices']=np.arange(0,len(R))
    var_dict['zb_indices']=np.arange(0,len(Rb))
    var_dict['zax_data']=R
    var_dict['zbax_data']=Rb
    var_dict['z_interfaces']=zi
    var_dict['dz']=h
    var_dict['Z']='Layer'
    var_dict['z']=zout
    var_dict['zunits']='kg m-3'
    var_dict['Ztype']='Isopycnal'


    self.add_field_from_array(temp,'temp_remap',var_dict=var_dict)
    self.add_field_from_array(salt,'salt_remap',var_dict=var_dict)            

#    if fields is not None:
#        for fld in fields:
#            expr='self.add_field_from_array('+fld+'_remap,\''+fld+'_remap\',var_dict=var_dict)'
#            print expr
#            exec(expr)


  def remap_vertical(self,fields=None,z_bounds=None,zbax_data=None,method='pcm',bndy_extrapolation=False):

    import pyale_mod
    """

    Re-mapping [fields] between generalized vertical coordinates:

    x1=self.var_dict['z_interfaces']  , and
    x2=z_bounds

    zbax_data is optional and is used only for plotting
    purposes.  For example, this could contain the coordinates
    of the underlying grid (e.g. sigma or isopycnal).

    Method is either,

    'pcm','plm','ppm',or 'pqm'

    following (White and Adcroft, JCP, 2008, vol 227, pp 7394-7422).

    
    
    """
    try:
        from midas import vertmap 
    except:
        print """ ALE/vertmap not installed """
        return 


    if z_bounds is None:
        print 'No output grid in call to vert_remap'
        return
    else:
        if z_bounds.ndim == 4:
            nx2=z_bounds.shape[1]-1
            ztype='Generalized'
        else:
            nx2=z_bounds.shape[0]-1
            ztype='Fixed'
    
    for fld in fields:

        # initialize an array for output data

        nt = vars(self)[fld].shape[0];nj=vars(self)[fld].shape[2];ni=vars(self)[fld].shape[3] 
        
        fld_out = np.ma.zeros((nt,nx2,nj,ni))

        vdict=self.var_dict[fld].copy()

        vdict['Z']='z_remap'

        vdict['zb_indices']=np.arange(0,nx2+1)

        if zbax_data is not None:
            vdict['zbax_data']=zbax_data
            zax_data=0.5*(zbax_data+np.roll(zbax_data,shift=-1))
            zax_data=zax_data[0:-1]
            vdict['zax_data']=zax_data            
        else:
            vdict['zbax_data']=None
            vdict['zax_data']=None
            
        vdict['z_indices']=np.arange(0,nx2+1)

        vdict['z_interfaces']=z_bounds.copy()


        if ztype == 'Fixed':
            vdict['dz']=np.zeros((nx2,nj,ni))            
            vdict['z']=np.zeros((nx2,nj,ni))
        else:
            vdict['z']=np.zeros((nx2,nj,ni))            
            vdict['dz']=np.zeros((nt,nx2,nj,ni))
            
        for n in np.arange(nt):

            # Convention is x is monotonic and increasing.
            # i.e. x(k+1)>x(k)
            # Midas convention is x(k)>x(k+1).
            # Output grid is further adjusted so that the
            # outer edges of x1,x2 are aligned. The modified
            # edges are stored in a temporary copy of the
            # global array t-slice here.  

        
            if self.var_dict[fld]['Ztype'] in ['Isopycnal','Generalized','Sigma']:
                xb1 = -self.var_dict[fld]['z_interfaces'][n,:]
                if ztype == 'Fixed':
                    xb2=-z_bounds
                else:
                    xb2=-np.take(z_bounds,[n],axis=0)                                
            else:
                xb1 = -self.var_dict[fld]['z_interfaces'][:]
                if ztype == 'Fixed':
                    xb2=-z_bounds.copy() # Force an array copy since we will be adjusting these
                                        # coordinates to match the outer edges of x1
                else:
                    xb2=-np.take(z_bounds,[n],axis=0)


            nx1=xb1.shape[0]-1
            xb2=sq(xb2)

            xb2[0,:,:]=xb1[0,:,:] # reset top interface to xb1[0]
            for k in np.arange(1,xb2.shape[0]-1):
                xb2[k,:,:]=np.maximum(xb2[k-1,:,:],xb2[k,:,:]) # avoid negative thicknesses 
                xb2[k,:,:]=np.minimum(xb2[k,:,:],xb1[nx1,:,:])

            xb2[nx2,:]=xb1[nx1,:]
            
            h2=np.roll(xb2,shift=-1,axis=0)-xb2
            h2=h2[:-1,:]

            


            if ztype == 'Fixed':
                vdict['z_interfaces'][:]=-xb2
                zout=xb2-np.roll(xb2,shift=-1,axis=0)
                vdict['z']=zout[:-1,:]
                vdict['dz'][:]=h2                
            else:
                vdict['z_interfaces'][n,:]=-xb2
                zout=xb2-np.roll(xb2,shift=-1,axis=0)
                vdict['z'][n,:]=zout[:-1,:]                
                vdict['dz'][n,:]=h2
                
            data=np.take(vars(self)[fld],[n],axis=0)
            data2=np.zeros((nx2,nj,ni))

            missing_value=1.e20
            data=np.ma.filled(data,missing_value)

            data=sq(data).T
            data2=data2.T
            xb1=xb1.T
            xb2=sq(xb2).T

            
            pyale_mod.pyale_mod.remap(data,data2,xb1,xb2,method,bndy_extrapolation=bndy_extrapolation,missing=missing_value)

            data2=data2.T
            xb2=xb2.T
            
            dz=np.roll(xb2,shift=-1,axis=0)-xb2
            dz=dz[:-1,:]
            mask=dz.copy()
            mask[mask>1.e-9]=1.0
            mask[mask<=1.e-9]=0.0

            data2=np.ma.masked_where(mask==0.0,data2)

            fld_out[n,:]=data2

        fnam=fld+'_remap'
        self.add_field_from_array(fld_out,fnam,var_dict=vdict)




    
  def adjust_thickness(self,field=None,min_thickness=0.0):
    """

    Adjust cell thicknesses based on grid.D 
    
    """

    if self.var_dict[field]['Z'] is None:
      return None
      

    if self.var_dict[field]['Ztype'] == 'Fixed':
        
        dz = self.var_dict[field]['dz']
        dz[dz<epsln]=epsln
        zb=np.cumsum(dz,axis=0)

        nz = vars(self)[field].shape[1]    
        D = np.tile(self.grid.D,(nz,1,1))
        floor_dz=D - np.roll(zb,shift=1,axis=0) # Distance from depth
        mask=np.logical_and(zb > D,np.roll(zb,shift=1,axis=0) <= D) # mask for bottom-most cell
        dz[mask]=floor_dz[mask]
        zb=np.cumsum(dz,axis=0)
        dz[zb>D]=epsln
        zb=np.cumsum(dz,axis=0)
        z=0.5*(zb+np.roll(zb,axis=0,shift=1))
        z[0,:]=0.5*dz[0,:]
        zb=np.concatenate((np.zeros((1,self.grid.jm,self.grid.im)),zb),axis=0)

        self.var_dict[field]['dz']=dz        
        if self.var_dict[field]['Zdir']==-1:
            self.var_dict[field]['z']=-z
            self.var_dict[field]['z_interfaces']=-zb
        else:
            self.var_dict[field]['z']=z
            self.var_dict[field]['z_interfaces']=zb
    else:
        
        dz = self.var_dict[field]['dz']

        zb=np.cumsum(dz,axis=1)
        
        nt= dz.shape[0]
        nz = vars(self)[field].shape[1]    
        D = np.tile(self.grid.D,(nt,nz,1,1))

        floor_dz=D - np.roll(zb,shift=1,axis=1) # Distance from depth
        mask=np.logical_and(zb > D,np.roll(zb,shift=1,axis=1) <= D) # mask for bottom-most cell
        dz[mask]=floor_dz[mask]
        dz= np.maximum(dz,min_thickness)
        zb=np.cumsum(dz,axis=1)
        z=0.5*(zb+np.roll(zb,axis=1,shift=1))
        z[:,0,:]=0.5*dz[:,0,:]
        zb=np.concatenate((np.zeros((nt,1,self.grid.jm,self.grid.im)),zb),axis=1)

        self.var_dict[field]['dz']=dz        
        if self.var_dict[field]['Zdir']==-1:
            self.var_dict[field]['z']=-z
            self.var_dict[field]['z_interfaces']=-zb
        else:
            self.var_dict[field]['z']=z
            self.var_dict[field]['z_interfaces']=zb                        
        
  def horiz_interp(self,field=None,target=None,src_modulo=False,method=1,PrevState=None,field_x=None,field_y=None,verbose=0):
    """
      Interpolate from a spherical grid to a general logically
      rectangular grid using a non-conservative \"bilinear\" interpolation
      algorithm (the default) or \"conservative\" area-weighted.
    """
    
    from midas import hinterp


    is_vector = False
    if field_x is not None:
        if field_y is not None:
            is_vector=True
        else:
            print """If field_x is not None then field_y should be not None as well"""
            return
    elif field_y is not None:
        print """If field_x is not None then field_y should be not None as well"""
        return

    if is_vector:
        if self.var_dict[field_x]['Z'] is not None:
            if self.var_dict[field_x]['Ztype'] is not 'Fixed':
                print """horiz_interp currently only configured for geopotential
              coordinate data """
                return None
    else:
        if self.var_dict[field]['Z'] is not None:
            if self.var_dict[field]['Ztype'] is not 'Fixed':
                print """horiz_interp currently only configured for geopotential
              coordinate data """
                return None


    
    open('input.nml','w+')
    
    deg_to_rad=np.pi/180.

    add_NP=True
    add_SP=True
    add_np=False
    add_sp=False


    nj_in = self.grid.lonh.shape[0]; ni_in = self.grid.lonh.shape[0]

    if method==0:
        if hasattr(self.grid,'x_T_bounds'):
            lon_in=self.grid.x_T_bounds
            lat_in=self.grid.y_T_bounds
        elif hasattr(self.grid,'lonq'):
            xax=self.grid.lonq
            yax=self.grid.latq
            lon_in,lat_in = np.meshgrid(xax,yax)
            lon_in=lon_in
            lat_in=lat_in
        else:
            print """ Unable to read grid cell boundaries on input grid"""
            return None


        if hasattr(target,'x_T'):
            nj=target.x_T.shape[0];ni=target.x_T.shape[1]
            lon_out = target.x_T_bounds
            lat_out = target.y_T_bounds
        elif hasattr(target,'x'):
            print """Conservative interpolation not available for
                  supergrids"""
            return None
        
    else:
        if hasattr(self.grid,'x_T'):
            lon_in=self.grid.x_T
            lat_in=self.grid.y_T
        elif hasattr(self.grid,'lonh'):
            xax1=self.grid.lonh
            yax1=self.grid.lath
            lon_in,lat_in = np.meshgrid(xax1,yax1)
        else:
            print """ Unable to read grid cell boundaries on input grid"""
            return None
            
        if hasattr(target,'x_T'):
            nj=target.x_T.shape[0];ni=target.x_T.shape[1]
            lon_out = target.x_T
            lat_out = target.y_T
        elif hasattr(target,'x'):
            nj=target.x.shape[0];ni=target.x.shape[1]
            lon_out = target.x
            lat_out = target.y

      
    max_lat_in = np.max(lat_in)
    if np.logical_and(max_lat_in < 90.0 - epsln,add_NP):
        add_np=True
        np_lat = np.reshape(np.tile([90.0],(ni_in)),(1,ni_in))
        np_lon = np.reshape(lon_in[-1,:],(1,ni_in))            
        lat_in=np.concatenate((lat_in,np_lat))
        lon_in=np.concatenate((lon_in,np_lon))
            
    min_lat_in = np.min(lat_in)

    if np.logical_and(min_lat_in > -90.0 + epsln,add_SP):
        add_sp=True
        sp_lat = np.reshape(np.tile([-90.0],(ni_in)),(1,ni_in))
        sp_lon = np.reshape(lon_in[0,:],(1,ni_in))                        
        lat_in=np.concatenate((sp_lat,lat_in))
        lon_in=np.concatenate((lon_in,sp_lon))

    
    ny=lat_in.shape[0];nx=lon_in.shape[1]
    
    lon_in=lon_in*deg_to_rad
    lat_in=lat_in*deg_to_rad    
    
    lon_out = lon_out*deg_to_rad
    lat_out = lat_out*deg_to_rad

    if hasattr(target,'wet'):
        mask_out=target.wet.copy()
    else:
        mask_out=np.ones((target.x_T.shape))

    
    
    if is_vector:
        if self.var_dict[field_x]['missing_value'] is not None:
            missing = np.float64(self.var_dict[field_x]['missing_value'])
        elif self.var_dict[field_x]['_FillValue'] is not None:
            missing = np.float64(self.var_dict[field_x]['_FillValue'])
        else:
            missing = np.float64(-1.e10)
            
        varin_x = vars(self)[field_x].astype('float64')
        varin_y = vars(self)[field_y].astype('float64')
        if np.ma.is_masked(varin_x):
            mask_in = np.ma.getmask(varin_x)
            mask=np.ones((mask_in.shape))
            mask[mask_in]=0
        else:
            mask=np.ones((varin_x.shape))

        nk=varin_x.shape[1];nt=varin_x.shape[0]
        varin_x[np.abs(varin_x-missing)<1.e-3*np.abs(varin_x)]=missing
        varin_y[np.abs(varin_y-missing)<1.e-3*np.abs(varin_y)]=missing        
    else:
        if self.var_dict[field]['missing_value'] is not None:
            missing = np.float64(self.var_dict[field]['missing_value'])
        elif self.var_dict[field]['_FillValue'] is not None:
            missing = np.float64(self.var_dict[field]['_FillValue'])
        else:
            missing = np.float64(-1.e10)
            
        varin = vars(self)[field].astype('float64')
        if np.ma.is_masked(varin):
            mask_in = np.ma.getmask(varin)
            mask=np.ones((mask_in.shape))
            mask[mask_in]=0
        else:
            mask=np.ones((varin.shape))
            
        varin[np.abs(varin-missing)<1.e-3*np.abs(varin)]=missing            

        nk=varin.shape[1];nt=varin.shape[0]
    
    
    if add_np:
      if is_vector:
          last_row=varin_x[:,:,-1,:]
          pole=np.ma.average(last_row,axis=2)
          pole=np.reshape(pole,(nt,nk,1,1))
          pole=np.tile(pole,(1,1,1,nx))
          varin_x=np.concatenate((varin_x,pole),axis=2)
          last_row=varin_y[:,:,-1,:]
          pole=np.ma.average(last_row,axis=2)
          pole=np.reshape(pole,(nt,nk,1,1))
          pole=np.tile(pole,(1,1,1,nx))
          varin_y=np.concatenate((varin_y,pole),axis=2)
          if hasattr(self.grid,'angle_dx'):
              angle_dx = self.grid.angle_dx
              angle_dx = angle_dx[np.newaxis,np.newaxis,:]
              angle_dx = np.tile(angle_dx,(nt,nk,1,1))              
              last_row=angle_dx[:,:,-1,:]
              pole=np.ma.average(last_row,axis=2)
              pole=np.reshape(pole,(nt,nk,1,1))
              pole=np.tile(pole,(1,1,1,nx))
              angle_dx=np.concatenate((angle_dx,pole),axis=2)              
      else:
          last_row=varin[:,:,-1,:]
          pole=np.ma.average(last_row,axis=2)
          pole=np.reshape(pole,(nt,nk,1,1))
          pole=np.tile(pole,(1,1,1,nx))
          varin=np.concatenate((varin,pole),axis=2)
      last_row=mask[:,:,-1,:]
      pole=np.ma.max(last_row,axis=2)
      pole=np.reshape(pole,(nt,nk,1,1))
      pole=np.tile(pole,(1,1,1,nx))
      mask=np.concatenate((mask,pole),axis=2)
    if add_sp:
      if is_vector:
          first_row=varin_x[:,:,0,:]
          pole=np.ma.average(first_row,axis=2)
          pole=np.reshape(pole,(nt,nk,1,1))
          pole=np.tile(pole,(1,1,1,nx))
          varin_x=np.concatenate((pole,varin_x),axis=2)
          first_row=varin_y[:,:,0,:]
          pole=np.ma.average(first_row,axis=2)
          pole=np.reshape(pole,(nt,nk,1,1))
          pole=np.tile(pole,(1,1,1,nx))
          varin_y=np.concatenate((pole,varin_y),axis=2)
          if hasattr(self.grid,'angle_dx'):
              first_row=angle_dx[:,:,0,:]
              pole=np.ma.average(first_row,axis=2)
              pole=np.reshape(pole,(nt,nk,1,1))
              pole=np.tile(pole,(1,1,1,nx))
              angle_dx=np.concatenate((pole,angle_dx),axis=2)              
      else:
          first_row=varin[:,:,0,:]
          pole=np.ma.average(first_row,axis=2)
          pole=np.reshape(pole,(nt,nk,1,1))
          pole=np.tile(pole,(1,1,1,nx))
          varin=np.concatenate((pole,varin),axis=2)
      first_row=mask[:,:,0,:]
      pole=np.ma.max(first_row,axis=2)
      pole=np.reshape(pole,(nt,nk,1,1))
      pole=np.tile(pole,(1,1,1,nx))
      mask=np.concatenate((pole,mask),axis=2)
      



    if src_modulo:
      mask=np.concatenate((mask,np.take(mask,[0],axis=3)),axis=3)
      mask=np.concatenate((np.take(mask,[-2],axis=3),mask),axis=3)                          
      if is_vector:
          varin_x=np.concatenate((varin_x,np.take(varin_x,[0],axis=3)),axis=3)
          varin_x=np.concatenate((np.take(varin_x,[-2],axis=3),varin_x),axis=3)                          
          varin_y=np.concatenate((varin_y,np.take(varin_y,[0],axis=3)),axis=3)
          varin_y=np.concatenate((np.take(varin_y,[-2],axis=3),varin_y),axis=3)
          if hasattr(self.grid,'angle_dx'):
              angle_dx=np.concatenate((angle_dx,np.take(angle_dx,[0],axis=3)),axis=3)
              angle_dx=np.concatenate((np.take(angle_dx,[-2],axis=3),angle_dx),axis=3)                                        
      else:
          varin=np.concatenate((varin,np.take(varin,[0],axis=3)),axis=3)
          varin=np.concatenate((np.take(varin,[-2],axis=3),varin),axis=3)                                    

      lat_in=np.concatenate((lat_in,np.take(lat_in,[0],axis=1)),axis=1)
      lat_in=np.concatenate((np.take(lat_in,[-2],axis=1),lat_in),axis=1)
      lon_in=np.concatenate((lon_in,np.take(lon_in,[0],axis=1)+2.0*np.pi),axis=1)
      lon_in=np.concatenate((np.take(lon_in,[-2],axis=1)-2.0*np.pi,lon_in),axis=1)            

    lon_in=np.float64(lon_in)
    lat_in=np.float64(lat_in)
    lon_out=np.float64(lon_out)
    lat_out=np.float64(lat_out)
      
    if is_vector:
        varin_x=np.ma.filled(varin_x,missing)
        varout_x=np.zeros((nt,nk,nj,ni))
        varin_y=np.ma.filled(varin_y,missing)
        varout_y=np.zeros((nt,nk,nj,ni))
        if hasattr(self.grid,'angle_dx'):
            x=varin_x.copy()
            y=varin_y.copy()
            varin_x = x*np.cos(angle_dx) + y*np.sin(angle_dx)
            varin_y = y*np.cos(angle_dx) - x*np.sin(angle_dx) 
#        hinterp.hinterp_mod.hinterp(lon_in.T,lat_in.T,mask.T,varin_x.T,lon_out.T,lat_out.T,mask_out.T,varout_x.T,False,method,missing)
#        hinterp.hinterp_mod.hinterp(lon_in.T,lat_in.T,mask.T,varin_y.T,lon_out.T,lat_out.T,mask_out.T,varout_y.T,False,method,missing)
        varin_x=np.float64(varin_x)
        varin_y=np.float64(varin_y)                
        
        
        hinterp.hinterp_mod.hinterp(lon_in.T,lat_in.T,varin_x.T,lon_out.T,lat_out.T,varout_x.T,False,method,missing)
        hinterp.hinterp_mod.hinterp(lon_in.T,lat_in.T,varin_y.T,lon_out.T,lat_out.T,varout_y.T,False,method,missing)        
        angle_dx = target.angle_dx[np.newaxis,np.newaxis,:]
        angle_dx = np.tile(angle_dx,(nt,nk,1,1))
        x_rot = varout_x*np.cos(angle_dx) - varout_y*np.sin(angle_dx)
        y_rot = varout_y*np.cos(angle_dx) + varout_x*np.sin(angle_dx)        
    else:
        varout=np.zeros((nt,nk,nj,ni))

        varin=np.float64(varin)
        
        hinterp.hinterp_mod.hinterp(lon_in.T,lat_in.T,varin.T,lon_out.T,lat_out.T,varout.T,False,method,missing,verbose)        
    
    if PrevState is not None:
      S=PrevState
    else:
      S = state(grid=target)


    if is_vector:
        vars(S)[field_x]=x_rot
        S.variables[field_x]=field_x
        vars(S)[field_y]=y_rot
        S.variables[field_y]=field_y
        var_dict_x=self.var_dict[field_x].copy()
        var_dict_y=self.var_dict[field_y].copy()
        var_dict_x['X']=str(var_dict_x['X']);var_dict_x['Y']=str(var_dict_x['Y'])
        var_dict_y['X']=str(var_dict_y['X']);var_dict_y['Y']=str(var_dict_y['Y'])
        var_dict_x['xax_data']=target.lonh
        var_dict_x['yax_data']=target.lath
        var_dict_y['xax_data']=target.lonh
        var_dict_y['yax_data']=target.lath        

        if var_dict_x['z_interfaces'] is not None:
            if var_dict_x['z_interfaces'].ndim == 4:
                zi=var_dict_x['z_interfaces'][0,:,0,0]
            elif var_dict_x['z_interfaces'].ndim == 3:
                zi=var_dict_x['z_interfaces'][:,0,0]
        
                
            zi=zi.reshape(len(zi),1,1)
            zi=np.tile(zi,(1,nj,ni))
            var_dict_x['z_interfaces']=zi
            var_dict_y['z_interfaces']=zi            

        if var_dict_x['Z'] is not None:        
            if var_dict_x['z'].ndim == 4:
                z=var_dict_x['z'][0,:,0,0]
            elif var_dict_x['z'].ndim == 3:
                z=var_dict_x['z'][:,0,0]
        
            z=z.reshape(len(z),1,1)
            z=np.tile(z,(1,nj,ni))
            var_dict['z']=z
            dz=var_dict_x['dz'][:,0,0]
            dz=dz.reshape(len(dz),1,1)
            dz=np.tile(dz,(1,nj,ni))
            var_dict_x['dz']=dz
            var_dict_y['dz']=dz                
    
    
    
        S.var_dict[field_x]=var_dict_x
        S.var_dict[field_y]=var_dict_y

        if S.var_dict[field_x]['masked']:
            vars(S)[field_x] = np.ma.masked_where(np.abs(varout_x - missing) < np.abs(missing)*1.e-3,x_rot)
            vars(S)[field_y] = np.ma.masked_where(np.abs(varout_y - missing) < np.abs(missing)*1.e-3,y_rot)            
  
    else:
        vars(S)[field]=varout
        S.variables[field]=field
        var_dict=self.var_dict[field].copy()
        var_dict['X']=str(var_dict['X']);var_dict['Y']=str(var_dict['Y'])
        var_dict['xax_data']=target.lonh
        var_dict['yax_data']=target.lath        

        if var_dict['z_interfaces'] is not None:
            if var_dict['z_interfaces'].ndim == 4:
                zi=var_dict['z_interfaces'][0,:,0,0]
            elif var_dict['z_interfaces'].ndim == 3:
                zi=var_dict['z_interfaces'][:,0,0]
        
                
            zi=zi.reshape(len(zi),1,1)
            zi=np.tile(zi,(1,nj,ni))
            var_dict['z_interfaces']=zi

        if var_dict['Z'] is not None:        
            if var_dict['z'].ndim == 4:
                z=var_dict['Z'][0,:,0,0]
                dz=var_dict['dz'][0,:,0,0]                
            elif var_dict['z'].ndim == 3:
                z=var_dict['z'][:,0,0]
                dz=var_dict['dz'][:,0,0]                                
            elif var_dict['z'].ndim == 1:
                z=var_dict['z'][:]
                try:
                    dz=var_dict['dz'][:]
                except:
                    dz=np.ones(z.shape)
                
            z=z.reshape(len(z),1,1)
            z=np.tile(z,(1,nj,ni))
            var_dict['z']=z

            dz=dz.reshape(len(dz),1,1)
            dz=np.tile(dz,(1,nj,ni))
            var_dict['dz']=dz    
    
    
    
        S.var_dict[field]=var_dict

        if S.var_dict[field]['masked']:
            vars(S)[field] = np.ma.masked_where(np.abs(varout - missing) < np.abs(missing)*1.e-3,varout)

            

    return S

  def grid_overlay(self,field=None,target=None,debug=0):

      """
      Use corner on target grid and centroids on source grid 
      to find the source grid indices which fall within the 
      target cell boundary defined by the corners.

      Calculate mean, max and min of source grid points on the current list.

      Fit the current list to a plane surface using optimize.leastsq method
      which minimizes

      J = np.sum(e(y,x)**2.0) = np.sum(z(y,x)-(a[0] + a[1]*x + a[2]*y)

      The grid roughness, std = J**0.5.

      Refer to horiz_interp or make_xgrid for conservative interpolation.
      
      X===Source grid centers
      O===Target grid corners
            
              n-w X     X    X     X     X       X n-e
                        
                       O----------------------O
                       |                      |
                  X    |X    X     X     X    |  X 
                       |                      |
                       |                      |
                       |                      |
                       O----------------------O
                  X     X    X     X     X       X
              s-w                                  s-e
      
                       """
      
      from scipy import optimize

      def linfit2d(x,y,z):
          """

      Fit a set of points in 2-d to a cartesian plane.

      z = a[0] + a[1]x + a[2]y

 """
          a_init = [0.0,0.0,0.0]

    
          fitfunc = lambda a, x, y: a[0] + a[1]*x + a[2]*y

          ff= fitfunc(a_init,x,y)
    
          errfunc = lambda a, x, y, z: z - fitfunc(a,x,y)

          err_res= errfunc(a_init,x,y,z)

          out = optimize.leastsq(errfunc,a_init,args=(x.flatten(),y.flatten(),z.flatten()),full_output=1)

          a_final=out[0]

          return a_final
      
      try:      
          from mpl_toolkits.basemap import interp as Interp
      except:
          print """ need mpl_toolkits.basemap for this function """
          return



      if self.grid.is_latlon :
          if target.is_cartesian is True:
              print """
                 grid_overlay does not work between cartesian grids and lat-lon grids
                    """
              return None

      if target.is_latlon :
          if self.grid.is_cartesian is True:
              print """
                 grid_overlay does not work between cartesian grids and lat-lon grids
                    """
              return None
          
      
      if self.grid.is_latlon:
          deg_to_rad=np.pi/180.
      else:
          deg_to_rad=1.0

      add_np=False
      add_sp=False    

      missing=-1.e20

      if field is not None:
          shape_out=vars(self)[field].shape
          if np.logical_or(shape_out[0]>1,shape_out[1]>1):
              print "grid_overlay is currently only written to handle lat-lon arrays without a time or vertical dimension"
              return None
          
      nj_in = self.grid.lath.shape[0]; ni_in = self.grid.lonh.shape[0]

      if hasattr(self.grid,'lonh'):
          lon_in=self.grid.lonh.copy()
          lat_in=self.grid.lath.copy()
          lon_edge_in=self.grid.lonq.copy()
          lat_edge_in=self.grid.latq.copy()          
      else:
          print """ Unable to read grid cell locations on input grid"""
          return None

      target_type = -1
      if hasattr(target,'x'):  # target grid is a super grid
          target_type=0
          nj=target.x.shape[0]-1;ni=target.x.shape[1]-1
          lon_out = target.x
          lat_out = target.y
      elif hasattr(target,'x_T_bounds'): # target grid is a model grid
          target_type=1
          nj=target.latq.shape[0]-1;ni=target.lonq.shape[0]-1
          lon_out = target.x_T_bounds.copy()
          lat_out = target.y_T_bounds.copy()

      if target_type < 0 :
          print 'unknown target grid type '
          raise
      
      if target.is_latlon is True:
          lon_out[lon_out<lon_edge_in[0]]=lon_out[lon_out<lon_edge_in[0]]+360.0
          lon_out[lon_out>lon_edge_in[-1]]=lon_out[lon_out>lon_edge_in[-1]]-360.0      

      i_indices = np.arange(0,ni_in) #.astype(int)
      i_indices = i_indices[np.newaxis,:]
      i_indices = np.tile(i_indices,(nj_in,1))
      j_indices = np.arange(0,nj_in) #.astype(int)
      j_indices = j_indices[:,np.newaxis]
      j_indices = np.tile(j_indices,(1,ni_in))

      x_out = Interp(i_indices,lon_in,lat_in,lon_out,lat_out)
      y_out = Interp(j_indices,lon_in,lat_in,lon_out,lat_out)      

      iwest = np.floor(x_out).astype(int)
      jsouth = np.floor(y_out).astype(int)      

      if field is None:
          return x_out, y_out
      else:
          data_in = np.ma.masked_invalid(sq(vars(self)[field]))
          ieast = np.roll(iwest,shift=-1,axis=1) # Left face of cell to the coordinate east
          jnorth = np.roll(jsouth,shift=-1,axis=0)
          
          meanval = np.ma.zeros((nj,ni))
          maxval  = np.ma.zeros((nj,ni))
          minval  = np.ma.zeros((nj,ni))
          std     = np.ma.zeros((nj,ni))
          count     = np.zeros((nj,ni)).astype(int)          

          
          
          for j in np.arange(0,nj):
              for i in np.arange(0,ni):
                  j_list=np.arange(jsouth[j,i],jnorth[j,i]+1)
                  if jnorth[j,i]<jsouth[j,i]:
                      j1 = jnorth[j,i];j2=jsouth[j,i]
                      j_list=np.arange(j1,j2+1)
                  i_list=np.arange(iwest[j,i],ieast[j,i]+1)
                  if debug == 1:
                      print 'j,i,iwest,ieast= ',j,i,iwest[j,i],ieast[j,i]
                  if ieast[j,i]<iwest[j,i]:
                      if debug == 1:
                          print 'ieast<iwest at j,i= ',j,i,ieast[j,i],iwest[j,i]
                      if target_type == 0:
                          if target.dict['cyclic_x'] and ieast[j,i]==0:
                              i_list=np.arange(iwest[j,i],ni_in)
                          else:
                              i1=ieast[j,i];i2=iwest[j,i]
                              i_list = np.arange(i1,i2+1)
#                              i_list = np.concatenate((np.arange(i2,ni_in),np.arange(0,i1+1)))                              
                      else:
                          if target.cyclic_x and ieast[j,i]==0:
                              i_list=np.arange(iwest[j,i],ni_in)
                          else:
                              i1=ieast[j,i];i2=iwest[j,i]
                              i_list = np.concatenate((np.arange(i2,ni_in),np.arange(0,i1+1)))

                      if debug == 1:
                          print 'i_list count=',len(i_list)
                          print 'i_list = ',i_list

                          
                  i_arr,j_arr = np.meshgrid(i_list,j_list)
                  if len(i_arr) > 1:
                      i_arr=i_arr.flatten()
                  if len(j_arr) > 1:
                      j_arr=j_arr.flatten()
                  b=data_in[j_arr,i_arr]
                  zp=np.zeros(b.shape)
                  x=np.linspace(0.,1.,len(i_list))
                  y=np.linspace(0.,1.,len(j_list))
                  xdata,ydata=np.meshgrid(x,y)

                  if np.prod(xdata.shape) > 2:
                      zcoef = linfit2d(xdata,ydata,b.reshape(xdata.shape))
                      zhat = zcoef[0]+zcoef[1]*xdata.flatten()+zcoef[2]*ydata.flatten()
                      zz = b.flatten()-zhat
                  else:
                      zz=np.zeros(b.shape)

                  if b.size > 0:
                      meanval[j,i] = b.mean()
                      maxval[j,i] = b.max()
                      minval[j,i] = b.min()
                      count[j,i]  = b.count()
                      std[j,i]    = zz.std()
                  else:
                      meanval[j,i] = missing
                      maxval[j,i] = missing
                      minval[j,i] = missing
                      count[j,i]  = 0
                      std[j,i]    = missing

          meanval=np.ma.masked_where(meanval==missing,meanval)
          maxval=np.ma.masked_where(maxval==missing,maxval)
          minval=np.ma.masked_where(minval==missing,minval)
          std=np.ma.masked_where(std==missing,std)  


          
          S = state(grid=target)
          var_dict=self.var_dict[field].copy()
          S.variables={}
          
          var_dict['_Fillvalue']=missing
          var_dict['Z'] = None
          var_dict['T'] = None          
          
          if hasattr(target,'lonh'):
              var_dict['xax_data']=target.lonh
              var_dict['yax_data']=target.lath
          else:
              xnew = target.grid_x.copy()
              xnew = 0.5*(xnew+np.roll(xnew,shift=-1))
              xnew = np.take(xnew,np.arange(0,ni))
              var_dict['xax_data']=xnew
              ynew = target.grid_y.copy()
              ynew = 0.5*(ynew+np.roll(ynew,shift=-1))
              ynew = np.take(ynew,np.arange(0,nj))              
              var_dict['yax_data']=ynew

          S.mean=meanval[np.newaxis,np.newaxis,:,:]
          S.var_dict['mean']=var_dict.copy()
          S.variables['mean']='mean'          
          S.max=maxval[np.newaxis,np.newaxis,:,:]
          S.var_dict['max']=var_dict.copy()
          S.variables['max']='max'          
          S.min=minval[np.newaxis,np.newaxis,:,:]
          S.var_dict['min']=var_dict.copy()
          S.variables['min']='min'                    
          S.count=count[np.newaxis,np.newaxis,:,:]
          S.var_dict['count']=var_dict.copy()
          S.variables['count']='count'                              
          S.std=std[np.newaxis,np.newaxis,:,:]
          S.var_dict['std']=var_dict.copy()
          S.variables['std']='std'                                        


          return S
      
  def pickle_it(self,file):
    
    pickle.dump(self,open(file,'wb'))

    return None



  def compress_field(self,field=None):

      if field is None:
          return None

      sout=None
      nt = vars(self)[field].shape[0]
      if self.var_dict[field]['masked']:
          for i in np.arange(0,nt):
              arr = np.ma.compressed(vars(self)[field][i,:])
              if sout is not None:
                  arr = arr[np.newaxis,:]
                  sout = np.concatenate((arr,sout),axis=0)
              else:
                  sout=arr
                  sout=sout[np.newaxis,:]
      else:
          return None
      
      return sout

  def uncompress_field(self,arr,field=None):

      
      if field is None:
          return None

      arr=np.ravel(arr)

      if self.var_dict[field]['masked']:
          sout = np.zeros(vars(self)[field].shape)
          mask=np.ma.getmask(vars(self)[field])
          sout[~mask]=arr
          vars(self)[field] = sout.copy()

      return None

  def eof(self,field=None,trunc=1.0):
# Calculate Eigenvectors and Eigenvalues
# of symmetric real covariance matrix of
# (field) and store this along with principal

      normalize = 1
      

      if field is None:
          return None

# Use field.mask to compress data

      arr=self.compress_field(field)
      nt=arr.shape[0];nv=arr.shape[1]

      if normalize == 1:
          v=np.max(np.var(arr,axis=0))
          arr=arr/v
          
      arr=arr-arr.mean(axis=0)
      


      if nt < nv:
# Convert eigenfunctions of time covariance matrix to
# eigenvectors of spatial covariance matrix

          cov=np.zeros((nt,nt))

# Compute upper part of time covariance matrix          
          for n in np.arange(0,nt):
              for m in np.arange(n,nt):
                  cov[m,n]=arr[n,:].dot(arr[m,:])/(nv-1.0)
              
          w,efunct=np.linalg.eigh(cov,UPLO='L')

          arg_sort = w.argsort()

          w=w[arg_sort]
          efunct=efunct[:,arg_sort]

              
          w=w[::-1]
          efunct=efunct[:,::-1]

          tv=np.sum(w)
          fcv=w/tv*1.e2

          rank=w.shape[0]

          cvv=0.0
          cutoff_percentage=trunc*100.0
          rank_cutoff=rank
          for n in np.arange(0,rank):
              cvv=cvv+fcv[n]
              print 'evec#=',n,' ; % ',fcv[n], ' cum % ',cvv
              if cvv > cutoff_percentage: 
                  rank_cutoff=n
                  break
              
          cv=np.cumsum(fcv)

          for n in np.arange(0,rank_cutoff):
              efunct[:,n]=efunct[:,n]/np.sqrt(w[n])

          for n in np.arange(1,rank_cutoff):
              print 'efunct.dot (0,',n,') = ', efunct[:,0].dot(efunct[:,n])
          

          efunc=np.zeros((nv,rank_cutoff))

          for j in np.arange(0,rank_cutoff):
              for i in np.arange(0,nv):
                  efunc[i,j]=arr[:,i].dot(efunct[:,j]) 

          for n in np.arange(1,rank_cutoff):
              print 'efunc.dot (0,',n,') = ', efunc[:,0].dot(efunc[:,n])
                  
      else:
          cov=np.zeros((nv,nv))
          for n in np.arange(0,nv):
              for m in np.arange(n,nv):
                  cov[m,n]=arr[:,n].dot(arr[:,m])/(nt-1.0)
              
          w,efunc=np.linalg.eigh(cov,UPLO='L')

          arg_sort = w.argsort()

          w=w[arg_sort]
          efunc=efunc[:,arg_sort]

              
          w=w[::-1]
          efunc=efunc[:,::-1]

          tv=np.sum(w)
          fcv=w/tv*1.e2

          rank=w.shape[0]

          cvv=0.0
          cutoff_percentage=trunc*100.0
          for n in np.arange(0,rank):
              cvv=cvv+fcv[n]
              print 'eigenvector=',n,' ; % ',fcv[n], ' cumulative % ',cvv
              if cvv > cutoff_percentage: 
                  rank_cutoff=n
                  break

          cv=np.cumsum(fcv)

          for n in np.arange(0,rank_cutoff):
              efunc[:,n]=efunc[:,n]/np.sqrt(w[n])

          for n in np.arange(1,rank_cutoff):
              print 'efunc.dot (0,',n,') = ', efunc[:,0].dot(efunc[:,n])

####

              
      for n in np.arange(0,rank_cutoff):
          norm=np.sqrt(efunc[:,n].dot(efunc[:,n]))
          if norm > epsln:
              rnorm=1.0/norm
              efunc[:,n]=efunc[:,n]*rnorm
          else:
              efunc[:,n]=0.0


      pc=np.zeros((nt,rank_cutoff))

      for n in np.arange(0,nt):
          for m in np.arange(0,rank_cutoff):
              pc[n,m]=arr[n,:].dot(efunc[:,m])

      for n in np.arange(0,rank_cutoff):
          norm = np.sqrt(pc[:,n].dot(pc[:,n])/(nt-1))
          rnorm = 1.0/max(norm,epsln)
          pc[:,n]=pc[:,n]*rnorm # normalized
          efunc[:,n]=efunc[:,n]*norm*v # data units
              
      for n in np.arange(1,rank_cutoff):
          print 'pc.dot (0,',n,') = ', pc[:,0].dot(pc[:,n])
              
      expression='self.'+field+'[0,0,:]'+'*0.0'
      nam=field+'_evec'

      var_dict=self.var_dict[field].copy()

      
      # eigenvectors of fields with a vertical
      # extent are not CURRENTLY calculated with
      # appropriate depth-weighting. [TO DO]
      var_dict['Z']='rank'
      var_dict['Ztype']='Fixed'
      var_dict['z']=np.arange(0,rank_cutoff)
      var_dict['dz']=None
      var_dict['z_interfaces']=None
      var_dict['zax_data']=np.arange(0,rank_cutoff)

      var_dict['T']=None



      
      self.create_field(expression,nam,var_dict=var_dict)

      vars(self)[nam]=np.tile(vars(self)[nam],(1,rank_cutoff,1,1))

      arr=self.compress_field(nam)

      self.uncompress_field(efunc[:,0:rank_cutoff].T,field=nam)

      cond=nam+'==0.0'
  
      self.mask_where(field=nam,condition=cond)

      pc_dict=self.var_dict[field].copy()
      pc_dict['Z']='Principal Component'
      pc_dict['Ztype']='Fixed'
      pc_dict['X']=None
      pc_dict['Y']=None
      pc_dict['z']=np.arange(0,rank_cutoff)
      pc_dict['dz']=None
      pc_dict['z_interfaces']=None
      pc_dict['zax_data']=np.arange(0,rank_cutoff)

      
      expression='self.'+field+'[:,0,0,0]*0.0'
      nam=field+'_pc'

      self.create_field(expression,nam,var_dict=pc_dict)

      vars(self)[nam]=np.tile(vars(self)[nam],(1,rank_cutoff,1,1))
      
      vars(self)[nam]=pc
              

              
              
  def sfc_buoyancy_production(self,sst=None, sss=None,heat_flux=None,fw_flux=None,salt_flux=None,p_ref=0.0,rho_bounds=None):

      import remap_sfc_fluxes
      
      rho0 = 1.e3
      Irho0 = 1.0/rho0
      cp=3989.0
      Icp=1.0/cp
      grav=9.8
      Igrav = 1.0/grav


      ny=vars(self)[sst].shape[2]
      
      SST=self.compress_field(sst)
      SSS=self.compress_field(sss)

      AREA=self.grid.Ah
      LAT=self.grid.y_T
      AREA=np.ma.masked_where(sq(vars(self)[sst][0,0,:].mask),AREA)
      LAT=np.ma.masked_where(sq(vars(self)[sst][0,0,:].mask),LAT)
      AREA = np.ma.compressed(AREA)
      LAT = np.ma.compressed(LAT)      
    
      nt=SST.shape[0];nv=SST.shape[1]
      p_ref=np.reshape(p_ref,(1,1))
      P_REF=np.tile(p_ref,(nt,nv))
      rho = wright_eos(SST,SSS,P_REF)    
      drho_dT = alpha_wright_eos(SST,SSS,P_REF)
      drho_dT = -Irho0*drho_dT
      drho_dS = Irho0*beta_wright_eos(SST,SSS,P_REF)
      HFLX=self.compress_field(heat_flux)
      bf_heat = drho_dT*HFLX*Icp*grav
      SALT_FLUX=self.compress_field(salt_flux)
      FW_FLUX = self.compress_field(fw_flux)
      bf_fw = -grav*drho_dS*(SALT_FLUX*1.e3 - FW_FLUX*SSS)

      nrhob = rho_bounds.shape[0]
      nrho=nrhob-1

      dR=np.zeros((nrho))
      IdR=np.zeros((nrho))
      R=np.zeros((nrho))
      
      for k in np.arange(0,nrho):
          dR[k]=rho_bounds[k+1]-rho_bounds[k]
          R[k]=rho_bounds[k]+0.5*dR[k]
          IdR[k]=1.0/dR[k]
          
      bf_heat_d = np.zeros((nt,nrho,ny))
      bf_fw_d = np.zeros((nt,nrho,ny))
      tmp_arr = np.zeros((nt,nrho))
      tmp_arr2 = np.zeros((nt,nrho))            

      area_mask=np.zeros(LAT.shape)

      for j_indx in np.arange(0,ny):
          area_mask=AREA.copy()
          area_mask[LAT<self.grid.latq[j_indx]]=0.0          
          for n in np.arange(0,nt):
              remap_sfc_fluxes.remap_sfc_fluxes.remap(bf_heat[n,:]*area_mask[:],rho[n,:],rho_bounds,tmp_arr[n,:])
              remap_sfc_fluxes.remap_sfc_fluxes.remap(bf_fw[n,:]*area_mask[:],rho[n,:],rho_bounds,tmp_arr2[n,:])
          bf_heat_d[:,:,j_indx]=tmp_arr[:,:].copy()
          bf_fw_d[:,:,j_indx]=tmp_arr2[:,:].copy()  

      for k in np.arange(0,nrho):
          bf_heat_d[:,k,:]=bf_heat_d[:,k,:]*Igrav*IdR[k]
          bf_fw_d[:,k,:]=bf_fw_d[:,k,:]*Igrav*IdR[k]

      vdict=self.var_dict[heat_flux].copy()
      vdict['units']='kg m-3 s-1'
      self.add_field_from_array(vars(self)[heat_flux],'thermal_buoyancy_flux',var_dict=vdict)
      self.uncompress_field(bf_heat,field='thermal_buoyancy_flux')
      self.add_field_from_array(vars(self)[heat_flux],'fw_buoyancy_flux',var_dict=vdict)
      self.uncompress_field(bf_fw,field='fw_buoyancy_flux')
      vdict['units']='m3 s-1'      
      vdict['X']=None
      vdict['Z']='potential_density'
      vdict['zax_data']=R
      vdict['zbax_data']=rho_bounds
      vdict['zunits']='kg m-3'
      vdict['Zdir']=-1
      vdict['Ztype']='Fixed'
      self.add_field_from_array(bf_heat_d,'thermal_buoyancy_flux_d',var_dict=vdict)
      self.add_field_from_array(bf_fw_d,'fw_buoyancy_flux_d',var_dict=vdict)            
      
      return None


  def write_nc(self,filename=None,fields=None,format='NETCDF3_CLASSIC',append=False):

    import os.path      
    """

    Write (fields) to an NC file. 

    """  
    if fields is None:
      return None

    if filename is None:
      return None

    file_exists = False
    dims=[]
    vars=[]

    
    if append == True:
        if os.path.exists(filename):
            print '%s exists, appending ' %(filename)
            file_exists = True
            f=nc.Dataset(filename,'a',format=format)
            dimlist = f.dimensions
            for d in dimlist:
                dims.append(str(d))
            varlist = f.variables
            for v in varlist:
                vars.append(str(v))
        else:
            f=nc.Dataset(filename,'w',format=format)            
    else:
        f=nc.Dataset(filename,'w',format=format)
        
    write_interfaces=False
    
    if file_exists == True:
        outv=[]
        for field in fields:
            if field not in vars:
                print " Attempting to append %(fld)s to existing file %(fn)s "%{'fld':field,'fn':filename}
                return

            dim_nam=str(self.var_dict[field]['T'])
            if string.lower(dim_nam).count('none') == 0:
                tdat_=self.var_dict[field]['tax_data']
                if np.isscalar(tdat_):
                    tdat=[]
                    tdat.append(tdat_)
                else:
                    tdat=tdat_

                tv  = f.variables[self.var_dict[field]['T']]
                tunits_out = tv.units
                tunits=self.var_dict[field]['tunits']

            tdat_=[]
            if tunits != tunits_out:
               dates_in = self.var_dict[field]['dates']
               t_delta=num2date(0,tunits,'standard')-num2date(0.,tunits_out,'standard')
#               print """ Time units are not identical in write_nc.append """
#               print 't_delta (days) tunits - tunits_out = ',t_delta.days
#               print 'original timestamp = ',tdat
               for t in tdat:
                   if tunits_out[0:3]=='min':
                       tdat_.append(t + t_delta.days * 1440.)
                   if tunits_out[0:3]=='sec':
                       tdat_.append(t + t_delta.days * 86400.)
                   if tunits_out[0:3]=='day':
                       tdat_.append(t + t_delta.days)
                   if tunits_out[0:3]=='hou':
                       tdat_.append(t + t_delta.days * 24.)                                                                     
                   
               tdat=tdat_

#               print "adjusted time stamps = ",tdat
            outv.append(f.variables[field])
            
        nt=len(tdat)
        tstart=len(tv[:])
    else:    
        for field in fields:
            dim_nam=str(self.var_dict[field]['T'])
            if dim_nam not in dims and string.lower(dim_nam).count('none') == 0:
                dims.append(dim_nam)
                tdat_=self.var_dict[field]['tax_data']
                if np.isscalar(tdat_):
                    tdat=[]
                    tdat.append(tdat_)
                else:
                    tdat=tdat_
                tdim=f.createDimension(dim_nam,None)
                tv=f.createVariable(dim_nam,'f8',(dim_nam,))
                if self.var_dict[field].has_key('tunits'):
                    tv.units= self.var_dict[field]['tunits']
                if self.var_dict[field].has_key('calendar'):              
                    tv.calendar= self.var_dict[field]['calendar']
                nt=len(tdat)
                tstart = 0
                tv.cartesian_axis = 'T'          
            dim_nam=str(self.var_dict[field]['Z'])
            if dim_nam not in dims and string.lower(dim_nam).count('none') == 0:
                dims.append(dim_nam)
                xdat=self.var_dict[field]['zax_data'][:]
                xdim=f.createDimension(dim_nam,len(xdat))
                xv=f.createVariable(dim_nam,'f8',(dim_nam,))
                xv[:]=xdat
                xv.units =   self.var_dict[field]['zunits']
                xv.direction = self.var_dict[field]['Zdir']
                xv.cartesian_axis = 'Z'
                if self.var_dict[field]['Ztype'] in ['Generalized','Isopycnal'] and write_interfaces is False:
                    if 'z_interfaces'  in self.var_dict[field]:
                        if self.var_dict[field]['z_interfaces'] is not None:
                            write_interfaces = True
                            ifield=field
                            zi=self.var_dict[field]['z_interfaces'][:]
                            ziax=self.var_dict[field]['zbax_data'][:]
            dim_nam=str(self.var_dict[field]['Y'])
            if dim_nam not in dims and string.lower(dim_nam).count('none') == 0:
                dims.append(dim_nam)
                xdat=self.var_dict[field]['yax_data'][:]
                xdim=f.createDimension(dim_nam,len(xdat))
                xv=f.createVariable(dim_nam,'f8',(dim_nam,))
                xv.units =   self.var_dict[field]['yunits']
                xv[:]=xdat
                xv.cartesian_axis='Y'
            dim_nam=str(self.var_dict[field]['X'])
            if dim_nam not in dims and string.lower(dim_nam).count('none') == 0:
                dims.append(dim_nam)
                xdat=self.var_dict[field]['xax_data'][:]
                xdim=f.createDimension(dim_nam,len(xdat))
                xv=f.createVariable(dim_nam,'f8',(dim_nam,))
                xv.units =   self.var_dict[field]['xunits']
                xv[:]=xdat
                xv.cartesian_axis='X'
            vars.append(field)
    
        outv=[]
        n=0

        if write_interfaces:
            dim_nam='interfaces'
            dims.append(dim_nam)
            xdim=f.createDimension(dim_nam,len(ziax))
            xv=f.createVariable(dim_nam,'f8',(dim_nam,))
            xv[:]=ziax

      
        for field in fields:
            dims=[]
            if self.var_dict[field]['T'] is not None:
                dims.append(str(self.var_dict[field]['T']))
            if self.var_dict[field]['Z'] is not None:
                dims.append(str(self.var_dict[field]['Z']))
            if self.var_dict[field]['Y'] is not None:
                dims.append(str(self.var_dict[field]['Y']))
            if self.var_dict[field]['X'] is not None:
                dims.append(str(self.var_dict[field]['X']))

            if DEBUG == 1:
                print 'field=',field,'dims= ',dims

            if self.var_dict[field]['_FillValue'] is not None:
                FillValue = self.var_dict[field]['_FillValue']          
                var=f.createVariable(field,'f4',dimensions=dims,fill_value=FillValue)
            else:
                var=f.createVariable(field,'f4',dimensions=dims)      
            if self.var_dict[field]['missing_value'] is not None:
                var.missing_value = self.var_dict[field]['missing_value']

            if 'longname' in self.var_dict[field].keys():
                var.longname = self.var_dict[field]['longname']
          
            if 'units' in self.var_dict[field].keys():
                var.units = self.var_dict[field]['units']
          
            outv.append(var)


        if write_interfaces:
            dims=[]
            if self.var_dict[ifield]['T'] is not None:
                dims.append(str(self.var_dict[ifield]['T']))
            if self.var_dict[ifield]['Z'] is not None:
                dims.append('interfaces')        
            if self.var_dict[ifield]['Y'] is not None:
                dims.append(str(self.var_dict[ifield]['Y']))
            if self.var_dict[ifield]['X'] is not None:
                dims.append(str(self.var_dict[ifield]['X']))

            var=f.createVariable('eta','f4',dims)
            outv.append(var)
      

        m=0
        for field in fields:
            if self.var_dict[field]['T'] is  None:
                outv[m][:]=sq(self.__dict__[field][:])
            m=m+1

    m=0;p=0
    for field in fields:
        if self.var_dict[field]['T'] is not None:
            for n in np.arange(tstart,nt+tstart):
                if self.var_dict[field]['X'] is None and self.var_dict[field]['Y'] is None and self.var_dict[field]['Z'] is None:
                    outv[m][n]=self.__dict__[field][n-tstart]
                else:
                    outv[m][n,:]=sq(self.__dict__[field][n-tstart,:])

                if write_interfaces and m == 0:
                    outv[-1][n,:]=sq(zi[n-tstart,:])
        

                tv[n]=tdat[n-tstart]
                
            m=m+1
        
    
    f.sync()
    f.close()


  def calculate_bias(self,path=None,varin=None,varout=None,monthly_clim=False,ann_clim=False):
      
      """
      Read a field from path/varin and calculate bias statistics
      with respect to varout

      """

      if path is None or varin is None or varout is None:
          print """ path, varin and varout need to be specified """
          return

      grid_obs = rectgrid(path,var=varin)
      O=state(path,grid=grid_obs,fields=[varin],default_calendar='noleap')

      
      if ~O.var_dict[varin].has_key('Ztype'):
          O.var_dict[varin]['Ztype'] = 'Fixed'

      if monthly_clim:
          O.monthly_avg(varin,vol_weight=True)
          var = varin+'_monthly'
          O.del_field(varin)
          O.rename_field(var,varin)
          Obs=O.horiz_interp(varin,target=self.grid)
      elif ann_clim:
          O.time_avg(varin,vol_weight=True)
          var = varin+'_tav'
          O.del_field(varin)
          O.rename_field(var,varin)          
          Obs=O.horiz_interp(varin,target=self.grid,src_modulo=True)          
      else:
          Obs=O.horiz_interp(varin,target=self.grid,src_modulo=True)
          
      self.obs = {}
      self.obs[varin]=Obs



      return


  
  def add_figure(self,name,commands):
    """
    
    Add (commands) to write figures from self

    """

    try:
      self.fig_dict[name] = commands
    except:
      
      self.fig_dict={}
      self.fig_dict[name] = commands
      

      
      
  def __getstate__(self):

    dict = self.__dict__.copy()

    try:
      del dict['rootgrp']
    except:
      pass

    for v in self.variables:
      try:
        del dict['var_dict'][v]['rootgrp']
      except:
        pass
    del dict['variables']
      
    return dict

  def __setstate__(self,dict):


#    if dict['is_MFpath'] :
#      dict['rootgrp']=nc.MFDataset(dict['path'])
#    else:
#      dict['rootgrp']=nc.Dataset(dict['path'])
      
    dict['variables'] = {}
    
    for v in dict['var_dict'].keys():
      dict['variables'][v] = v
      
    self.__dict__.update(dict)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
