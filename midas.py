#==============================
#
#  This Python module contains a set of
# routines for handling data on a 2-dimensional
# lattice of grid cells lying in a
# cartesian space or on a sphere.  The vertical
# coordinate accomodates generailzed layered data
# z(x,y,k,t) or fixed surfaces z(k). 
#
# MIDAS (Modular Isopycnal Data Analysis System)
# was first developed by Matthew Harrison
# in 2011-2012 while working in the GFDL Oceans and Climate
# Group.
#
#===============================

import numpy as np
import netCDF4 as nc
from netCDF4 import num2date
import string
import copy
from datetime import *
from mpl_toolkits.basemap import shiftgrid,cm
import types
import matplotlib.pyplot as plt
try:
    import matplotlib.animation as animation
    HAVE_MATPLOTLIB_ANIMATION=True
except:
    HAVE_MATPLOTLIB_ANIMATION=False
    pass
import pickle
import scipy as sp
from midas_grid_gen import *

epsln = 1.e-2
R_earth = 6371.e3
DEBUG = 1
Omega=7.295e-5

def sq(arr):
    """
    Shorthand for np.squeeze()
    """
    res=np.squeeze(arr)

    return res
  
#==============================
# Working with Images
#==============================

def image_from_array(arr):

  from PIL import Image

  arr_min=np.min(arr)
  arr=arr+arr_min
  arr_max=np.max(arr)
  arr=255*arr/arr_max
  arr=arr.astype('uint8')

  n=arr.shape[0];m=arr.shape[1]
  img=Image.fromarray(arr)

  return img


def array_from_image(img,flipud=False):

  from PIL import Image
  from scipy.misc import fromimage

  if flipud:
    img=img.transpose(Image.FLIP_TOP_BOTTOM)
    
  arr=fromimage(img)

  return arr

#==============================
# Modify some image attributes and
# append to a pre-existing list
#==============================

def store_frame(im,fig,ims):

  def setvisible(self,vis):
    for c in self.collections: c.set_visible(vis)


  im.set_visible = types.MethodType(setvisible,im,None)
  im.axes = plt.gca()
  im.figure = fig
  ims.append([im])


def unpickle(file):

#==============================
# Unpickle a state
#==============================

  S=pickle.load(open(file,'rb'))
  
  return S


#==============================
# Equation of State 
#==============================

def wright_eos(T,S,p):
  """
  
 **********************************************************************
   The subroutines in this file implement the equation of state for   *
   sea water using the formulae given by  Wright, 1997, J. Atmos.     *
   Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
 ***********************************************************************
    Converted to Python from F90 by M Harrison 10/11.

 Calculate seawater equation of state, given T[degC],S[PSU],p[Pa]
 Returns density [kg m-3]
 
 ***********************************************************************
 
 """

  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  lam = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  I_denom = 1.0 / (lam + al0*(p+p0))  
  rho = (p + p0) * I_denom

  return rho 

def alpha_wright_eos(T,S,p):
  """

**********************************************************************
   The subroutines in this file implement the equation of state for   *
   sea water using the formulae given by  Wright, 1997, J. Atmos.     *
   Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
 ***********************************************************************
   Converted to Python from F90 by M Harrison 10/11.

 Calculate seawater equation of state, given T[degC],S[PSU],p[Pa]
 Returns partial derivative of density with respect to temperature [kg m-3 C-1]

 ***********************************************************************
 
 """
  
  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  lam = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  I_denom = 1.0 / (lam + al0*(p+p0))  
  I_denom2 = I_denom*I_denom
  drho_dT =  I_denom2*(lam*(b1+T*(2*b2 + 3*b3*T) + b5*S) - (p+p0)*((p+p0)*a1 + (c1+T*(2*c2 + 3*c3*T) + c5*S)))

  return drho_dT

def beta_wright_eos(T,S,p):
  """

 **********************************************************************
   The subroutines in this file implement the equation of state for   *
   sea water using the formulae given by  Wright, 1997, J. Atmos.     *
   Ocean. Tech., 14, 735-740.  Coded by R. Hallberg, 7/00.            *
 ***********************************************************************
   Converted to Python from F90 by M Harrison 10/11.

 Calculate seawater equation of state, given T[degC],S[PSU],p[Pa]
 Returns partial derivative of density with respect to salinity [kg m-3 PSU-1]
 
 ***********************************************************************
 
 """
  
  a0 = 7.057924e-4; a1 = 3.480336e-7; a2 = -1.112733e-7;
  b0 = 5.790749e8;  b1 = 3.516535e6;  b2 = -4.002714e4;
  b3 = 2.084372e2;  b4 = 5.944068e5;  b5 = -9.643486e3;
  c0 = 1.704853e5;  c1 = 7.904722e2;  c2 = -7.984422;
  c3 = 5.140652e-2; c4 = -2.302158e2; c5 = -3.079464;

  al0 = a0 + a1*T +a2*S
  p0  = b0 + b4*S + T * (b1 + T*(b2 + b3*T) + b5*S)
  lam = c0 +c4*S + T * (c1 + T*(c2 + c3*T) + c5*S)
  I_denom = 1.0 / (lam + al0*(p+p0))  
  I_denom2 = I_denom*I_denom
  drho_dS =  I_denom2*(lam*(b4+b5*T) - (p+p0)*((p+p0)*a2 + (c4+c5*T)))

  return drho_dS


def get_axis_cart(dimension,dimname=None):
  """

 **********************************************************************
 Determine axis orientation using a restricted list of detectable
 parameters, returns [X,Y,Z,T]
 **********************************************************************

 
 """

  valid_x_units = ['cm','m','meters','km','degrees_east','degrees_e','degree_e','deg_e']
  valid_y_units = ['cm','m','meters','km','degrees_north','degrees_n','degree_m','deg_n']
  valid_z_units = ['cm','m','meters','km','interface','layer']
  valid_t_units = ['seconds','minutes','hours','days','months','years']

  cart = None

  # first try to retreive the cartesian attribute

  try: 
    cart  = getattr(dimension,'cartesian_axis')
  except:
    pass

  try: 
    cart  = getattr(dimension,'axis')
  except:
    pass

  if cart is None:
    try:
      ax_units = getattr(dimension,'units')
      for units in valid_t_units:
        if string.lower(ax_units).count(units) > 0:
          cart = 'T'
      for units in valid_y_units:
        if string.lower(ax_units).count(units) > 0:
          cart = 'Y'
      for units in valid_x_units:
        if string.lower(ax_units).count(units) > 0:
          cart = 'X'
      for units in valid_z_units:
        if string.lower(ax_units).count(units) > 0:
          cart = 'Z'          
    except:
       pass

  
  if cart is None and dimname is not None:
      if string.lower(dimname).count('latitude') >0:
          cart = 'Y'
      if string.lower(dimname).count('longitude') >0:          
          cart = 'X'
          
      
  if cart == 'Z':
    try:
      orient = getattr(dimension,'positive')
      if string.lower(orient) == 'down':
        orientation = -1
      elif orient == -1:
        orientation = -1
    except:
      pass


  return cart

def get_axis_direction(dimension):

  dir = 1

  try:
    orient = getattr(dimension,'positive')
    if string.lower(orient) == 'down':
      dir = -1
    elif orient == -1:
      dir = -1
  except:
    pass

  return dir

def min_resolution(grid=None,x=None,y=None):
  if grid is not None:
    dx = grid.lonh-np.roll(grid.lonh,1)
    dy = grid.lath-np.roll(grid.lath,1)
    return np.min(dx[1:]),np.min(dy[1:])
  elif np.logical_or(x is None,y is None):
    if x is not None:
      dx = x-np.roll(x,1)
      return np.min(dx[1:]),None
    else:
      dy = y-np.roll(y,1)
      return None,np.min(dy[1:])
  else:
    print """
      Either x or y or grid must exist in call to min_resolution"""
    raise

def max_resolution(grid=None,x=None,y=None):
  if grid is not None:
    dx = grid.lonh-np.roll(grid.lonh,1)
    dy = grid.lath-np.roll(grid.lath,1)
    return np.max(dx[1:]),np.max(dy[1:])
  elif np.logical_or(x is None,y is None):
    if x is not None:
      dx = x-np.roll(x,1)
      return np.max(dx[1:]),None
    else:
      dy = y-np.roll(y,1)
      return None,np.max(dy[1:])
  else:
    print """
      Either x or y or grid must exist in call to min_resolution"""
    raise
    
  
def find_geo_bounds(grid,x=None,y=None):
  [min_dx,min_dy] = min_resolution(grid)
  [max_dx,max_dy] = max_resolution(grid)

  xs=None;xe=None
  if x is not None:
    xmin=x[0];xmax=x[1]
    xs = np.nonzero(np.abs(grid.lonh - xmin) <= max_dx)[0][0]
    xe = np.nonzero(np.abs(grid.lonh - xmax) <= max_dx)[0][0]

  ys=None;ye=None
  if y is not None:
    ymin=y[0];ymax=y[1]
    ys = np.nonzero(np.abs(grid.lath - ymin) <= max_dy)[0][0]
    ye = np.nonzero(np.abs(grid.lath - ymax) <= max_dy)[0][0]

  return xs,xe,ys,ye

def find_axis_bounds(axis,x=None):
  [min_dx,junk] = min_resolution(x=axis)
  [max_dx,junk] = max_resolution(x=axis)  

  xs=None;xe=None
  if x is not None:
    xmin=x[0];xmax=x[1]
    xs = np.nonzero(np.abs(axis - xmin) <= max_dx)[0][0]
    xe = np.nonzero(np.abs(axis - xmax) <= max_dx)[0][0]

  return xs,xe

def find_date_bounds(dates_in,tmin,tmax):
  if type(dates_in[0]) is not datetime:
    dates = []
    for i in np.arange(0,len(dates_in)):
      mon=int(dates_in[i].strftime()[5:7])
      year=np.maximum(int(dates_in[i].strftime()[0:4]),1)
      day=int(dates_in[i].strftime()[8:10])
      date=datetime(year,mon,day)
      dates.append(date)

  else:
    dates = dates_in


  
  ts=-1;te=-1
  for i in np.arange(0,np.maximum(1,len(dates)-1)):
    if ts == -1 and dates[i] >= tmin:
      ts = i
    if te == -1 and dates[i+1] > tmax:
      te = i
      break

  if te == -1:
    te=len(dates_in)-1
    
  return ts,te

def get_months(dates_in):
  if type(dates_in[0]) is not datetime:
    months = []
    for i in np.arange(0,len(dates_in)):
      mon=int(dates_in[i].strftime()[5:7])
      months.append(mon)
  else:
    months = []
    for i in np.arange(0,len(dates_in)):
      mon.append=dates_in[i].month

  return months


def make_monthly_axis(year=1900):
    dates=[];delta=[]
    for i in np.arange(1,13):
        dates.append(datetime(year,i,1))

    dates.append(datetime(year+1,1,1))
    
    for i in np.arange(0,12):
        delta.append((dates[i+1]-dates[i])/2)

    for i in np.arange(0,12):
        dates[i]=dates[i]+delta[i]

    return dates[0:12]
        

   
class generic_grid(object):

  """A grid object is for the horizontal lattice of points lying on
  the sphere. A generic_grid class is constructed by reading the
  lat-lon vectors of points associated with (var) from file (path).
  In order to construct an array of cells, an additional lattice
  of points are required in order to define the primary cell perimeters.
  The primary lattice defines the (T) cell centers, and the perimeter
  lattice of (Q) points are located at coordinates (y_T,x_T) and (y_T_bounds,x_T_bounds)
  respectively.  There are (jm,im) T cell points and (jm+1,im+1) perimeter
  locations. The cell face lengths and cell areas are constructed using
  lines of constant latitude and longitude.  

     
           +-------+-------+-------+
           :       :       :       :
           :       :       :       :
           :       :       :       :           
           +-------+-------Q-------:< Y_T_bounds
           :       :       :       :
           :       :   T   :       :
           :       :       :       :           
           +-------Q-------+-------+< Y_T_bounds
                   ^       ^
                  X_T_bounds
  """
  
  def __init__(self,path=None,cyclic=False,tripolar_n=False,var=None,simple_grid=False,supergrid=None,refine=1,lon=None,lat=None):



    if supergrid is not None:
        var_dict = {}
        x=supergrid.x; y=supergrid.y
        grid_x=supergrid.grid_x; grid_y=supergrid.grid_y
        self.lonh=grid_x[1::refine+1]
        self.lath=grid_y[1::refine+1]
        self.lonq=grid_x[0::refine+1]
        self.latq=grid_y[0::refine+1]
        self.x_T=x[1::refine+1,1::refine+1]
        self.y_T=y[1::refine+1,1::refine+1]

        self.x_T_bounds,self.y_T_bounds = np.meshgrid(self.lonq,self.latq)
        return

    if lon is not None and lat is not None:
        lon_range=lon.max()-lon.min()
        delta_lon = lon_range/lon.shape[1]
        lat_range=lat.max()-lat.min()
        delta_lat = lat_range/lat.shape[0]        
        self.lonh=np.arange(lon.min(),lon.min()+lon_range,delta_lon)
        self.lath=np.arange(lat.min(),lat.max()+lat_range,delta_lat)
        self.x_T=lon.copy()
        self.y_T=lat.copy()
    if lon is None and lat is None:
        f=nc.Dataset(path)
    else:
        f=None
    
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
              self.lonh = sq(f.variables[var_dict['X']][:])
          elif lon is not None:
              self.lonh = sq(lon[0,:])
          else:
              print "Longitude axis not detected "
              raise              
          if var_dict['Y'] is not None and lat is None:
              lat_axis = f.variables[var_dict['Y']]
              self.lath = sq(f.variables[var_dict['Y']][:])
          elif lat is not None:
              self.lath=sq(lat[:,0])
          else:
              print "Latitude axis not detected "
              raise
      


        

      self.lonq = 0.5*(self.lonh + np.roll(self.lonh,-1))

      self.lonq[-1] = 2.0*self.lonq[-2] -self.lonq[-3]
      lon_end = self.lonq[-1]+2.0*(self.lonh[-1]-self.lonq[-1])
      
      self.lonq=np.hstack((self.lonq,[lon_end]))
      self.latq = 0.5*(self.lath + np.roll(self.lath,-1))

      self.latq[-1] = 2.0*self.latq[-2]-self.latq[-3]
      lat_end=self.latq[-1]+2.0*(self.lath[-1]-self.latq[-1])
      self.latq=np.concatenate((self.latq,[lat_end]))      
      
      self.im = len(self.lonh)
      self.jm = len(self.lath)
      
      if simple_grid is True:
          self.simple_grid = True
          return
      else:
          self.simple_grid = False

      self.x_T,self.y_T = np.meshgrid(self.lonh,self.lath)
      
      xtb=0.5*(self.x_T + np.roll(self.x_T,shift=1,axis=1))
      xtb0=2.0*xtb[:,-1]-xtb[:,-2]
      xtb0=xtb0[:,np.newaxis]
      xbt=np.hstack((xtb,xtb0))
      xtb0=2.0*xtb[:,1]-xtb[:,2]
      xtb[:,0]=xtb0
      ytb=0.5*(self.y_T + np.roll(self.y_T,shift=1,axis=1))
      ytb0=2.0*ytb[-1,:]-ytb[-2,:]
      ytb0=ytb0[np.newaxis,:]
      ybt=np.vstack((ytb,ytb0))
      ytb0=2.0*ytb[1,:]-ytb[2,:]
      ytb[0,:]=ytb0
      self.x_T_bounds=xtb.copy()
      self.y_T_bounds=ytb.copy()
          
#      xtb0=2.0*self.x_T_bounds[:,0]-self.x_T_bounds[:,1]
#      xtb0=xtb0[:,np.newaxis]
#      self.x_T_bounds = np.hstack((xtb0,self.x_T_bounds))
#      xtb0=self.x_T_bounds[0,:]
#      self.x_T_bounds = np.vstack((xtb0,self.x_T_bounds))          
#      ytb0=2.0*self.y_T_bounds[0,:]-self.y_T_bounds[1,:]
#      self.y_T_bounds = np.vstack((ytb0,self.y_T_bounds))
#      ytb0=self.y_T_bounds[:,0]
#      ytb0=ytb0[:,np.newaxis]
#      self.y_T_bounds = np.hstack((ytb0,self.y_T_bounds))
      dx = (self.x_T_bounds - np.roll(self.x_T_bounds,axis=1,shift=1))
      dx=dx*np.pi/180.
      dx[:,0] = dx[:,1]
      self.dxh = dx*R_earth*np.cos(self.y_T*np.pi/180.)
      dy = (self.y_T_bounds - np.roll(self.y_T_bounds,axis=0,shift=1))
      dy = dy*np.pi/180.
      dy[0,:]=dy[1,:]
      
      self.dyh = dy*R_earth

      self.Ah=self.dxh*self.dyh


      self.cyclic_x = cyclic

      if self.cyclic_x:
        self.xmod_len = self.x_T_bounds[0,-1] - self.x_T_bounds[0,0] 


  def geo_region(self,y=None,x=None,name=None):
     [min_dx,min_dy] = min_resolution(self)
     [max_dx,max_dy] = max_resolution(self)  
  
     section={}

     if y is not None:
         ys = np.nonzero(np.abs(self.lath - y[0]) <= max_dy)[0][0]
         ye = np.nonzero(np.abs(self.lath - y[1]) <= max_dy)[0][0]  
         section['y']=np.arange(ys,ye+1)
         section['yax_data']= self.lath[section['y']]
     else:
         section['y']=None

     x_T,y_T = np.meshgrid(self.lonh,self.lath)
    
     if x is not None:
         if x[0] >= self.lonh[0] and x[1] <= self.lonh[-1]:
             xs,xe,ys,ye = find_geo_bounds(self,x,y)
             section['x']=np.arange(xs,xe+1)
             section['lon0']=x[0]
             section['x_offset'] = 0
             section['shifted']=False
             section['xax_data']= self.lonh[section['x']]
         elif x[0] < self.lonh[0]:
             result=find_axis_bounds(self.lonh-360.,x=[x[0],x[0]])
             section['x_offset']=result[0]        
             x_T_shifted,lon_shifted = shiftgrid(x[0],x_T,self.lonh)
             result=find_axis_bounds(lon_shifted,x=x)
             xs=result[0];xe=result[1]
             section['lon0']=x[0]
             section['shifted']=True
             section['x']=np.arange(xs,xe+1)        
             section['xax_data']= lon_shifted[section['x']][xs:xe+1]
         elif x[0] > self.lonh[-1]:
             result=find_axis_bounds(self.lonh+360.,x=[x[0],x[0]])        
             section['x_offset']=result[0]        
             x_T_shifted,lon_shifted = shiftgrid(x[0],x_T,self.lonh)
             result=find_axis_bounds(lon_shifted,x=x)
             xs=result[0];xe=result[1]
             section['lon0']=x[0]
             section['shifted']=True
             section['x']=np.arange(xs,xe+1)        
             section['xax_data']= lon_shifted[section['x']][xs:xe+1]
         else:
             result=find_axis_bounds(self.lonh,x=[x[0],x[0]])        
             section['x_offset']=result[0]        
             x_T_shifted,lon_shifted = shiftgrid(x[0],x_T,self.lonh)
             result=find_axis_bounds(lon_shifted,x=x)
             xs=result[0];xe=result[1]
             section['lon0']=x[0]
             section['shifted']=True
             section['x']=np.arange(xs,xe+1)        
             section['xax_data']= lon_shifted[section['x']][xs:xe+1]
     else:
         im=self.lonh.shape[0]
         section['x']=np.arange(0,im)
         section['lon0']=self.lonh[0]
         section['x_offset'] = 0
         section['shifted']=False
         section['xax_data']= self.lonh

     section['name'] = name
     section['parent_grid'] = self

     return section

  def indexed_region(self,i=None,j=None,name=None):
  
    section={}

    if j is not None:
      ys=j[0];ye=j[1]
      section['y']=np.arange(ys,ye+1)
      section['yax_data']= self.lath[section['y']]
    else:
      section['y']=None
      
    if i is not None:
        xs=i[0];xe=i[1]
        section['x']=np.arange(xs,xe+1)
        section['xax_data']= self.lonh[section['x']]
        section['lon0']=section['xax_data'][0]
        section['x_offset'] = 0
        section['shifted']=False
    else:
      section['x']=None

    section['name'] = name
    section['parent_grid'] = self

    return section
      
  def extract(self,geo_region=None):
    if geo_region is None:
        grid=copy.copy(self)
        return grid
    else:
      grid = copy.copy(self)

      x_section = geo_region['x'];y_section = geo_region['y']      
      if x_section is not None:
          xb_section = np.hstack((x_section,x_section[-1]+1))
      if y_section is not None:
          yb_section = np.hstack((y_section,y_section[-1]+1))        

      grid.lath = np.take(self.lath,y_section,axis=0)
      grid.latq = np.take(self.latq,yb_section,axis=0)

      if self.simple_grid:
          x,y=np.meshgrid(self.lonh,self.lath)
          if geo_region['shifted']:
              x_shifted,lon_shifted = shiftgrid(geo_region['lon0'],x,self.lonh)
              grid.lonh=np.take(lon_shifted,x_section,axis=0)
              x_shifted,lonq_shifted = shiftgrid(geo_region['lon0'],x,self.lonq)
              grid.lonq=np.take(lonq_shifted,xb_section,axis=0)              
          else:
              grid.lonh=np.take(grid.lonh,x_section,axis=0)
              grid.lonq=np.take(grid.lonq,xb_section,axis=0)
          return grid
      
      if geo_region['shifted']:
        x_T_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.x_T,grid.lonh)
        grid.x_T = np.take(np.take(x_T_shifted,y_section,axis=0),x_section,axis=1)
        y_T_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.y_T,grid.lonh)        
        grid.y_T = np.take(np.take(y_T_shifted,y_section,axis=0),x_section,axis=1)        
        lonb = grid.x_T_bounds[0,:]
        x_T_bounds_shifted,lon_bounds_shifted = shiftgrid(geo_region['lon0'],grid.x_T_bounds,lonb)
        grid.x_T_bounds = np.take(np.take(x_T_bounds_shifted,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1)
        grid.lonq = np.take(lon_bounds_shifted,xb_section,axis=0)        
        y_T_bounds_shifted,lon_bounds_shifted = shiftgrid(geo_region['lon0'],grid.y_T_bounds,lonb)        
        grid.y_T_bounds = np.take(np.take(y_T_bounds_shifted,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1)

        dxh_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dxh,grid.lonh)                
        grid.dxh = np.take(np.take(dxh_shifted,y_section,axis=0),x_section,axis=1)
        dyh_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dyh,grid.lonh)                
        grid.dyh = np.take(np.take(dyh_shifted,y_section,axis=0),x_section,axis=1)

        try:
          grid.mask=np.roll(grid.mask,axis=1,shift=-geo_region['x_offset'])
        except:
          pass
      else:
        grid.x_T = np.take(np.take(self.x_T,y_section,axis=0),x_section,axis=1)
        grid.y_T = np.take(np.take(self.y_T,y_section,axis=0),x_section,axis=1)
        grid.x_T_bounds = np.take(np.take(self.x_T_bounds,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1)
        grid.y_T_bounds = np.take(np.take(self.y_T_bounds,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1) 
        grid.lonh = np.take(self.lonh,x_section,axis=0)
        grid.lath = np.take(self.lath,y_section,axis=0)
        grid.lonq = np.take(self.lonq,xb_section,axis=0)
        grid.latq = np.take(self.latq,yb_section,axis=0)
        grid.dxh = np.take(np.take(self.dxh,y_section,axis=0),x_section,axis=1)
        grid.dyh = np.take(np.take(self.dyh,y_section,axis=0),x_section,axis=1)
          
      grid.im = np.shape(grid.lonh)[0]
      grid.jm = np.shape(grid.lath)[0]

      try:
        grid.mask=np.take(np.take(grid.mask,geo_region['y'],axis=0),geo_region['x'],axis=1)
      except:
        pass
                
      
      return grid

  def add_mask(self,field,path=None):
    """Add a 2-D mask to the grid. The mask can have any values, e.g.
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

class gold_grid(object):

  
  """A gold grid object is for the horizontal lattice of points lying
     on the sphere. A gold_grid class is constructed by reading the
     contents of an \"ocean_geometry\" file produced by the GFDL
     General Ocean Layered Dynamics (GOLD) model. Alternatively, the
     cell positions can be constructed from a MOM4 grid_spec file.

     The cell perimeter locations, face lengths and areas are
     read from the geometry file. 

           +-------+-------+-------+
           :       :       :       :
           :       :       :       :
           :       :       :       :           
           +-------+-------Q-------:< Y_T_bounds
           :       :       :       :
           :       :   T   :       :
           :       :       :       :           
           +-------Q-------+-------+< Y_T_bounds
                   ^       ^
                  X_T_bounds

     A grid object describes the horizontal discretization of the model
     along with its topology (for instance how the edges are connected in
     to a mosaic)

                  
  """


  
  
  def __init__(self,path='ocean_geometry.nc',cyclic=False,tripolar_n=False,grid_type='gold_geometry'):
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
        self.dxv    = f.variables['dxv'][:]
        self.dyu    = f.variables['dyu'][:]
        self.dxu    = f.variables['dxu'][:]
        self.dyv    = f.variables['dyv'][:]
        self.dxh    = f.variables['dxh'][:]
        self.dyh    = f.variables['dyh'][:]
        self.dxq    = f.variables['dxq'][:]
        self.dyq    = f.variables['dyq'][:]
        self.Ah    = f.variables['Ah'][:]
        self.Aq    = f.variables['Aq'][:]
        self.wet    = f.variables['wet'][:]
        self.cyclic_x = cyclic
        self.tripolar_n = tripolar_n
        self.im = np.shape(self.lonh)[0]
        self.jm = np.shape(self.lath)[0]
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
        self.dxq    = f.variables['dxu'][:]
        self.dyq    = f.variables['dyu'][:]
        self.Ah    = self.dxh*self.dyh
        self.Aq    = self.dxq*self.dyq
        self.wet    = f.variables['wet'][:]
        self.cyclic_x = cyclic
        self.tripolar_n = tripolar_n
        self.im = np.shape(self.lonh)[0]
        self.jm = np.shape(self.lath)[0]
        
  def geo_region(self,y=None,x=None,name=None):
    [min_dx,min_dy] = min_resolution(self)
    [max_dx,max_dy] = max_resolution(self)  
  
    section={}

    if y is not None:
      ys = np.nonzero(np.abs(self.lath - y[0]) <= max_dy)[0][0]
      ye = np.nonzero(np.abs(self.lath - y[1]) <= max_dy)[0][0]  
      section['y']=np.arange(ys,ye+1)
      section['yax_data']= self.lath[section['y']]
    else:
      section['y']=None
      
    if x is not None:
      if x[0] >= self.lonh[0] and x[1] <= self.lonh[-1]:
        xs,xe,ys,ye = find_geo_bounds(self,x,y)
        section['x']=np.arange(xs,xe+1)
        section['lon0']=x[0]
        section['x_offset'] = 0
        section['shifted']=False
        section['xax_data']= self.lonh[section['x']]
      elif x[0] < self.lonh[0]:
        result=find_axis_bounds(self.lonh-360.,x=[x[0],x[0]])
        section['x_offset']=result[0]        
        x_T_shifted,lon_shifted = shiftgrid(x[0],self.x_T,self.lonh)
        result=find_axis_bounds(lon_shifted,x=x)
        xs=result[0];xe=result[1]
        section['lon0']=x[0]
        section['shifted']=True
        section['x']=np.arange(xs,xe+1)        
        section['xax_data']= lon_shifted[section['x']][xs:xe+1]
      elif x[0] > self.lonh[-1]:
        result=find_axis_bounds(self.lonh+360.,x=[x[0],x[0]])        
        section['x_offset']=result[0]        
        x_T_shifted,lon_shifted = shiftgrid(x[0],self.x_T,self.lonh)
        result=find_axis_bounds(lon_shifted,x=x)
        xs=result[0];xe=result[1]
        section['lon0']=x[0]
        section['shifted']=True
        section['x']=np.arange(xs,xe+1)        
        section['xax_data']= lon_shifted[section['x']][xs:xe+1]
      else:
        result=find_axis_bounds(self.lonh,x=[x[0],x[0]])        
        section['x_offset']=result[0]        
        x_T_shifted,lon_shifted = shiftgrid(x[0],self.x_T,self.lonh)
        result=find_axis_bounds(lon_shifted,x=x)
        xs=result[0];xe=result[1]
        section['lon0']=x[0]
        section['shifted']=True
        section['x']=np.arange(xs,xe+1)        
        section['xax_data']= lon_shifted[section['x']][xs:xe+1]
    else:
      section['x']=None

    section['name'] = name
    section['parent_grid'] = self

    return section
      
  def indexed_region(self,i=None,j=None,name=None):
  
    section={}

    if j is not None:
      ys=j[0];ye=j[1]
      section['y']=np.arange(ys,ye+1)
      section['yax_data']= self.lath[section['y']]
    else:
      section['y']=None
      
    if i is not None:
        xs=i[0];xe=i[1]
        section['x']=np.arange(xs,xe+1)
        section['xax_data']= self.lonh[section['x']]
        section['lon0']=section['xax_data'][0]
        section['x_offset'] = 0
        section['shifted']=False
    else:
      section['x']=None

    section['name'] = name
    section['parent_grid'] = self

    return section
      
  def extract(self,geo_region=None):
    if geo_region is None:
        grid=copy.copy(self)
        return grid
    else:
      grid = copy.copy(self)

      x_section = geo_region['x'];y_section = geo_region['y']
      if x_section is not None:
        xb_section = np.hstack((x_section,x_section[-1]+1))
      if y_section is not None:
        yb_section = np.hstack((y_section,y_section[-1]+1))        

      if not geo_region['shifted']:
        grid.x_T = np.take(np.take(self.x_T,y_section,axis=0),x_section,axis=1)
        grid.y_T = np.take(np.take(self.y_T,y_section,axis=0),x_section,axis=1)
        grid.x_T_bounds = np.take(np.take(self.x_T_bounds,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1)
        grid.y_T_bounds = np.take(np.take(self.y_T_bounds,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1) 
        grid.lonh = np.take(self.lonh,x_section,axis=0)
        grid.lath = np.take(self.lath,y_section,axis=0)
        grid.lonq = np.take(self.lonq,xb_section,axis=0)
        grid.latq = np.take(self.latq,yb_section,axis=0)
        grid.D = np.take(np.take(self.D,y_section,axis=0),x_section,axis=1)
        grid.f = np.take(np.take(self.f,y_section,axis=0),x_section,axis=1)

        if hasattr(grid,'dxv'):
            grid.dxv = np.take(np.take(self.dxv,y_section,axis=0),x_section,axis=1)
            grid.dyu = np.take(np.take(self.dyu,y_section,axis=0),x_section,axis=1)
            grid.dxu = np.take(np.take(self.dxu,y_section,axis=0),x_section,axis=1)
            grid.dyv = np.take(np.take(self.dyv,y_section,axis=0),x_section,axis=1)
            
        grid.dxh = np.take(np.take(self.dxh,y_section,axis=0),x_section,axis=1)
        grid.dyh = np.take(np.take(self.dyh,y_section,axis=0),x_section,axis=1)
        grid.dxq = np.take(np.take(self.dxq,y_section,axis=0),x_section,axis=1)
        grid.dyq = np.take(np.take(self.dyq,y_section,axis=0),x_section,axis=1)
        grid.Ah = np.take(np.take(self.Ah,y_section,axis=0),x_section,axis=1)
        grid.Aq= np.take(np.take(self.Aq,y_section,axis=0),x_section,axis=1)
        grid.wet = np.take(np.take(self.wet,y_section,axis=0),x_section,axis=1)
        grid.im = np.shape(grid.lonh)[0]
        grid.jm = np.shape(grid.lath)[0]

        try:
          grid.mask=np.take(np.take(grid.mask,y_section,axis=0),x_section,axis=1)
        except:
          pass
        
      else:
        x_T_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.x_T,grid.lonh)
        grid.x_T = np.take(np.take(x_T_shifted,y_section,axis=0),x_section,axis=1)
        y_T_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.y_T,grid.lonh)        
        grid.y_T = np.take(np.take(y_T_shifted,y_section,axis=0),x_section,axis=1)        
        lonb = grid.x_T_bounds[0,:]
        x_T_bounds_shifted,lon_bounds_shifted = shiftgrid(geo_region['lon0'],grid.x_T_bounds,lonb)
        grid.x_T_bounds = np.take(np.take(x_T_bounds_shifted,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1)
        grid.lonq = np.take(lon_bounds_shifted,xb_section,axis=0)        
        y_T_bounds_shifted,lon_bounds_shifted = shiftgrid(geo_region['lon0'],grid.y_T_bounds,lonb)        
        grid.y_T_bounds = np.take(np.take(y_T_bounds_shifted,np.hstack((y_section,y_section[-1]+1)),axis=0),np.hstack((x_section,x_section[-1]+1)),axis=1)
        grid.lath = np.take(grid.lath,yb_section,axis=0)
        grid.latq = np.take(self.latq,yb_section,axis=0)
        D_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.D,grid.lonh)
        grid.D = np.take(np.take(D_shifted,y_section,axis=0),x_section,axis=1)
        grid.f = np.take(np.take(self.f,y_section,axis=0),x_section,axis=1)
        dxv_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dxv,grid.lonh)
        grid.dxv = np.take(np.take(dxv_shifted,y_section,axis=0),x_section,axis=1)
        dyu_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dyu,grid.lonh)
        grid.dyu = np.take(np.take(dyu_shifted,y_section,axis=0),x_section,axis=1)
        dxu_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dxu,grid.lonh)
        grid.dxu = np.take(np.take(dxu_shifted,y_section,axis=0),x_section,axis=1)
        dyv_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dyv,grid.lonh)        
        grid.dyv = np.take(np.take(dyv_shifted,y_section,axis=0),x_section,axis=1)
        dxh_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dxh,grid.lonh)                
        grid.dxh = np.take(np.take(dxh_shifted,y_section,axis=0),x_section,axis=1)
        dyh_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dyh,grid.lonh)                
        grid.dyh = np.take(np.take(dyh_shifted,y_section,axis=0),x_section,axis=1)
        dxq_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dxq,grid.lonh)                
        grid.dxq = np.take(np.take(dxq_shifted,y_section,axis=0),x_section,axis=1)
        dyq_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.dyq,grid.lonh)                        
        grid.dyq = np.take(np.take(dyq_shifted,y_section,axis=0),x_section,axis=1)
        Ah_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.Ah,grid.lonh)                
        grid.Ah = np.take(np.take(Ah_shifted,y_section,axis=0),x_section,axis=1)
        Aq_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.Aq,grid.lonh)                
        grid.Aq= np.take(np.take(Aq_shifted,y_section,axis=0),x_section,axis=1)
        wet_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.wet,grid.lonh)                
        grid.wet = np.take(np.take(wet_shifted,y_section,axis=0),x_section,axis=1)
        try:
          mask_shifted,lon_shifted = shiftgrid(geo_region['lon0'],grid.mask,grid.lonh)                
          grid.mask = np.take(np.take(mask_shifted,y_section,axis=0),x_section,axis=1)
        except:
          pass
        
        grid.lonh = np.take(lon_shifted,x_section,axis=0)        
        grid.im = np.shape(grid.lonh)[0]
        grid.jm = np.shape(grid.lath)[0]        


          
      return grid

  def add_mask(self,field,path=None):
    """Add a 2-D mask to the grid. The mask can have arbitrary values, e.g.
       a different value for separate ocean basins. This can be used to mask
       fields using the (state.mask_where) method.
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

  (TODO:***zlevels***) is not yet implemented.

  (NOTE:The layer (interfaces) variable must be stored in order to calculate
  accurate finite volume integrals.
  
     """
  
  def __init__(self,path=None,grid=None,geo_region=None,time_indices=None,date_bounds=None,z_indices=None,zlevels=None,fields=None,default_calendar=None,MFpath=None,interfaces=None,path_interfaces=None,MFpath_interfaces=None):

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

# Loop through variables in list and bind them to the state

    if zlevels is not None:
      print """
       Use of zlevels to bracket a region is not yet implemented."""
      raise

    z_indices_in = z_indices

    if fields is None:
      if grid is not None:
        new_grid = grid.extract(geo_region)
        self.grid = new_grid
      return
        
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

       var_dict['is_MFpath']=self.is_MFpath

       for n in range(0,self.variables[v].ndim):

# For each dimension , determine its Cartesian attribute
# This is a critical step. Along with the Cartesian
# attribute, the distance along each axis is determined at
# this stage.


         dimnam = self.variables[v].dimensions[n]
         dim = self.rootgrp.variables[dimnam]
         cart = get_axis_cart(dim,dimnam)

         if cart is not None:
           var_dict[cart]=self.variables[v].dimensions[n]
           
         if cart == 'Z':
           var_dict['Zdir'] = get_axis_direction(dim)
           zunits = getattr(dim,'units')
           var_dict['zunits']=string.lower(zunits)
           if var_dict['zunits'] in ('cm','m','meters','km','pa','hpa'):
             var_dict['Ztype'] = 'Fixed'
           else:
             var_dict['Ztype'] = 'Lagrangian'

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
           var_dict['tunits'] = self.rootgrp.variables[var_dict['T']].units
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
               var_dict['date_bounds']=var_dict['date_bounds'][tb_indices,0]
         elif date_bounds is not None:
           if var_dict['calendar'] is not None:
             ts,te = find_date_bounds(var_dict['date_bounds'][:],date_bounds[0],date_bounds[1])
             t_indices=np.arange(ts,te+1)
             tb_indices=np.arange(ts,te+2)             
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
         var_dict['zunits'] =   self.rootgrp.variables[var_dict['Z']].units           
         if z_indices is not None:
           tmp = np.array([z_indices[-1]+1])
           z_interfaces = np.hstack((z_indices,tmp))
           nz = len(z_indices)
         else:
           nz = len(self.rootgrp.variables[var_dict['Z']][:])
           z_indices=np.arange(0,nz)
           z_interfaces = np.arange(0,nz+1)

         var_dict['zax_data']= self.rootgrp.variables[var_dict['Z']][z_indices]
       
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

       
       if var_dict['Y'] is not None and var_dict['X'] is not None:
         var_dict['xunits'] =   self.rootgrp.variables[var_dict['X']].units
         var_dict['yunits'] =   self.rootgrp.variables[var_dict['Y']].units                                                       
         if geo_region is not None:
           var_dict['yax_data'] = geo_region['yax_data'][:]
           y_indices=geo_region['y']
           if geo_region['shifted']:
             nx = len(geo_region['x'][:])
             nx_g = len(self.rootgrp.variables[var_dict['X']][:])
             x_indices_read = np.arange(0,nx_g)
             x_indices = np.arange(0,nx)             
             var_dict['xax_data'] = geo_region['xax_data'][:]
           else:
             var_dict['xax_data'] = geo_region['xax_data'][:]
             nx = len(geo_region['x'][:])
             x_indices = np.arange(0,nx)
             x_indices_read = geo_region['x'][:]
         else:
           ny = len(self.rootgrp.variables[var_dict['Y']][:])
           y_indices = np.arange(0,ny)
           var_dict['yax_data'] = self.rootgrp.variables[var_dict['Y']][y_indices]
           nx = len(self.rootgrp.variables[var_dict['X']][:])
           x_indices = np.arange(0,nx)
           var_dict['xax_data'] = self.rootgrp.variables[var_dict['X']][x_indices]
           x_indices_read = x_indices
           
       if var_dict['X'] is None:
         x_indices = None
         x_indices_read = None         

       if var_dict['Y'] is None:
         y_indices = None         
         y_indices_read = None

       slice_read = [];shape_read = []; shape_out=[]; slice_out = []
       for s in [t_indices,z_indices,y_indices]:
         if s is not None:
           slice_read.append(s)
           slice_out.append(np.arange(0,s.shape[0]))             
           shape_read.append(s.shape[0])
           shape_out.append(s.shape[0])           
         else:
           slice_out.append([0])
           shape_out.append(1)
           shape_read.append(1)

       if geo_region is not None and x_indices is not None:
         slice_read.append(x_indices_read)
         slice_out.append(x_indices)
         shape_read.append(x_indices_read.shape[0])
         shape_out.append(x_indices.shape[0])
       elif x_indices is not None:
         slice_read.append(x_indices_read)
         slice_out.append(x_indices)         
         shape_read.append(x_indices_read.shape[0])
         shape_out.append(x_indices.shape[0])
       else:
         slice_out.append([0])
         shape_out.append(1)
         shape_read.append(1)

       slice_int_read = [];shape_int_read = []; shape_int_out=[]; slice_int_out = []
       for s in [t_indices,z_interfaces,y_indices]:
         if s is not None:
           slice_int_read.append(s)
           slice_int_out.append(np.arange(0,s.shape[0]))                      
           shape_int_read.append(s.shape[0])
           shape_int_out.append(s.shape[0])           
         else:
           slice_int_out.append([0])
           shape_int_out.append(1)
           shape_int_read.append(1)

       if geo_region is not None and x_indices is not None:
         slice_int_read.append(x_indices_read)
         slice_int_out.append(x_indices)
         shape_int_read.append(x_indices_read.shape[0])
         shape_int_out.append(x_indices.shape[0])         
       elif x_indices is not None:
         slice_int_read.append(x_indices_read)
         slice_int_out.append(x_indices)             
         shape_int_read.append(x_indices_read.shape[0])
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


       data_read = np.reshape(np.array(self.rootgrp.variables[v][slice_read]),(shape_read))

       if geo_region is not None:
         if geo_region['shifted']:
           vars(self)[v] = np.roll(data_read,axis=3,shift=-geo_region['x_offset'])[:,:,:,x_indices]
         else:
           vars(self)[v] = data_read
       else:
         vars(self)[v] = data_read


       if DEBUG == 1:
         print " Successfully extracted data named %(nam)s from %(fil)s "%{'nam':v,'fil':self.path}
         print " Resulting shape = ",vars(self)[v].shape
         print " Dictionary keys = ", var_dict.keys()

                  


       var_dict['masked']=False
       var_dict['FillValue'] = None
       var_dict['missing_value'] = None

         
       for i in f.variables[v].ncattrs():
         if i == 'units':
           var_dict['units']=f.variables[v].units
         if i == '_FillValue':
           var_dict['FillValue'] = f.variables[v]._FillValue
         if i == 'missing_value':
           var_dict['missing_value'] = f.variables[v].missing_value

       if var_dict['FillValue'] is not None or var_dict['missing_value']  is not None:
         if var_dict['FillValue'] is not None:
             vars(self)[v] = np.ma.masked_where(np.abs(vars(self)[v] - var_dict['FillValue']) < epsln,vars(self)[v])
         else:
             vars(self)[v] = np.ma.masked_where(np.abs(vars(self)[v] - var_dict['missing_value']) < np.abs(var_dict['missing_value']*0.01),vars(self)[v])


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
           
           if geo_region is not None:
             if geo_region['shifted']:
               vars(self)[interfaces] = np.roll(data_int_read,axis=3,shift=-geo_region['x_offset'])[:,:,:,x_indices]
             else:
               vars(self)[interfaces] = data_int_read
           else:
             vars(self)[interfaces] = data_int_read
         
           self.interfaces=interfaces
           
           if DEBUG == 1:
             print " Successfully extracted interface data named %(nam)s from %(fil)s "%{'nam':interfaces,'fil':path_interfaces}
             print " Resulting shape = ",vars(self)[interfaces].shape


       else:
           self.interfaces = None
           
       if var_dict['Z'] is not None and var_dict['Ztype'] is not 'Lagrangian':
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
             
       if var_dict['Z'] is not None and var_dict['Ztype'] is 'Lagrangian' and interfaces is not None:           
         var_dict['z_interfaces']  = vars(self)[interfaces]
         tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=1,shift=-1))             
         var_dict['z'] = tmp[:,0:-1,:,:]
         var_dict['dz'] = var_dict['z_interfaces'] - np.roll(var_dict['z_interfaces'],axis=1,shift=-1)
         var_dict['dz'] = var_dict['dz'][:,0:-1,:,:]
         
       try:
         time_avg_info = getattr(self.rootgrp.variables[v],'time_avg_info')
       except:
         time_avg_info = ""
         
       if 'average_DT' in time_avg_info:
         dt = self.rootgrp.variables['average_DT'][t_indices]
         units = self.rootgrp.variables['average_DT'].units
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

       var_dict['rootgrp'] = self.rootgrp
          
       self.var_dict[v] = var_dict
       self.geo_region = geo_region

       del t_indices,y_indices,x_indices
       z_indices = z_indices_in

    
    if grid is not None:
        new_grid = grid.extract(geo_region)
        self.grid = new_grid

         
  def add_field(self,field,path=None,MFpath=None):
    """Add a field to the existing state (e.g. more tracers).
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
      if cart == 'Z':
        var_dict['Zdir'] = get_axis_direction(dim)
        zunits = getattr(dim,'units')
        if zunits in ('cm','m','meters','km','Pa','HPa'):
          var_dict['Ztype'] = 'Fixed'
        else:
          var_dict['Ztype'] = 'Lagrangian'
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
           
        if t_indices is not None:
            tb_indices = np.hstack((t_indices,t_indices[-1]+1))
            if var_dict['tbax_data'] is not None:
                var_dict['tbax_data']=var_dict['tbax_data'][tb_indices]
            if var_dict['tbax_data'] is not None:
                var_dict['date_bounds']=var_dict['date_bounds'][tb_indices]
      
    if var_dict['Z'] is not None:
        z_indices = self.slice_read[1]
        nz=len(z_indices)
        z_interfaces = self.slice_int_read[1]
        var_dict['zax_data']= f.variables[var_dict['Z']][z_indices]
        if var_dict['Zb'] is not None:
            var_dict['zbax_data'] = f.variables[var_dict['Zb']][z_interfaces]
        else:
            zdat=f.variables[var_dict['Z']][:]
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

    if geo_region is not None:
      if geo_region['shifted']:
        x_indices = geo_region['x_indices'][:]
        vars(self)[field] = np.roll(data_read,axis=3,shift=-geo_region['x_offset'])[:,:,:,x_indices]
      else:
        vars(self)[field] = data_read
    else:
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

         
    if DEBUG == 1:
      print " Successfully extracted data named %(nam)s from %(fil)s "%{'nam':field,'fil':var_dict['path']}
      print " Resulting shape = ",vars(self)[field].shape                 


    var_dict['masked']=False
    var_dict['FillValue'] = None
    var_dict['missing_value'] = None

    
    for i in f.variables[field].ncattrs():
      if i == 'units':
        var_dict['units']=f.variables[field].units      
      if i == '_FillValue':
        var_dict['FillValue'] = f.variables[field]._FillValue
      if i == 'missing_value':
        var_dict['missing_value'] = f.variables[field].missing_value

    if var_dict['FillValue'] is not None or var_dict['missing_value']  is not None:
      if var_dict['FillValue'] is not None:
        vars(self)[field] = np.ma.masked_where(np.abs(vars(self)[field] - var_dict['FillValue']) < epsln,vars(self)[field])
      else:
        vars(self)[field] = np.ma.masked_where(np.abs(vars(self)[field] - var_dict['missing_value']) < epsln,vars(self)[field])             

      var_dict['masked']=True

    
    var_dict['z_interfaces'] = None
       
    if var_dict['Z'] is not None:
      if var_dict['zbax_data'] is not None and var_dict['Ztype'] != 'Lagrangian':
        tmp = np.reshape(var_dict['zbax_data'][z_interfaces],(nz+1,1,1))
        if geo_region is not None:
          ny = len(geo_region['y'])
          nx = len(geo_region['x'])          
        else:
          ny = len(f.variables[var_dict['Y']][:])
          nx = len(f.variables[var_dict['X']][:])
          
        var_dict['z_interfaces']  = var_dict['Zdir']*np.tile(tmp,(1,ny,nx))
        tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))             
        var_dict['z'] = tmp[0:-1,:,:]
        var_dict['dz'] = var_dict['z_interfaces'] - np.roll(var_dict['z_interfaces'],axis=0,shift=-1)
        var_dict['dz'] = var_dict['dz'][0:-1,:,:]
      if var_dict['Z'] is not None and var_dict['Ztype'] is 'Lagrangian' and self.interfaces is not None:           
        var_dict['z_interfaces']  = vars(self)[self.interfaces]
        tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=1,shift=-1))             
        var_dict['z'] = tmp[:,0:-1,:,:]
        var_dict['dz'] = var_dict['z_interfaces'] - np.roll(var_dict['z_interfaces'],axis=1,shift=-1)
        var_dict['dz'] = var_dict['dz'][:,0:-1,:,:]
         
               
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
    """Define a new field and add to the existing state.
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
    cmd = string.join(['self.',field_new,'=self.',field,'.copy()'],sep='')
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

    vars(self)[field].mask = False

  def add_interface_bounds(self,field=None):
    """
    Add interfaces at cell boundaries (ny+1,nx+1).

    ***NOTE: Assumes cyclic x at this point.
    ***TODO
    
    """
    
    if field is None:
      return None

    e=self.var_dict[field]['z_interfaces']
    eb = 0.5*(e+np.roll(e,shift=-1,axis=3))
    self.var_dict[field]['z_interfaces_ew']=np.concatenate((np.take(eb,[-1],axis=3),eb),axis=3)
    eb = 0.5*(e+np.roll(e,shift=-1,axis=2))
    self.var_dict[field]['z_interfaces_ns']=np.concatenate((np.take(eb,[0],axis=2),eb),axis=2)    
    
  def fill_interior(self,field=None):
    """

    Fill interior above the topography .
    
    """

    import vertmap as vmap

    FVal_=-1.e10

  
    if field is None:
      print """
       Please specify (all) : field, numpass,crit and relc
      """
      return None

    if self.var_dict[field]['masked'] is False:
      print """
       Input field needs to be masked for sor_fill
      """
      return None
    
    val = vars(self)[field]
    val_prev = np.zeros([val.shape[2],val.shape[3]])    
    wet = self.grid.wet


    for i in np.arange(0,val.shape[0]):
      for j in np.arange(0,val.shape[1]):
        zbot = self.var_dict[field]['z_interfaces'][j+1,:]
        ztop = self.var_dict[field]['z_interfaces'][j,:]        
        tmp = np.squeeze(np.take(np.take(val,[i],axis=0),[j],axis=1))
        mask_in = np.ma.getmask(tmp)

        # fill has a value of one over points
        # which are in the interior at the current depth
        # and zero otherwise
        
        fill = np.zeros([tmp.shape[0],tmp.shape[1]])
        fill[np.logical_and(np.logical_and(wet == 1.0,-ztop < self.grid.D),mask_in) ]=1.0

        mask_out=np.logical_or(wet == 0.0,-ztop > self.grid.D)

        good = np.zeros([tmp.shape[0],tmp.shape[1]])
        good[~mask_in]=1


        v_filled = np.zeros([tmp.shape[1],tmp.shape[0]])

        if j>0:
          # initialize with nearest fill or value at previous level
          v_filled=vmap.midas_vertmap.fill_miss_2d(tmp.T,good.T,fill.T,val_prev.T,cyclic_x=True,tripolar_n=True)
        else:
          v_filled=vmap.midas_vertmap.fill_miss_2d(tmp.T,good.T,fill.T,cyclic_x=True,tripolar_n=True)

        
        v_filled=v_filled.T

        v_filled[mask_out==1]=FVal_
        v_filled=np.ma.masked_where(mask_out==1,v_filled)

        val_prev = v_filled.copy()
        

        val[i,j,:]=v_filled[:]


        

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
    
    
    cmd = string.join(['sout=self.',field],sep='')
    exec(cmd)
    
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
      var_dict['zax_data']=[np.mean(var_dict['zax_data'])]
      var_dict['zbax_data']=[var_dict['zbax_data'][0],var_dict['zbax_data'][-1]]

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
        var_dict['dz'] = np.sum(dz_masked*dy_masked*dx_masked,axis=2)/np.sum(dy_masked*dx_masked,axis=2)      
        var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],sout.shape[1],1,sout.shape[3]))      

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
        
        var_dict['dz'] = np.sum(dz_masked*dy_masked*dx_masked,axis=3)/np.sum(dy_masked*dx_masked,axis=3)      
        var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],sout.shape[1],sout.shape[2],1))      
        z0 = np.take(var_dict['z_interfaces'],[0],axis=0)
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
        var_dict['dz'] = np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=2)/np.sum(np.sum(dy_masked*dx_masked,axis=3),axis=2)
        var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],sout.shape[1],1,1))      
        z0 = np.take(var_dict['z_interfaces'],[0],axis=1)
        dy0 = np.take(dy_masked,[0],axis=1)
        dx0 = np.take(dx_masked,[0],axis=1)      
        result = np.sum(np.sum(z0*dy0*dx0,axis=3),axis=2)/np.sum(np.sum(dy0*dx0,axis=3),axis=2)
        result = np.reshape(result,(sout.shape[0],1,1,1))
        var_dict['z_interfaces']=-np.cumsum(var_dict['dz'],axis=1)
        var_dict['z_interfaces']=np.concatenate((result,var_dict['z_interfaces']),axis=1)
        tmp=0.5*(var_dict['z_interfaces'] + np.roll(var_dict['z_interfaces'],axis=1,shift=-1))
        tmp = tmp[:,0:-1,:,:]        
        var_dict['rootgrp']=None
      

        if var_dict['Ztype'] is 'Fixed':
          var_dict['z']=np.take(tmp,[0],axis=0)
          var_dict['dz']=np.take(var_dict['dz'],[0],axis=0)
          var_dict['z_interfaces']=np.take(var_dict['z_interfaces'],[0],axis=0)                
        else:
          var_dict['z']=tmp

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
      
      var_dict['dz'] = np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=1)/np.sum(np.sum(dy_masked*dx_masked,axis=3),axis=1)
      var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],1,sout.shape[2],1))            

      z0 = np.take(var_dict['z_interfaces'],[0],axis=1)
      dy0 = np.take(dy_masked,[0],axis=1)
      dx0 = np.take(dx_masked,[0],axis=1)      
      result = np.sum(z0*dy0*dx0,axis=3)/np.sum(dy0*dx0,axis=3)
      result = np.reshape(result,(sout.shape[0],1,sout.shape[2],1))
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
      
      var_dict['dz'] = np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=2),axis=1)/np.sum(np.sum(dy_masked*dx_masked,axis=2),axis=1)
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
      
      var_dict['dz'] = np.sum(np.sum(np.sum(dz_masked*dy_masked*dx_masked,axis=3),axis=2),axis=1)/np.sum(np.sum(np.sum(dy_masked*dx_masked,axis=3),axis=2),axis=1)
      var_dict['dz']=np.reshape(var_dict['dz'],(sout.shape[0],1,1,1))                  

      z0 = np.take(var_dict['z_interfaces'],[0],axis=1)
      dy0 = np.take(dy_masked,[0],axis=1)
      dx0 = np.take(dx_masked,[0],axis=1)      
      result = np.sum(np.sum(z0*dy0*dx0,axis=3),axis=2)/np.sum(np.sum(dy0*dx0,axis=3),axis=2)
      result = np.reshape(result,(sout.shape[0],1,1,1))
      var_dict['z_interfaces']=-np.reshape(np.sum(var_dict['dz'],axis=1),(sout.shape[0],1,1,1))
    
      var_dict['z_interfaces']=np.concatenate((result,var_dict['z_interfaces']),axis=1)

      var_dict['z_indices']=[0]
      var_dict['zax_data']=[np.mean(var_dict['zax_data'])]
      var_dict['zbax_data']=[var_dict['zbax_data'][0],var_dict['zbax_data'][-1]]      
      if var_dict['Ztype'] == 'Fixed':
        var_dict['z_interfaces']=np.take(var_dict['z_interfaces'],[0],axis=0)
        
      var_dict['Z']=None
      var_dict['X']=None
      var_dict['xax_data']=[np.mean(var_dict['xax_data'])]
      var_dict['Y']=None
      var_dict['yax_data']=[np.mean(var_dict['yax_data'])]      
      var_dict['rootgrp']=None
      

      self.var_dict[name]=var_dict
      self.variables[name]=name

  def time_avg(self,field=None,vol_weight=True):
    """

    Calculate a finite-volume-weighted average of (field)
    along the time dimension.
    
    """

    if self.var_dict[field]['T'] is None:
      return None
      

    cmd = string.join(['sout=self.',field],sep='')
    exec(cmd)

    var_dict = dict.copy(self.var_dict[field]) # inherit variable dictionary from parent 


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
            w_masked = w
    else:
        w = np.ones((sout.shape))      
        w_masked = w
            
    dt=self.var_dict[field]['dt'][:]

    shape=sout.shape
    w_masked = w_masked*np.tile(np.reshape(dt,(shape[0],1,1,1)),(1,shape[1],shape[2],shape[3]))
    result = np.ma.zeros((np.hstack((1,sout.shape[1:]))))
    nz=shape[1]+1
    result2 = np.ma.zeros((np.hstack((1,nz,shape[2:]))))
    avg_time = 0.
    for i in  np.arange(0,sout.shape[0]):
      result[0,:,:,:]=sout[i,:,:,:]*w_masked[i,:,:,:] + result[0,:,:,:]

      if var_dict['z_interfaces'] is not None:
          if var_dict['Ztype'] is not 'Fixed':
              result2[0,:,:,:]=var_dict['z_interfaces'][i,:,:,:]*dt[i] + result2[0,:,:,:]
          else:
              result2[0,:,:,:]=var_dict['z_interfaces'][:,:,:]*dt[i] + result2[0,:,:,:]
              
      avg_time = avg_time+dt[i]*var_dict['tax_data'][i]

    result = result/np.sum(w_masked,axis=0)

    if var_dict['z_interfaces'] is not None:
      result2 = result2/np.sum(dt)
      
    avg_time = avg_time/np.sum(dt)

    

    dz = np.ma.zeros((np.hstack((1,sout.shape[1:]))))

    if var_dict['z_interfaces'] is not None:
      for k in np.arange(0,sout.shape[1]):
        dz[0,k,:,:]=result2[0,k,:,:]-result2[0,k+1,:,:]

    if var_dict['Z'] is not None:
        if var_dict['Ztype'] is not 'Fixed':        
            var_dict['z_interfaces']=result2
            var_dict['dz']=np.squeeze(dz)
            tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))
            var_dict['z'] = tmp[:,0:-1,:,:]
            var_dict['z'] = np.reshape(var_dict['z'],(1,shape[1],shape[2],shape[3]))
            var_dict['dz'] = np.reshape(var_dict['dz'],(1,shape[1],shape[2],shape[3]))            

        else:
            var_dict['z_interfaces']=np.squeeze(result2)
            var_dict['dz']=np.squeeze(dz)
            tmp = 0.5*(var_dict['z_interfaces']+np.roll(var_dict['z_interfaces'],axis=0,shift=-1))
            var_dict['z'] = tmp[0:-1,:,:]
            var_dict['z'] = np.reshape(var_dict['z'],(shape[1],shape[2],shape[3]))
            var_dict['dz'] = np.reshape(var_dict['dz'],(shape[1],shape[2],shape[3]))                        
    var_dict['tax_data'] = avg_time
    tbax_data = [var_dict['tbax_data'][0],var_dict['tbax_data'][-1]]
    var_dict['tbax_data']=tbax_data
    date_bounds = [var_dict['date_bounds'][0],var_dict['date_bounds'][-1]]
    var_dict['dates'] = num2date(avg_time,var_dict['tunits'],var_dict['calendar'])
    var_dict['T'] = 'time'
    var_dict['dt']=[np.sum(dt)]
    var_dict['rootgrp']=None

    
    

    name = field+'_tav'
    vars(self)[name]=result
    self.var_dict[name]=var_dict
    self.variables[name]=name      

  def monthly_avg(self,field=None,vol_weight=True):
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
        if self.var_dict[field]['Z'] is not None:
            if self.var_dict[field]['Ztype'] is 'Fixed':

                w=self.var_dict[field]['dz'][:]
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
            w = np.ones((sout.shape))      
            w_masked = w
    else:
        w = np.ones((sout.shape))      
        w_masked = w

        
    dt=self.var_dict[field]['dt'][:]

    shape=sout.shape
    result = np.ma.zeros((np.hstack((12,sout.shape[1:]))))
    nz=shape[1]+1
    result2 = np.ma.zeros((np.hstack((12,nz,shape[2:]))))
    num_samples=np.zeros((12))
    months=get_months(self.var_dict[field]['dates'])
    weights=np.zeros((12,shape[1],shape[2],shape[3]))

    for i in  np.arange(0,sout.shape[0]):
      result[months[i]-1,:,:,:]=sout[i,:,:,:]*w_masked[i,:,:,:] + result[months[i]-1,:,:,:]
      weights[months[i]-1,:]=weights[months[i]-1,:]+w_masked[i,:]

      if var_dict['Ztype'] is not 'Fixed' and var_dict['z_interfaces'] is not None:
        result2[months[i]-1,:,:,:]=var_dict['z_interfaces'][i,:,:,:] + result2[months[i]-1,:,:,:]
        
      num_samples[months[i]-1]=num_samples[months[i]-1]+1

    for i in np.arange(0,12):
        result[i,:] = result[i,:]/weights[i,:]

        if var_dict['Ztype'] is not 'Fixed' and var_dict['z_interfaces'] is not None:
            result2[i,:] = result2[i,:]/num_samples[i]
      



    if var_dict['Ztype'] is not 'Fixed' and var_dict['z_interfaces'] is not None:
        dz = np.ma.zeros((np.hstack((12,sout.shape[1:]))))        
        for i in np.arange(0,12):
            for k in np.arange(0,sout.shape[1]):
                dz[i,k,:,:]=result2[i,k,:,:]-result2[i,k+1,:,:]

    if var_dict['Ztype'] is not 'Fixed' and var_dict['Z'] is not None:
      var_dict['z_interfaces']=result2
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
    self.var_dict[name]=var_dict
    self.variables[name]=name      

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

  def remap_Z_to_layers(self,temp_name='temp',salt_name='salt',R=None,p_ref=2.e7,wet=None,nkml=None,nkbl=None,hml=None,fit_target=None):
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

    import vertmap as vmap

    if self.var_dict[temp_name]['Ztype'] is not 'Fixed':
      print """
        Need to provide z-space temperature and salinity in call to remap_Z_to_layers
        """
      return None

    print '...Filling interior points in z-space'
    
    self.fill_interior(temp_name)
    self.fill_interior(salt_name)

    
    self.create_field('wright_eos(vars(self)[\''+temp_name+'\'],vars(self)[\''+salt_name+'\'],2.e7)','sigma2',var_dict=self.var_dict[temp_name])

    
    Rb=np.zeros(len(R)+1)
    Rb[0]=0.0
    Rb[1:-1]=0.5*(R[0:-1]+R[1:])
    Rb[-1]=2.0*Rb[-2]

    rho=self.sigma2.T

    nlevs=self.sigma2.count(axis=1)

    if nlevs.ndim > 2:
        nlevs=nlevs[0,:]
    nlevs=np.squeeze(nlevs.T)
    depth=self.grid.D.T

    zax=self.var_dict[temp_name]['zax_data']
    zbax=self.var_dict[temp_name]['zbax_data']
    
    nt=self.sigma2.shape[0];nz=self.sigma2.shape[1];ny=self.sigma2.shape[2];nx=self.sigma2.shape[3]

    if wet is not None:
      wet_=wet.T
    else:
      wet_=np.ones([nx,ny])

    print '...Finding interface positions '

    zax=zax.astype('float64')
    Rb=Rb.astype('float64')
    depth=depth.astype('float64')
    
    zi=vmap.midas_vertmap.find_interfaces(rho,zax,Rb,depth,nlevs,nkml,nkbl,hml)

    ptemp=vars(self)[temp_name].T
    salt=vars(self)[salt_name].T

    land_fill = -1.e10

    print '...Remapping temperature '
    temp=vmap.midas_vertmap.tracer_z_init(ptemp,-zbax,zi,nkml,nkbl,land_fill,wet_,len(R),nlevs)
    print '...Remapping salinity '
    salt=vmap.midas_vertmap.tracer_z_init(salt,-zbax,zi,nkml,nkbl,land_fill,wet_,len(R),nlevs)

    h=zi-np.roll(zi,axis=2,shift=-1)
    h=h[:,:,0:-1]
    
    if fit_target is not None:
      if fit_target:
        print '...Adjusting temp/salt to fit target densities '
        temp=temp.astype('float64')
        salt=salt.astype('float64')
        R=R.astype('float64')
        h=h.astype('float64')
        vmap.midas_vertmap.determine_temperature(temp,salt,R,p_ref,10,land_fill,h,nkml+nkbl+1)

    temp=np.ma.masked_where(temp==land_fill,temp)
    salt=np.ma.masked_where(salt==land_fill,salt)

    temp=temp.T
    salt=salt.T
    h=h.T
    zi=zi.T
    
    shape=temp.shape

    temp=np.reshape(temp,(1,shape[0],shape[1],shape[2]))
    salt=np.reshape(salt,(1,shape[0],shape[1],shape[2]))
    h=np.reshape(h,(1,shape[0],shape[1],shape[2]))

    zint=np.reshape(zi,(1,shape[0]+1,shape[1],shape[2]))            

        
    zout = 0.5*(zint + np.roll(zint,axis=1,shift=-1))
    zout = zout[:,0:-1,:]

    var_dict=dict.copy(self.var_dict[temp_name])
    var_dict['z_indices']=np.arange(0,len(R))
    var_dict['zb_indices']=np.arange(0,len(Rb))
    var_dict['zax_data']=R
    var_dict['zbax_data']=Rb
    var_dict['z_interfaces']=zint
    var_dict['dz']=h
    var_dict['Z']='potential_density'
    var_dict['z']=zout
    var_dict['zunits']='kg m-3'
    var_dict['Ztype']='Lagrangian'


    self.add_field_from_array(temp,'temp_remap',var_dict=var_dict)
    self.add_field_from_array(salt,'salt_remap',var_dict=var_dict)            


  def adjust_thickness(self,field=None):
    """

    Adjust cell thicknesses based on grid.D 
    
    """

    if self.var_dict[field]['Z'] is None:
      return None
      

    dz = self.var_dict[field]['dz']
    dz_cumsum=np.cumsum(dz,axis=0)

    nz = vars(self)[field].shape[1]    
    D = np.tile(self.grid.D,(nz,1,1))

    floor_dz=D - np.roll(dz_cumsum,shift=1,axis=0)
    mask=np.logical_and(dz_cumsum > D,np.roll(dz_cumsum,shift=1,axis=0) <= D)
    dz[mask]=floor_dz[mask]
    dz_cumsum=np.cumsum(dz,axis=0)
    dz[dz_cumsum>D]=0.0

    z=np.cumsum(dz,axis=0)
    z=0.5*(z+np.roll(z,axis=0,shift=1))
    z[0,:]=0.5*dz[0,:]
    
    self.var_dict[field]['dz']=dz
    self.var_dict[field]['z']=z
    
  def horiz_interp(self,field=None,target=None,src_modulo=False,method='bilinear',PrevState=None,add_NP=True,add_SP=True):
    """
      Interpolate from a spherical grid to a general logically
      rectangular grid using a non-conservative \"bilinear\" interpolation
      algorithm (the default) or \"conservative\" area-weighted.
    """
    
    import hinterp_mod 

    if self.var_dict[field]['Ztype'] is not 'Fixed':
        print """horiz_interp currently only configured for geopotential
              coordinate data """
        return None
    

    open('input.nml','w+')
    
    deg_to_rad=np.pi/180.

    add_np=False
    add_sp=False    

    nj_in = self.grid.lonh.shape[0]; ni_in = self.grid.lonh.shape[0]

    if method=="conservative":
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
    
        if method=="conservative":
            print """Conservative interpolation not available for
                  supergrids"""
            return None

      
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
      
#    lon_in_shifted,xax_shifted = shiftgrid(lon_out.min(),lon_in,xax)

    varin = vars(self)[field].copy()

    if np.ma.is_masked(varin):
        mask_in = np.ma.getmask(varin)
        mask=np.ones((mask_in.shape))
        mask[mask_in]=0.0
    else:
        mask=np.ones((varin.shape))


    nk=varin.shape[1];nt=varin.shape[0]
    
    
    if add_np:
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
      


    missing=-1.e10
    varin=np.ma.filled(varin,missing)


    varout=np.zeros((nt,nk,nj,ni))

    print 'lon_in max/min= ',np.max(lon_in),np.min(lon_in)
    print 'lat_in max/min= ',np.max(lat_in),np.min(lat_in)

    print 'lon_out max/min= ',np.max(lon_out),np.min(lon_out)
    print 'lat_out max/min= ',np.max(lat_out),np.min(lat_out)        

    hinterp_mod.hinterp_mod.hinterp(lon_in.T,lat_in.T,mask.T,varin.T,lon_out.T,lat_out.T,varout.T,src_modulo,method,missing)

        
    varout=np.ma.masked_where(varout==missing,varout)


    
    if PrevState is not None:
      S=PrevState
    else:
      S = state(grid=target)

    
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
            z=var_dict['z'][0,:,0,0]
        elif var_dict['z'].ndim == 3:
            z=var_dict['z'][:,0,0]
        
        z=z.reshape(len(z),1,1)
        z=np.tile(z,(1,nj,ni))
        var_dict['z']=z
        dz=var_dict['dz'][:,0,0]
        dz=dz.reshape(len(dz),1,1)
        dz=np.tile(dz,(1,nj,ni))
        var_dict['dz']=dz    
    
    
    
    S.var_dict[field]=var_dict
    

    return S

  def horiz_interp_refined(self,field=None,target=None,src_modulo=False,method='bilinear',PrevState=None,add_NP=True,add_SP=True):

    deg_to_rad=np.pi/180.

    add_np=False
    add_sp=False    

    nj_in = self.grid.lonh.shape[0]; ni_in = self.grid.lonh.shape[0]

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
            
    if hasattr(target,'x_T_bounds'):  # target grid is a model grid
        nj=target.x_T.shape[0];ni=target.x_T.shape[1]
        lon_out = target.x_T
        lat_out = target.y_T
    elif hasattr(target,'x'): # target grid is a supergrid
        nj=target.x.shape[0];ni=target.x.shape[1]
        lon_out = target.x
        lat_out = target.y

    i_indices = np.arange(0,ni_in).astype(int)
    j_indices = np.arange(0,nj_in).astype(int)
            
    
      
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
              

  def svd(self,fields_left=None,fields_right=None,trunc=1.0):
#
# perform SVD based on multivariate left and right
# side vectors, i.e. Given 2 timeseries A (n1 x m) and
# B (n2 x m), calculate AB^T = UDV^T.
#
# Note: calculate EOFS of A and B then calculate SVD
# of principal components.  This allows for considerable
# memory savings where (n1,n2) > m.

      normalize = 1
      

      if fields_left is None or fields_right is None:
          print """ Required at least 1 field for left and right input vector """
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


  def write_nc(self,filename=None,fields=None,format='NETCDF3_CLASSIC'):
    """

    Write (fields) to an NC file. 

    """  
    if fields is None:
      return None

    if filename is None:
      return None


    f=nc.Dataset(filename,'w',format=format)
    write_interfaces=False
    
    dims=[]
    vars=[]

    
    for field in fields:
      dim_nam=str(self.var_dict[field]['T'])
      if dim_nam not in dims and  dim_nam !=  'None':
          dims.append(dim_nam)
          tdat_=self.var_dict[field]['tax_data']
          if np.isscalar(tdat_):
              tdat=[]
              tdat.append(tdat_)
          else:
              tdat=tdat_
          tdim=f.createDimension(dim_nam,None)
          tv=f.createVariable(dim_nam,'f8',(dim_nam,))
          tv.units= self.var_dict[field]['tunits']
          tv.calendar= self.var_dict[field]['calendar']
          nt=len(tdat)
      dim_nam=str(self.var_dict[field]['Z'])

      if dim_nam not in dims and dim_nam !=  'None':
        dims.append(dim_nam)
        xdat=self.var_dict[field]['zax_data'][:]
        xdim=f.createDimension(dim_nam,len(xdat))
        xv=f.createVariable(dim_nam,'f8',(dim_nam,))
        xv[:]=xdat
        xv.units =   self.var_dict[field]['zunits']
        xv.direction = self.var_dict[field]['Zdir']
        if self.var_dict[field]['Ztype'] is 'Lagrangian' and write_interfaces is False:
          if 'z_interfaces'  in self.var_dict[field]:
              if self.var_dict[field]['z_interfaces'] is not None:
                  write_interfaces = True
                  ifield=field
                  zi=self.var_dict[field]['z_interfaces'][:]
                  ziax=self.var_dict[field]['zbax_data'][:]
      dim_nam=str(self.var_dict[field]['Y'])
      if dim_nam not in dims and dim_nam != 'None':
        dims.append(dim_nam)
        xdat=self.var_dict[field]['yax_data'][:]
        xdim=f.createDimension(dim_nam,len(xdat))
        xv=f.createVariable(dim_nam,'f8',(dim_nam,))
        xv.units =   self.var_dict[field]['yunits']
        xv[:]=xdat
      dim_nam=str(self.var_dict[field]['X'])
      if dim_nam not in dims and dim_nam != 'None':
        dims.append(dim_nam)
        xdat=self.var_dict[field]['xax_data'][:]
        xdim=f.createDimension(dim_nam,len(xdat))
        xv=f.createVariable(dim_nam,'f8',(dim_nam,))
        xv.units =   self.var_dict[field]['xunits']
        xv[:]=xdat

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

      if self.var_dict[field]['FillValue'] is not None:
          FillValue = self.var_dict[field]['FillValue']          
          var=f.createVariable(field,'f4',dimensions=dims,fill_value=FillValue)
      else:
          var=f.createVariable(field,'f4',dimensions=dims)      
      if self.var_dict[field]['missing_value'] is not None:
          var.missing_value = self.var_dict[field]['missing_value']

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
        outv[m][:]=self.__dict__[field][:]
      m=m+1

    m=0;p=0
    for field in fields:
        if self.var_dict[field]['T'] is not None:
            for n in np.arange(0,nt):
                if self.var_dict[field]['X'] is None and self.var_dict[field]['Y'] is None and self.var_dict[field]['Z'] is None:
                    outv[m][n]=self.__dict__[field][n]
                else:
                    outv[m][n,:]=sq(self.__dict__[field][n,:])

                if write_interfaces and m == 0:
                    outv[-1][n,:]=zi[n,:]
        

                tv[n]=tdat[n]
                
            m=m+1
        
    
    f.sync()
    f.close()
    
  def add_figure(self,name,commands):
    """
    
    Add (commands) to write figures from self

    """

    try:
      self.fig_dict[name] = commands
    except:
      
      self.fig_dict={}
      self.fig_dict[name] = commands
      

  def plot_vsect(self,field,title,x=None,y=None,z=None,t=None):

      """

      Make a vertical section of (field) along specified
      coordinate axes.

      """

      xs=0;xe=self.grid.im+1
      ys=0;ye=self.grid.jm+1
      ks=0;ke=len(self.var_dict[field]['zax_data'])+1
      ts=0;te=len(self.var_dict[field]['tax_data'])+1
      
      if x is not None and y is not None:
          xs,xe,ys,ye=find_geo_bounds(self.grid,x=x,y=y)
      else:
          if x is not None:
              xs,xe=find_axis_bounds(self.grid.lonh,x=x)
          if y is not None:
              ys,ye=find_axis_bounds(self.grid.lath,x=y)

      if t is not None:
          ts,te=find_date_bounds(self.var_dict[field]['dates'],t[0],t[1])
          
      if z is not None:
          ks,ke = find_axis_bounds(self.var_dict[field]['zax_data'],x=z)

          
      sout = vars(self)[field][ts:te,ks:ke,ys:ye,xs:xe]
      xax  = self.grid.lonh[xs:xe]
      yax  = self.grid.lath[ys:ye]
      zax  = self.var_dict[field]['zax_data'][ks:ke]

      if len(xax)>1 and len(yax) > 1:
          print """Inconsistent array shape in call to plot_vsect, add section information when calling
                """
          return

      if len(xax)==1:
          xout=yax
      else:
          xout=xax

      plt.pcolormesh(xout,zax,sq(sout))

      plt.show()


  def cfill_hsect(self,field,title,x=None,y=None,z=None,t=None,crange=None,grid=True,axis_labels=True,colormap='jet'):

      """

      Contourfill a section of along specified
      coordinate axes.

      """

      xs=0;xe=self.grid.im
      ys=0;ye=self.grid.jm
      ks=0;ke=len(self.var_dict[field]['zax_data'])
      ts=0;te=len(self.var_dict[field]['tax_data'])
      
      if x is not None and y is not None:
          xs,xe,ys,ye=find_geo_bounds(self.grid,x=x,y=y)
      else:
          if x is not None:
              xs,xe=find_axis_bounds(self.grid.lonh,x=x)
          if y is not None:
              ys,ye=find_axis_bounds(self.grid.lath,x=y)
              
      if t is not None:
          ts,te=find_date_bounds(self.var_dict[field]['dates'],t[0],t[0])
          te=te+1
      if z is not None:
          ks,ke = find_axis_bounds(self.var_dict[field]['zax_data'],x=z)


      sout = vars(self)[field][ts:te,ks:ke,ys:ye,xs:xe]
      xax  = self.grid.lonh[xs:xe]
      yax  = self.grid.lath[ys:ye]
      zax  = self.var_dict[field]['zax_data'][ks:ke]

      if crange is None:
          crange=np.arange(np.min(sout),np.max(sout),(np.max(sout)-np.min(sout))/20.)

      cmap = 'plt.cm.'+colormap

      plt.contourf(xax,yax,sq(sout),levels=crange,cmap=cmap)

      plt.show()
      
      
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



      


      

            




    

  







