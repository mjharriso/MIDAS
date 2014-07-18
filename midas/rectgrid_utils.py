"""
==============================

 This work is licensed under the Creative Commons
 Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 To view a copy of this license, visit   
 http://creativecommons.org/licenses/by-nc-sa/3.0/
 or send a letter to Creative Commons, 444 Castro Street,
 Suite 900, Mountain View, California, 94041, USA.

===============================
"""

import numpy as np


PI_180 = np.pi/180.
R_earth = 6371.e3

def min_resolution(grid=None,x=None,y=None):
  """
  Returns the minimum separation between a point and
  the adjacent point to its left along a 1-dimensional
  axis or 2-d regular grid.

  >>> import rectgrid_utils
  >>> min_dx=rectgrid_utils.min_resolution(x=[0.9999,1.0,2.0])
  >>> print min_dx[0]
  0.0001
  """

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
  """
  Returns the maximum separation between a point and
  the adjacent point to its left along a 1-dimensional
  axis or 2-d regular grid.

  >>> import rectgrid_utils
  >>> max_dx=rectgrid_utils.max_resolution(x=[0.9999,1.0,2.0,10.0])
  >>> print max_dx[0]
  8.0
  """

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
    
  

def find_axis_bounds(axis,x=None,modulo_360=False):
  """
  Returns the bounding indices for axis in the
  range x.

  >>> import rectgrid_utils
  >>> import numpy as np
  >>> axis=np.arange(0.,10.,1.)
  >>> xs,xe=rectgrid_utils.find_axis_bounds(axis,x=[3.,5.])
  >>> print xs,xe
  3 4
  """

  [max_dx,junk] = max_resolution(x=axis)  

  xs=None;xe=None
  if x is not None:
    xmin=x[0];xmax=x[1]
    if modulo_360:
      if x[0]<axis[0]:
        xmin=xmin+360.
      if x[1]>axis[-1]:
        xmax=xmax-360.
    res = np.nonzero(np.abs(axis - xmin) < max_dx)
    xs = res[0][0]
    res = np.nonzero(np.abs(axis - xmax) <= max_dx)
    xe = res[0][0]


  return xs,xe


def cartesian_dist(x1,y1,x2,y2,metric):
    """
    Calculate the distance between (x1,y1) and (x2,y2) on a cartesian grid

    >>> import rectgrid_utils
    >>> x1=0.0;x2=1.0;y1=0.0;y2=1.0
    >>> metric=1.0
    >>> d=rectgrid_utils.cartesian_dist(x1,y1,x2,y2,metric)
    >>> print d
    1.41421356237
    """
    
    dist = metric*np.sqrt((x1-x2)**2.0 + (y1-y2)**2.0)
    
    return dist

def spherical_dist_latlon(x1,y1,x2,y2,metric):
    """
    Calculate the distance between (x1,y1) and (x2,y2) on a sphere along lines
    of constant latitude or longitude

    >>> import rectgrid_utils
    >>> x1=0.0;x2=1.0;y1=0.0;y2=0.0
    >>> metric=R_earth
    >>> d=rectgrid_utils.spherical_dist_latlon(x1,y1,x2,y2,metric)
    >>> print d
    6371000.0
    """


    
    if y1-y2 != 0.:
        dist = metric*np.abs(y1-y2)
    elif x1-x2 != 0.:
        dist = metric*np.cos(y1*PI_180)*np.abs(x1-x2)
    else:
        print """
          This is not a spherical grid"""
        raise
    
    return dist

def mdist(x1,x2):
  """
  Returns positive distance modulo 360.

  >>> import rectgrid_utils
  >>> x1=0.0;x2=730.
  >>> d=rectgrid_utils.mdist(x1,x2)
  >>> print d
  10.0
  """

    
  a=np.mod(x1-x2+720.,360.)
  b=np.mod(x2-x1+720.,360.)

  d=np.minimum(a,b)

  return d

def shiftgrid(lon0,datain,lonsin,start=True,cyclic=360.0):

    import numpy.ma as ma
    """
    Shift global lat/lon grid east or west.
    copied directly from mpl_toolkits v1.0.2 by mjh

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Arguments        Description
    ==============   ====================================================
    lon0             starting longitude for shifted grid
                     (ending longitude if start=False). lon0 must be on
                     input grid (within the range of lonsin).
    datain           original data.
    lonsin           original longitudes.
    ==============   ====================================================

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Keywords         Description
    ==============   ====================================================
    start            if True, lon0 represents the starting longitude
                     of the new grid. if False, lon0 is the ending
                     longitude. Default True.
    cyclic           width of periodic domain (default 360)
    ==============   ====================================================

    returns ``dataout,lonsout`` (data and longitudes on shifted grid).
    """
    if np.fabs(lonsin[-1]-lonsin[0]-cyclic) > 1.e-4:
        # Use all data instead of raise ValueError, 'cyclic point not included'
        start_idx = 0
    else:
        # If cyclic, remove the duplicate point
        start_idx = 1
    if lon0 < lonsin[0] or lon0 > lonsin[-1]:
        msg = 'lon0 outside of range of lonsin %(l0)4.1f %(st)4.1f %(ed)4.1f'%{'l0':lon0,'st':lonsin[0],'ed':lonsin[-1]}
        raise ValueError(msg)
    i0 = np.argmin(np.fabs(lonsin-lon0))
    i0_shift = len(lonsin)-i0
    if ma.isMA(datain):
        dataout  = ma.zeros(datain.shape,datain.dtype)
    else:
        dataout  = np.zeros(datain.shape,datain.dtype)
    if ma.isMA(lonsin):
        lonsout = ma.zeros(lonsin.shape,lonsin.dtype)
    else:
        lonsout = np.zeros(lonsin.shape,lonsin.dtype)
    if start:
        lonsout[0:i0_shift] = lonsin[i0:]
    else:
        lonsout[0:i0_shift] = lonsin[i0:]-cyclic
    dataout[:,0:i0_shift] = datain[:,i0:]
    if start:
        lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]+cyclic
    else:
        lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]
    dataout[:,i0_shift:] = datain[:,start_idx:i0+start_idx]
    return dataout,lonsout

def addcyclic(arrin,lonsin):
    """
    ``arrout, lonsout = addcyclic(arrin, lonsin)``
    adds cyclic (wraparound) point in longitude to ``arrin`` and ``lonsin``.
    """
    nlats = arrin.shape[0]
    nlons = arrin.shape[1]
    if ma.isMA(arrin):
        arrout  = ma.zeros((nlats,nlons+1),arrin.dtype)
    else:
        arrout  = np.zeros((nlats,nlons+1),arrin.dtype)
    arrout[:,0:nlons] = arrin[:,:]
    arrout[:,nlons] = arrin[:,0]
    if ma.isMA(lonsin):
        lonsout = ma.zeros(nlons+1,lonsin.dtype)
    else:
        lonsout = np.zeros(nlons+1,lonsin.dtype)
    lonsout[0:nlons] = lonsin[:]
    lonsout[nlons]  = lonsin[-1] + lonsin[1]-lonsin[0]
    return arrout,lonsout

if __name__ == "__main__":
    import doctest
    doctest.testmod()
