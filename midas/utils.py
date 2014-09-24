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


import numpy
from netCDF4 import num2date
import  datetime 
import string
import copy
import types
import pickle

DEBUG = 0


def sq(arr):
    """
    Shorthand for numpy.squeeze()
    """
    res=numpy.squeeze(arr)

    return res
  
#==============================
# Working with Images
#==============================

def image_from_array(arr):
  """
  Convert an array to a greyscale image (experimental).
  """

  try:
    from PIL import Image
  except:
    print "PIL is not installed."
    return None

  arr_min=numpy.min(arr)
  arr=arr+arr_min
  arr_max=numpy.max(arr)
  arr=255*arr/arr_max
  arr=arr.astype('uint8')

  n=arr.shape[0];m=arr.shape[1]
  img=Image.fromarray(arr)

  return img


def array_from_image(img,flipud=False):
  """
  Convert an image to an array (experimental)
  """

  try:
    from PIL import Image
  except:
    print "PIL is not installed."
    return None
 
  from scipy.misc import fromimage

  if flipud:
    img=img.transpose(Image.FLIP_TOP_BOTTOM)
    
  arr=fromimage(img)

  return arr


def store_frame(im,fig,ims):
  """
  Append an image to a pre-existing list
  """

  def setvisible(self,vis):
    for c in self.collections: c.set_visible(vis)


  im.set_visible = types.MethodType(setvisible,im,None)
  im.axes = plt.gca()
  im.figure = fig
  ims.append([im])


def unpickle(file):

  """
 Unpickle a state
  """
    
  S=pickle.load(open(file,'rb'))
  
  return S




def get_axis_cart(dimension,dimname=None):
  """

 **********************************************************************
 Determine axis orientation using a restricted list of detectable
 parameters, returns [X,Y,Z,T]
 **********************************************************************

 
 """


  valid_x_units = ['cm','meters','km','degrees_east','degrees_e','degree_e','deg_e']
  valid_y_units = ['cm','meters','km','degrees_north','degrees_n','degree_n','deg_n']
  valid_z_units = ['cm','meters','km','interface','layer']
  valid_t_units = ['seconds','minutes','hours','days','months','years']

  cart = None

  # first try to retreive the cartesian attribute

  try: 
      cart  = getattr(dimension,'cartesian_axis')
      return cart
  except:
      pass

  try: 
      cart  = getattr(dimension,'axis')
      return cart
  except:
      pass

  try:
      cart=str(dimension.dimensions[0]).upper()
      if cart in ['X','Y','Z','T']:
          return cart
  except:
      pass

  if cart not in ['X','Y','Z','T']:
      cart = None
      
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
      if string.lower(dimname).count('time') >0:          
          cart = 'T'          
          
      
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

  if len(dimension) == 1:
      return dir
  
  if dimension[0] > dimension[1]:
      dir=-1
      
  try:
    orient = getattr(dimension,'positive')
    if string.lower(orient) == 'down':
      dir = -1
    elif orient == -1:
      dir = -1
  except:
    pass

  try:
    orient = getattr(dimension,'direction')
    if orient == -1:
      dir = -1
  except:
    pass
  
  

  return dir



def instance_to_datetime(dates_in):

    fmt='%Y-%m-%d %H:%M:%S'
    dates=[]
    for d in dates_in:
        yr=str(d)[0:4]
        y='%(y)04i'%{'y':int(yr)}
        d=str(d)
        d=y+d[4:]
        dates.append(datetime.datetime.strptime(str(d),fmt))

        
#    dates=[datetime.datetime.strptime(str(d),fmt) for d in dates_in]

    return dates

#def instance_to_datetime(date_in):
#    
#    mon=int(date_in.strftime()[5:7])
#    year=numpy.maximum(int(date_in.strftime()[0:4]),1)
#    day=int(date_in.strftime()[8:10])
#    hr=int(date_in.strftime()[11:13])
#    mn=int(date_in.strftime()[14:16])
#    sec=int(date_in.strftime()[17:19])      
#    date=datetime(year,mon,day,hr,mn,sec)

#    return date

def find_date_bounds(dates_in,tmin,tmax):
  if type(dates_in[0]) is not datetime:
     dates=instance_to_datetime(dates_in)
  else:
     dates = dates_in

  
  ts=-1;te=-1
  for i in numpy.arange(0,numpy.maximum(1,len(dates)-1)):
      if ts == -1 and tmin <= dates[i+1] and tmin >= dates[i]:
          ts = i
      if ts > -1 and tmax <= dates[i+1]:
          te = i+1
          break


  if DEBUG == 1:
      print tmin,tmax,dates[ts],dates[te],ts,te
  
  return ts,te

def time_interp_weights(dates_in,target_in):

    try:
        nt = len(target_in)
    except:
        nt = 1

    if type(dates_in[0]) is not datetime:
        dates=instance_to_datetime(dates_in)
    else:
        dates=dates_in

    if nt > 1:
        if type(target_in[0]) is not datetime:
            target=instance_to_datetime(target_in)
        else:
            target=target_in
    else:
        if type(target_in) is not datetime:
            target=instance_to_datetime(target_in)
        else:
            target=target_in
            
    if nt > 1:
        t1=numpy.zeros(nt,dtype=numpy.int); t2=numpy.zeros(nt,dtype=numpy.int)
        w1=numpy.zeros(nt); w2=numpy.zeros(nt)        
        for i in numpy.arange(0,nt):
            t1[i],t2[i]=find_date_bounds(dates,target[i],target[i])
            date1=dates[t1[i]]
            date2=dates[t2[i]]
            dt=date2-date1;dt=dt.total_seconds()
            dt1=target[i]-date1;dt1=dt1.total_seconds()
            if dt>0:
                w1[i]=1.0-dt1/dt
            else:
                w1[i]=1.0
            w2[i]=1.0-w1[i]
            
    else:
        t1,t2=find_date_bounds(dates,target,target)
        date1=dates[t1];date2=dates[t2]
        dt=date2-date1;dt=dt.total_seconds()
        dt1=target-date1;dt1=dt1.total_seconds()
        if dt>0:
            w1=1.0-dt1/dt
        else:
            w1=1.0
        w2=1.0-w1
        
    return t1,t2,w1,w2

def get_months(dates_in):
  # if type(dates_in[0]) is not datetime:
  #   months = []
  #   for i in numpy.arange(0,len(dates_in)):
  #     mon=int(dates_in[i].strftime()[5:7])
  #     months.append(mon)
  # else:
    months = []
    for i in numpy.arange(0,len(dates_in)):
      months.append(dates_in[i].month)

    return months


def make_monthly_axis(year=1900):
    dates=[];delta=[]
    for i in numpy.arange(1,13):
        dates.append(datetime.datetime(year,i,1))

    dates.append(datetime.datetime(year+1,1,1))
    
    for i in numpy.arange(0,12):
        delta.append((dates[i+1]-dates[i])/2)

    for i in numpy.arange(0,12):
        dates[i]=dates[i]+delta[i]

    return dates[0:12]
        

