import numpy as np
import netCDF4 as nc
import copy

PI_180 = np.pi/180.

def min_resolution(grid=None,x=None,y=None):
  if grid is not None:
    dx = grid.grid_x-np.roll(grid.grid_x,1)
    dy = grid.grid_y-np.roll(grid.grid_y,1)
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
    dx = grid.grid_x-np.roll(grid.grid_x,1)
    dy = grid.grid_y-np.roll(grid.grid_y,1)
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
    xs = np.nonzero(np.abs(grid.grid_x - xmin) <= max_dx)[0][0]
    xe = np.nonzero(np.abs(grid.grid_x - xmax) <= max_dx)[0][0]

  ys=None;ye=None
  if y is not None:
    ymin=y[0];ymax=y[1]
    ys = np.nonzero(np.abs(grid.grid_y - ymin) <= max_dy)[0][0]
    ye = np.nonzero(np.abs(grid.grid_y - ymax) <= max_dy)[0][0]

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

def cartesian_dist(x1,y1,x2,y2,metric):
    
# Calculate the distance between (x1,y1) and (x2,y2) on a cartesian grid
    
    dist = metric*np.sqrt((x1-x2)**2.0 + (y1-y2)**2.0)
    
    return dist

def spherical_dist(x1,y1,x2,y2,metric):

# Calculate the distance between (x1,y1) and (x2,y2) on a sphere along lines
# of constant latitude or longitude
    
    if np.all(y1-y2):
        dist = metric*np.abs(y1-y2)
    elif np.all(x1-x2):
        dist = metric*np.cos(y1*PI_180)*np.abs(x1-x2)
    else:
        print """
          This is not a spherical grid"""
        raise
    
    return dist

def mdist(x1,x2):

# Returns positive distance modulo 360
    
    a=np.mod(x1-x2+720.,360.)
    b=np.mod(x2-x1+720.,360.)

    d=np.minimum(a,b)

    return d

class supergrid(object):

  def __init__(self,nxtot=0,nytot=0,config=None,axis_units=None,ystart=0.0,leny=0.0,xstart=0.0,lenx=0.0,xdat=None,ydat=None,cyclic_x=False,cyclic_y=False,tripolar_n=False,displace_pole=False,r0_pole=0.0,lon0_pole=0.0,doughnut=0.0,radius=6.378e6,file=None):
      """

 A super-grid for FMS-based models

"""
      
      self.yDir = 1
      
      if file is not None:
        f=nc.Dataset(file)
        vdict={}
        self.x=f.variables['x'][:]
        self.y=f.variables['y'][:]
        try:
          self.dx=f.variables['dx'][:]
          self.have_metrics=True
        except:
          self.have_metrics=False
        if self.have_metrics:
          self.dy=f.variables['dy'][:]
          self.area=f.variables['area'][:]
          self.angle_dx=f.variables['angle_dx'][:]
        self.grid_x=self.x[0,:]
        self.grid_y=self.y[:,0]
        units=getattr(f.variables['x'],'units')
        if units == 'meters':
          self.is_cartesian=True
          self.is_latlon=False          
        else:
          self.is_cartesian=False
          self.is_latlon=True
        self.im=self.grid_x.shape[0]
        self.jm=self.grid_y.shape[0]        
        return 
      
      if xdat is not None:
        if ydat is not None:
          vdict={}
          vdict['nxtot']=xdat.shape[1]
          vdict['nytot']=ydat.shape[0]
          vdict['config']=config
          vdict['axis_units']=axis_units
          vdict['ystart']=ydat.min()
          vdict['leny']=ydat.max()-ydat.min()
          vdict['xstart']=xdat.min()
          vdict['lenx']=xdat.max()-xdat.min()
          vdict['cyclic_x']=cyclic_x
          vdict['cyclic_y']=cyclic_y
          vdict['tripolar_n']=tripolar_n
          vdict['displaced_pole']=displace_pole
          vdict['r0_pole']=r0_pole
          vdict['lon0_pole']=lon0_pole
          vdict['doughnut']=doughnut
          vdict['radius']=radius

          jind=np.arange(vdict['nytot']);iind=np.arange(vdict['nxtot'])
          jindp=np.arange(vdict['nytot']+1);iindp=np.arange(vdict['nxtot']+1)                  
          self.grid_y=vdict['ystart']+jindp*vdict['leny']/vdict['nytot']
          self.grid_x=vdict['xstart']+iindp*vdict['lenx']/vdict['nxtot']          
          self.grid_x=vdict['xstart']+iindp*vdict['lenx']/vdict['nxtot']
          self.x=xdat
          self.y=ydat
          self.is_cartesian=False
          self.is_latlon=False
          if vdict['config'] == 'cartesian':
            self.is_cartesian = True
          else:
            self.is_latlon = True

          self.im=self.grid_x.shape[0]-1
          self.jm=self.grid_y.shape[0]-1

          self.dict=dict.copy(vdict)
          
          self.have_metrics = False
          
          return
        
        else:
          
          print """ Both xdat,ydat need to be defined if one is"""
          return None
        
      vdict={}
      vdict['nxtot']=nxtot
      vdict['nytot']=nytot
      vdict['config']=config
      vdict['axis_units']=axis_units
      vdict['ystart']=ystart
      vdict['leny']=leny
      vdict['xstart']=xstart
      vdict['lenx']=lenx
      vdict['cyclic_x']=cyclic_x
      vdict['cyclic_y']=cyclic_y
      vdict['tripolar_n']=tripolar_n
      vdict['radius']=radius
      vdict['displaced_pole']=displace_pole
      vdict['r0_pole']=r0_pole
      vdict['lon0_pole']=lon0_pole
      vdict['doughnut']=doughnut
      
      if config=='mercator':
        vdict['isotropic']=True
      else:
        vdict['isotropic']=False
        
      self.dict=dict.copy(vdict)
        
      jind=np.arange(nytot);iind=np.arange(nxtot)
      jindp=np.arange(nytot+1);iindp=np.arange(nxtot+1)        

      if config == 'cartesian':
        self.is_cartesian=True
        self.is_latlon=False        
        self.grid_y=ystart+jindp*leny/nytot
        self.grid_x=xstart+iindp*lenx/nxtot
        self.x=np.tile(self.grid_x,(nytot+1,1))
        self.y=np.tile(self.grid_y.reshape((nytot+1,1)),(1,nxtot+1))
      elif config == 'spherical':
        self.is_cartesian=False
        self.is_latlon=True
        self.grid_y=ystart+jindp*leny/nytot
        self.grid_x=xstart+iindp*lenx/nxtot
        self.x=np.tile(self.grid_x,(nytot+1,1))
        self.y=np.tile(self.grid_y.reshape((nytot+1,1)),(1,nxtot+1))
      elif config == 'mercator':
        self.is_cartesian=False
        self.is_latlon=True             
        jRef=np.floor(-nytot*ystart/leny)
        fnRef=self.Int_dj_dy(0.0)
        y0=ystart*PI_180
        self.grid_y=np.zeros(nytot+1)
        itt=0
        for j in np.arange(nytot+1):
          jd = fnRef + (j-jRef)
          self.grid_y[j]=self.find_root_y(jd,y0,-2.0*np.pi,2.0*np.pi,itt)
          y0=self.grid_y[j]
        self.grid_x=xstart+iindp*lenx/nxtot
        self.grid_y=self.grid_y/PI_180
        self.x=np.tile(self.grid_x,(nytot+1,1))
        self.y=np.tile(self.grid_y.reshape((nytot+1,1)),(1,nxtot+1))
      if tripolar_n:
        self.lon_bpnp=self.grid_x[0]
        self.join_lat=self.grid_y[0]
        self.rp=np.tan(0.5*(0.5*np.pi - (self.grid_y[0])*PI_180))            
        lon,lat=self.tp_trans()
        self.x=lon
        self.y=lat
      if displace_pole:
        r,phi = self.displaced_pole(r0_pole,lon0_pole,excluded_fraction=doughnut)
        self.x=phi.copy()
        self.y=r.copy()
        self.grid_x = self.x[-1,:]
        self.grid_x[-1]=self.grid_x[-1]+360.0
        self.grid_y=self.y[:,nxtot/4]
        
      self.im=self.grid_x.shape[0]-1
      self.jm=self.grid_y.shape[0]-1
        
      self.have_metrics = False
        
  def dy_dj(self,y):
# Returns the grid increment in the j direction given y in radians
        
    nxtot=self.dict['nxtot'];nytot=self.dict['nytot']
    lenx=self.dict['lenx'];leny=self.dict['leny']

    if self.dict['isotropic']:
      C0=PI_180*lenx/nxtot
      dydj=C0*np.cos(y)
    else:
      dydj=PI_180*leny/nytot
      
    return dydj

  def ds_dj(self,x,y):
# Returns the grid spacing in the logical y direction given (x,y) in radians
    dsdj=self.dict['radius']*self.dy_dj(y)

        

  def dx_di(self,x):
# Returns the grid spacing in the logical x direction given (x,y) in radians
        
    dxdi=PI_180*self.dict['lenx']/self.dict['nxtot']
    
    return dxdi

  def ds_di(self,x,y):
# Returns the grid spacing in the logical x direction given (x,y) in radians
    dsdi=self.dict['radius']*np.cos(y)*self.dx_di(x)

    return dsdi

  def dL(self,x1,x2,y1,y2):
        
# Calculate the contribution from the line integral along one of the four sides of a cell
# face to the area of a cell, assuming that the sides follow a linear path in latitude and
# longitude.
        
    dy=y2-y1

    if np.abs(dy) > 2.5e-8:
      r=((1.0-np.cos(dy))*np.cos(y1) + np.sin(dy)*np.sin(y1))/dy
    else:
      r=(0.5*dy*np.cos(y1) + np.sin(y1))

    dl=r*(x2-x1)

    return dl

    
    

  def Int_di_dx(self,x):
# Calculates and returns the integral of the inverse
# of dx/di to the point x, in radians.
    Int_didx=x*self.dict['nxtot']/self.dict['lenx']/PI_180

    return Int_didx

  def Int_dj_dy(self,y):
# Calculates and returns the integral of the inverse
# of dy/dj to the point y, in radians.        
    if self.dict['isotropic']:
      I_C0=self.dict['nxtot']/self.dict['lenx']/PI_180


      if y>=0.0:
        r=I_C0*np.log((1.0+np.sin(y))/np.cos(y))
      else:
        r=-1.0*I_C0*np.log((1.0-np.sin(y))/np.cos(y))
    else:
      I_C0=self.dict['nytot']/self.dict['leny']/PI_180
      r=I_C0*y
            
    return r
              
                                           
  def find_root_y(self,fnval,y1,ymin,ymax,ittmax):
# Finds and returns the value of y at which the
# monotonic function self.Int_dj_dy takes the value fnval, also returning
# in ittmax the number of iterations of Newton's method that were
# used to polish the root.
    y=y1;ybot=y1
    fnbot=self.Int_dj_dy(ybot)-fnval
    itt=0

#        print 'fnbot=',fnbot

    while fnbot > 0.0:
      if (ybot-2.0*self.dy_dj(ybot)) < 0.5*(ybot+ymin):
        ybot=ybot-2.0*self.dy_dj(ybot)
      else:
        ybot=0.5*(ybot+ymin)
        itt=itt+1
      fnbot=self.Int_dj_dy(ybot) - fnval

      if np.logical_and(itt > 50,fnbot>0.0):
        print """
              Unable to find bottom bound for grid function"""
        raise
        
    if y+2.0*self.dy_dj(y) < 0.5*(y+ymax):
      ytop  = y + 2.0*self.dy_dj(y)
    else:
      ytop = 0.5*(y+ymax)

    fntop = self.Int_dj_dy(ytop) - fnval
    itt=0
    while fntop < 0.0:
      if ytop + 2.0*self.dy_dj(ytop) < 0.5*(ytop+ymax):
        ytop = ytop + 2.0*self.dy_dj(ytop)
      else:
        ytop = 0.5*(ytop+ymax)
        itt=itt+1
      fntop=self.Int_dj_dy(ytop)-fnval

      if np.logical_and(itt>50,fntop<0.0):
        print """
              Unable to find top bound for grid function"""
        raise

    for itt in np.arange(10):
      y=0.5*(ybot+ytop)
      fny=self.Int_dj_dy(y)-fnval
      if fny < 0.0:
        fnbot=fny
        ybot=y
      else:
        fntop=fny
        ytop=y

    for itt in np.arange(10):
      dy_dfn=self.dy_dj(y)
      fny=self.Int_dj_dy(y)-fnval
      dy=-1.0*fny*dy_dfn
      y=y+dy
      if y>ytop:
        y=ytop
      if y<ybot:
        y=ybot
      if np.abs(dy) < 8.0e-15*np.abs(y)+1.e-20:
        break

    if np.abs(y) < 1.e-12:
      y=0.0

    ittmax=itt

    return y

  def bp_colat(self):


    d=mdist(self.x,self.lon_bpnp)
    phi=d*PI_180


    return phi
    
  def bp_lon(self):
        
# Co-latitudes of rotated grid sphere
    chi=2.0*np.arctan(np.tan(0.5*(0.5*np.pi-self.y*PI_180))/self.rp)
# Base grid longitude
    lam = 0.5*np.pi - chi
    nxtot=self.dict['nxtot']
    lam[:,:nxtot/2]=lam[:,:nxtot/2]-np.pi/2
    lam[:,nxtot/2:]=np.pi/2-lam[:,nxtot/2:]
    
    return lam

  def tp_trans(self):
    
    nytot=self.dict['nytot']
    nxtot=self.dict['nxtot']        
        
    lamc=self.bp_lon()
    phic=self.bp_colat()
        
    chic=np.arccos(np.sin(phic)*np.cos(lamc))
    gamma=np.arctan(np.tan(phic)*np.sin(lamc))        
    lat=2.0*np.arctan(self.rp*np.tan(chic/2.0))
    lat=lat/PI_180
    lat=90.-lat
    lon=gamma/PI_180

    lon[:,:nxtot/4]=-lon[:,:nxtot/4]
    lon[:,nxtot/4]=90.0
    lon[:,nxtot/4+1:nxtot/2]=180.0-lon[:,nxtot/4+1:nxtot/2]
    lon[:,nxtot/2]=180.
    lon[:,nxtot/2+1:3*nxtot/4]=180.0-lon[:,nxtot/2+1:3*nxtot/4]
    lon[:,3*nxtot/4]=270.
    lon[:,3*nxtot/4+1:]=360.-lon[:,3*nxtot/4+1:]

    lon=lon+self.lon_bpnp
        
    return lon,lat

  def displaced_pole(self,ra2,phia,excluded_fraction=None,pole=-1,):

    if pole == -1:
      radius = 90.+self.grid_y[-1]
      r=(90.0+self.y)/radius
      print 'ending latitude = ',self.grid_y[-1]      
    else:  # Did not test this option yet (probably does not work)
      radius = 90.-self.grid_y[0]
      r=(90.0-self.y)/radius
      print 'ending latitude = ',self.grid_y[0]
      
    ra=ra2

    print 'applying a conformal remapping of the pole, original  radius = ',radius, ' degrees'
    print 'displaced pole location (relative to unit sphere) = ',ra
    print 'displaced pole angle ( clockwise degrees relative to Greenwich) = ', -phia+180.
    if excluded_fraction is not None:
      print 'excluding inner ',excluded_fraction*100.,' percent of the grid'

    a=np.complex(ra2*np.cos(PI_180*-phia),ra2*np.sin(PI_180*-phia))

    phi=self.x
    
    lon0 = phi[0,0]

    ny,nx=self.x.shape

    theta_chg=np.arccos(2.0*ra/(1.0+ra*ra))
    r_out=np.zeros(phi.shape)
    phi_out=np.zeros(phi.shape)

    for j in np.arange(phi.shape[0]):
      for i in np.arange(phi.shape[1]):

        th_1 = (phi[j,i]+phia)*PI_180
        if th_1 > np.pi:
          th_1=th_1-2.0*np.pi
        elif th_1 < -np.pi:
          th_1 = th_1+2.0*np.pi

        theta=th_1.copy()
        if np.abs(np.sin(th_1)) <= 0.5*np.abs(np.cos(th_1)):
          tan_tgt = np.tan(th_1)
          if th_1 >= 0.0:
            if tan_tgt >= 0.0:
              th_min = 0.0; th_max = theta_chg
            else:
              th_min = theta_chg ; th_max = np.pi
          else:
            if tan_tgt <= 0.0:
              th_max = 0.0 ; th_min = -theta_chg
            else:
              th_max = -theta_chg ; th_min = -np.pi
          if np.logical_or(theta>th_max,theta<th_min):
            theta = 0.5*(th_max+th_min)

          for ii in np.arange(20):
            denom =  np.cos(theta)*(1.+ra*ra) - 2.0*ra;
            val = np.sin(theta)*(1.-ra*ra) / denom;
            dval_dth = ((1.-ra**4) - 2.0*ra*(1.-ra**2)*np.cos(theta) ) / denom**2;
            err = val - tan_tgt;
            if np.abs(err) < 1e-12:
              break
            theta_prev = theta.copy();
            theta = theta - err / dval_dth;
            if theta > th_max:
              theta = 0.5*(theta_prev+th_max)
            if theta < th_min:
              theta = 0.5*(theta_prev+th_min)
        else:
          if th_1==0.0:
            cot_tgt=0.0
          else:
            cot_tgt = 1.0/np.tan(th_1)
                    
          if th_1>=0.0:
            if cot_tgt >=0.0:
              th_min = 0.0; th_max = theta_chg.copy()
            else:
              th_min = theta_chg.copy() ; th_max = np.pi
          else:
            if cot_tgt<=0.0:
              th_max = 0.0
              th_min = -theta_chg
            else:
              th_max = -theta_chg
              th_min = -np.pi
          if np.logical_or(theta>th_max,theta<th_min):
            theta = 0.5*(th_max+th_min)


          for ii in np.arange(20):
            denom = np.sin(theta)*(1.0-ra*ra)
            val=(np.cos(theta)*(1.0+ra*ra) - 2.0*ra)/denom
            dval_dth = (-(1.0-ra**4.0) + 2.0*ra*(1.0-ra**2.0)*np.cos(theta)) / denom**2.
            err = val -  cot_tgt
            if np.abs(err) < 1.e-12:
              break
            theta_prev=theta.copy()
            theta = theta - err / dval_dth
            if theta>th_max:
              theta=0.5*(theta_prev+th_max)
            if theta<th_min:
              theta=0.5*(theta_prev+th_min)

        th_cor = theta-phia*PI_180
        z=np.complex(r[j,i]*np.cos(th_cor),r[j,i]*np.sin(th_cor))
        w = (z-a)/(1.0-z*np.conj(a))
        r_out[j,i]=np.abs(w)
        phi_out[j,i]=np.angle(w)/PI_180

    r_out=-90.0+r_out*radius
    phi_out[phi_out<lon0]=phi_out[phi_out<lon0]+360.0
    phi_out[phi_out>=lon0+360.]=phi_out[phi_out>=lon0+360.]-360.0
    phi_out[np.abs(phi_out[:,:-1]-(lon0+360.0))<1.e-4]=lon0
    
    jmin=0;jmax=ny
    if excluded_fraction is not None:
      if pole == -1:
        jmin=np.ceil(ny*excluded_fraction)
        jmin=jmin-np.mod(jmin,2)
        return r_out[jmin:,:], phi_out[jmin:,:]        
      else:
        jmax=np.ceil(ny*(1.0 - excluded_fraction))
        return r_out[:jmax,:], phi_out[:jmax,:]                
        



        
        
  def merge_grid(self,grid2):

    self.x=np.concatenate((self.x,grid2.x),axis=0)
    self.y=np.concatenate((self.y,grid2.y),axis=0)
    self.dx=np.concatenate((self.dx,grid2.dx),axis=0)                
    self.dy=np.concatenate((self.dy,grid2.dy),axis=0)
    self.area=np.concatenate((self.area,grid2.area),axis=0)
    self.angle_dx=np.concatenate((self.angle_dx,np.take(grid2.angle_dx,[0],axis=0)),axis=0)    
    self.angle_dx=np.concatenate((self.angle_dx,grid2.angle_dx),axis=0)
    

  def create_grid(refine=1):

    from midas import generic_grid

    grid=generic_grid(supergrid=self)
      
    
  def grid_metrics(self):


    nytot=self.dict['nytot']
    nxtot=self.dict['nxtot']
        
    if  self.dict['axis_units'] == 'm':
      metric=1.0
    if  self.dict['axis_units'] == 'km':            
      metric=1.e3
    if  self.dict['axis_units'] == 'degrees':                        
      metric=self.dict['radius']*PI_180


    Re = self.dict['radius']
    ymid_j = 0.5*(self.y+np.roll(self.y,shift=-1,axis=0))
    ymid_i = 0.5*(self.y+np.roll(self.y,shift=-1,axis=1))      
    dy_j = np.roll(self.y,shift=-1,axis=0) - self.y
    dy_i = np.roll(self.y,shift=-1,axis=1) - self.y
    dx_i = mdist(np.roll(self.x,shift=-1,axis=1),self.x)
    dx_j = mdist(np.roll(self.x,shift=-1,axis=0),self.x)
    self.dx = metric*metric*(dy_i*dy_i + dx_i*dx_i*np.cos(ymid_i*PI_180)*np.cos(ymid_i*PI_180))
    self.dx = np.sqrt(self.dx)
    self.dy = metric*metric*(dy_j*dy_j + dx_j*dx_j*np.cos(ymid_j*PI_180)*np.cos(ymid_j*PI_180))
    self.dy = np.sqrt(self.dy)
    self.dx=self.dx[:,:-1]
    self.dy=self.dy[:-1,:]

    self.area=self.dx[:-1,:]*self.dy[:,:-1]
    self.angle_dx=np.zeros((nytot,nxtot))

    self.angle_dx = np.arctan2(dy_i,dx_i)*180.0/np.pi
#    self.angle_dx = self.angle_dx[:-1,:-1]
                  
      
      
    self.have_metrics = True


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
            
    

  def extract(self,geo_region=None):
    if geo_region is None:
      grid=copy.copy(self)
      return grid
    else:
      print """Supergrid extraction not supported yet """
      return None
      
  def write_nc(self,fnam=None,format='NETCDF3_CLASSIC'):
    import netCDF4 as nc

    if fnam is None:
      fnam='supergrid.nc'

    f=nc.Dataset(fnam,'w',format=format)

    dims=[]
    vars=[]

    nyp=f.createDimension('nyp',len(self.grid_y))
    nxp=f.createDimension('nxp',len(self.grid_x))
    ny=f.createDimension('ny',len(self.grid_y)-1)
    nx=f.createDimension('nx',len(self.grid_x)-1)    
    
    yv=f.createVariable('y','f8',('nyp','nxp'))
    xv=f.createVariable('x','f8',('nyp','nxp'))    
    yv.units=self.dict['axis_units']
    xv.units=self.dict['axis_units']    
    yv[:]=self.y
    xv[:]=self.x

    if self.have_metrics:
      dyv=f.createVariable('dy','f8',('ny','nxp'))
      dyv.units='meters'
      dyv[:]=self.dy
      dxv=f.createVariable('dx','f8',('nyp','nx'))
      dxv.units='meters'
      dxv[:]=self.dx
      areav=f.createVariable('area','f8',('ny','nx'))
      areav.units='m2'
      areav[:]=self.area
      anglev=f.createVariable('angle_dx','f8',('nyp','nxp'))
      anglev.units='degrees'
      anglev[:]=self.angle_dx            

    f.sync()
    f.close()
      

        
