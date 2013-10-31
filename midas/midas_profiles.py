from midas_utils import *
import copy
import matplotlib.pyplot as plt
import scipy as sp


class  profile(object):

  """ A generic data structure for column profile data.  """

  def __init__(self):
    self.data = None
    self.next = None

class profile_list(object):

  def __init__(self,path=None,format='nodc'):

    self.pr=[]

    pr = profile()

    if path is not None:
        if type(path) is list:

          flist=[];prof_dict={}
          for file in path:
            flist.append(nc.Dataset(file))
          i=0
          for f in flist:
            if format == 'nodc':
              pr.data={}
              pr.data['DataType']=f.variables['data_type'][:].tostring()
              pr.data['ref_dateTime']=f.variables['reference_date_time'][:]
              pr.data['wmoid']=f.variables['platform_number'][:].tostring()
              pr.data['cycle_number']=f.variables['cycle_number'][:][0]
              pr.data['direction']=f.variables['direction'][:].tostring()
              pr.data['time']=f.variables['time'][:][0]
              pr.data['latitude']=f.variables['latitude'][:][0]
              pr.data['longitude']=f.variables['longitude'][:][0]

              pr.data['pressure']=sq(f.variables['pressure'][:] )
              pr.data['temp']=sq(f.variables['temperature'][:])
              pr.data['salt']=sq(f.variables['salinity'][:])

              # Avoid confusion with direction attribute
              dp = pr.data['pressure']-np.roll(pr.data['pressure'],shift=-1)
              dp = dp[:-1]
              if np.all(dp>0):
                pr.data['direction']='A'
              elif np.all(dp<0):
                pr.data['direction']='D'
              else:
                print "Problem with profile data. Depths not monotonic"
                raise

              pr_new = copy.copy(pr)

              self.pr.append(pr_new)

  def sort_by_time(self):
    
    self.pr=sorted(self.pr,key=lambda profile: profile.data['time'])
      

  def first_deriv(self,var,bc='interior',positive='down'):

    for pr in self.pr:

      try:
        data=np.squeeze(pr.data[var])
        z=np.squeeze(pr.data['pressure'])
      except:
        print 'unable to find ',var
        raise

      if pr.data['direction'] == 'A':
        dVar = np.roll(data,shift=-1)-data
        Idz  = 1.0/(np.roll(z,shift=-1)-z)
        dVar_dz = dVar*Idz 
        if bc == 'interior':
          dVar_dz[-1]=dVar_dz[-2]
        elif bc == 'zero':
          dVar_dz[-1]=0.


      if pr.data['direction'] == 'D':
        dVar = data-np.roll(data,shift=-1)
        Idz  = 1.0/(z-np.roll(z,shift=-1))
        dVar_dz = dVar*Idz 
        if bc == 'interior':
          dVar_dz[-1]=dVar_dz[-2]
        elif bc == 'zero':
          dVar_dz[-1]=0.
 
      if positive is 'down':
        pr.data['d'+var+'_dZ']=dVar_dz
      else:
        pr.data['d'+var+'_dZ']=-dVar_dz

  def expression(self,expr,name=None,units=None):

    for pr in self.pr:
      cmd = string.join(['pr.data[\'',name,'\']=',expr],sep='')
      exec(cmd)
                   
  def Turner_angle(self):


    # Use positive up convention. Pressure is the vertical
    # coordinate (positive down). 
    
    self.first_deriv(var='temp',bc='interior',positive='up')
    self.first_deriv(var='salt',bc='interior',positive='up')
    
    self.expression('-alpha_wright_eos(pr.data[\'temp\'],pr.data[\'salt\'],pr.data[\'pressure\'])*pr.data[\'dtemp_dZ\']',name='alpha_dTdZ')

    self.expression('beta_wright_eos(pr.data[\'temp\'],pr.data[\'salt\'],pr.data[\'pressure\'])*pr.data[\'dsalt_dZ\']',name='beta_dSdZ')

    self.expression('np.arctan2(pr.data[\'alpha_dTdZ\'] + pr.data[\'beta_dSdZ\'],pr.data[\'alpha_dTdZ\'] - pr.data[\'beta_dSdZ\'])',name='Tu')

  

