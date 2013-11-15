from midas_utils import *
from wright_eos import *
import copy
import matplotlib.pyplot as plt
import scipy as sp
import netCDF4 as nc
from dateutil import parser

try:
  import gsw
except:
  pass

class profile_index(object):

  def __init__(self,index_file=None,lon_bounds=None,lat_bounds=None,months=None,index_file_format='nodc'):

    import os

    if index_file is None:
      return None

    f_list = []

    full_list = os.listdir(index_file)
    for f in full_list:
      if f.find('.txt') > 0:
        f_list.append(f)

    f_list.sort()

    self.f_list = f_list

    self.dict = []

    for f in f_list:
      fo=open(f)
      line=fo.readline()
      while True:
        dict={}
        line=fo.readline()
        line=line.split(',')
        if len(line) >= 14: 
          dict['id']=line[0]
          dict['path']=line[2]
          dict['date_string']=line[4]
          dict['time_string']=line[5]
          dict['ncDateTime']=parser.parse(dict['date_string'])
          dict['lat']=np.float(line[7])
          dict['lon']=np.float(line[8])
          dict['data_mode']=line[11]
          dict['num_lev']=np.int(line[12])


          passed=0;failed=0

          if lon_bounds is not None or lat_bounds is not None:
            if lon_bounds is not None:
              if np.logical_and((dict['lon']>=lon_bounds[0]),(dict['lon']<=lon_bounds[1])):
                passed=passed+1
              else:
                failed=failed+1
            if lat_bounds is not None:
              is_n = dict['lat'] >= lat_bounds[0]
              is_s = lat_bounds[1] >= dict['lat'] 

              if np.logical_and(is_n,is_s):
                passed=passed+1
              else:
                failed=failed+1

          if months is not None:
            m=dict['ncDateTime'].month - 1
            if m in months:
              passed=passed+1
            else:
              failed=failed+1
            
          if np.logical_and(passed>=0,failed==0):              
            self.dict.append(dict)
        else:
          break
          

        
        
  def write_ascii(self,file='profile_index_data.asc'):

    f=open(file,'w+')

    for L in self.dict:
      f.write('===============================\n')
      for key,val in L.items():
        txt= '%(k)s = %(v)s \n'%{'k':key,'v':val}
        f.write(txt)

    
    f.close()


    
class profile(object):

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
            f=nc.Dataset(file['path'])
            date_string = file['date_string']
            i=0
            if format == 'nodc':
              pr.data={}
              pr.data['DataType']=f.variables['data_type'][:].tostring()
              pr.data['ref_dateTime']=f.variables['reference_date_time'][:].tostring()
              pr.data['time']=f.variables['time'][:][0]
              reft = pr.data['ref_dateTime']
              reft = 'days since '+reft[0:4]+'-'+reft[4:6]+'-'+reft[6:9]
              pr.data['ncDateTime_indx']=parser.parse(date_string)
              pr.data['ncDateTime']=num2date(pr.data['time'],reft)
              pr.data['wmoid']=f.variables['platform_number'][:].tostring()
              pr.data['cycle_number']=f.variables['cycle_number'][:][0]
              pr.data['direction']=f.variables['direction'][:].tostring()

              pr.data['latitude']=f.variables['latitude'][:][0]
              pr.data['longitude']=f.variables['longitude'][:][0]

              pr.data['pressure']=sq(f.variables['pressure'][:] )
              pr.data['temp']=sq(f.variables['temperature'][:])
              pr.data['salt']=sq(f.variables['salinity'][:])
              pr.data['positioning_system']=string.join(f.variables['positioning_system'][:])
              pr.data['position_qc']=string.join(f.variables['position_qc'][:])
              pr.data['inst_reference']=string.join(f.variables['inst_reference'][:])

                # Avoid confusion with direction attribute
              dp = np.roll(pr.data['pressure'],shift=-1)-pr.data['pressure']
              dp = dp[:-2]

              if len(dp) > 1:
                if np.all(dp>0):
                  pr.data['direction']='A'
                  pr.data['min_dp']=np.min(dp)
                elif np.all(dp<0):
                  pr.data['direction']='D'
                  pr.data['min_dp']=np.min(-dp)
                else:
                  continue
              else:
                continue
              pr_new = copy.copy(pr)
              
              self.pr.append(pr_new)

      

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
                   
  def sort_by_time(self):
    
    self.pr=sorted(self.pr,key=lambda profile: profile.data['time'])

  def sort_by_id(self):
    
    self.pr=sorted(self.pr,key=lambda profile: profile.data['wmoid'])

  def show_profile_ids(self):

    print 'UNIQUE WMOID LIST'
    print '================='

    prev_id = ''
    for pr in self.pr:
      if pr.data['wmoid'] != prev_id:
        print pr.data['wmoid']
        prev_id = pr.data['wmoid'] 

        
  

  def Turner_angle(self):


    # Use positive up convention. Pressure is the vertical
    # coordinate (positive down). 
    
    self.first_deriv(var='temp',bc='interior',positive='up')
    self.first_deriv(var='salt',bc='interior',positive='up')
    
    self.expression('-alpha_wright_eos(pr.data[\'temp\'],pr.data[\'salt\'],pr.data[\'pressure\'])*pr.data[\'dtemp_dZ\']',name='alpha_dTdZ')

    self.expression('beta_wright_eos(pr.data[\'temp\'],pr.data[\'salt\'],pr.data[\'pressure\'])*pr.data[\'dsalt_dZ\']',name='beta_dSdZ')

    self.expression('np.arctan2(pr.data[\'alpha_dTdZ\'] + pr.data[\'beta_dSdZ\'],pr.data[\'alpha_dTdZ\'] - pr.data[\'beta_dSdZ\'])',name='Tu')

  
  
