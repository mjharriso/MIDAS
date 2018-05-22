from midas.rectgrid import *
from midas.profiles import *
import matplotlib.pyplot as plt
import numpy as np
import datetime

# This is a simple plot routine.


fpath='https://opendap.jpl.nasa.gov:443/opendap/OceanTemperature/amsre/L3/sst_1deg_1mo/tos_AMSRE_L3_v7_200206-201012.nc'
grid_obs=quadmesh(path=fpath,var='tos')
S=state(path=fpath,grid=grid_obs,geo_region=None,fields=['tos'],date_bounds=[datetime.datetime(2005,8,14),datetime.datetime(2005,8,14)],default_calendar='noleap')
P=state(path=fpath,grid=grid_obs,geo_region=None,fields=['tos'],date_bounds=[datetime.datetime(2003,8,14),datetime.datetime(2003,8,14)],default_calendar='noleap')

xax=S.grid.x_T
yax=S.grid.y_T
sout=np.squeeze(S.tos[0,0,:])
pout=np.squeeze(P.tos[0,0,:])

fig=plt.figure(1,figsize=(5,4))
cf=plt.contourf(xax,yax,sout,np.arange(273,303,.5),extend='both',cmap=plt.cm.spectral)
plt.colorbar(cf)
plt.grid()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
tit='AMSRE L3 SST '+str(S.var_dict['tos']['dates'][0])
plt.title(tit,fontsize=8)

fig2=plt.figure(2,figsize=(5,4))
cf=plt.contourf(xax,yax,sout-pout,np.arange(-5,5,.25),extend='both',cmap=plt.cm.bwr)
plt.colorbar(cf)
plt.grid()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
tit='AMSRE L3 SST '+str(S.var_dict['tos']['dates'][0])+' minus '+str(P.var_dict['tos']['dates'][0])
plt.title(tit,fontsize=8)

plt.show()
