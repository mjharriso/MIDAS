from midas.rectgrid import *
import matplotlib.pyplot as plt
import numpy as np

# This is a simple plot routine.


k_plt = 10
grid_obs=quadmesh(path='http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',var='t_an')
S=state(path='http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',grid=grid_obs,geo_region=None,fields=['t_an'],z_indices=np.arange(k_plt,k_plt+1))

xax=S.grid.x_T
yax=S.grid.y_T
sout=np.squeeze(S.t_an[0,0,:])

fig=plt.figure(1,figsize=(5,4))
cf=plt.contourf(xax,yax,sout,np.arange(-2,25,1),extend='both',cmap=plt.cm.jet)
plt.colorbar(cf)
plt.grid()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
zplt=S.var_dict['t_an']['zax_data'][0]
print zplt
tit='ANN Temperature from NODC WOA09 at Z= %(zplt)s  meters '%{'zplt':zplt}
plt.title(tit,fontsize=8)

plt.show()
