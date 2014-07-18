from midas.rectgrid import *
import matplotlib.pyplot as plt
import numpy as np


# Construct a 5 degree lat-lon grid

lon=np.arange(2.5,360.,5.)
lat=np.arange(-77.5,80.,5.)
X,Y=np.meshgrid(lon,lat)
grid=quadmesh(lon=X,lat=Y,cyclic=True)

# Load a single time slice and the first 10 vertical levels of WOA09 data from NODC

S=state(path='http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',fields=['t_an'],z_indices=np.arange(0,10))
S.volume_integral('t_an','XY',normalize=True)
sout=np.squeeze(S.t_an[0,0,:])
notes='Original MEAN/MAX/MIN= %(me)4.2f %(mx)4.2f %(mi)4.2f'%{'me':sq(S.t_an_xyav[0,0]),'mx':sq(np.max(sout)),'mi':sq(np.min(sout))}

T=S.horiz_interp('t_an',target=grid,src_modulo=True,method='bilinear') 
T.volume_integral('t_an','XY',normalize=True)

sout=np.squeeze(T.t_an[0,0,:])
xax=T.grid.x_T
yax=T.grid.y_T

fig,ax=plt.subplots(1)
cf=ax.contourf(xax,yax,sout,np.arange(-2,31,1),extend='both')
plt.colorbar(cf)
ax.grid()
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')


zplt=S.var_dict['t_an']['zax_data'][0]
tit='ANN Temperature from NODC WOA09 at Z= %(zplt)s  meters interpolated to 5 deg'%{'zplt':zplt}

notes2='Regridded MEAN/MAX/MIN= %(me)4.2f %(mx)4.2f %(mi)4.2f'%{'me':sq(T.t_an_xyav[0,0]),'mx':sq(np.max(sout)),'mi':sq(np.min(sout))}

ax.set_title(tit,fontsize=8)

ax.text(0.05,0.95,notes,transform=ax.transAxes,fontsize=6)
ax.text(0.05,0.90,notes2,transform=ax.transAxes,fontsize=6)

plt.show()
