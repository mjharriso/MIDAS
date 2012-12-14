from midas import *


grid_obs=generic_rectgrid(path='http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',var='t_an')
S=state(path='http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',grid=grid_obs,geo_region=None,fields=['t_an'])

xax=S.grid.x_T
yax=S.grid.y_T
sout=np.squeeze(S.t_an[0,0,:])

fig=plt.figure(1,figsize=(5,4))
cf=plt.contourf(xax,yax,sout,np.arange(-2,25,1),extend='both')
plt.colorbar(cf)
plt.grid()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Jan SST from NODC WOA09')
plt.savefig('contourf.png')
plt.show()
