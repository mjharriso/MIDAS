from midas import *


xlon=np.arange(0.5,360.5,1.)
ylat=np.arange(-39.5,40.5,1.)
lon,lat = np.meshgrid(xlon,ylat)
grid_out = generic_rectgrid(lon=lon,lat=lat,is_latlon=True)
grid=generic_rectgrid('http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/8day',var='AAssta')
S=state(grid=grid,fields=['AAssta'],path='http://oceanwatch.pfeg.noaa.gov/thredds/dodsC/satellite/AA/ssta/8day',date_bounds=[datetime(2011,1,1),datetime(2011,1,2)],default_calendar='gregorian')
S.rename_field('AAssta','sst')
R=S.grid_overlay(field='sst',target=grid_out)

plt.figure(1)
xax=S.grid.lonh
yax=S.grid.lath
cf=plt.pcolormesh(xax,yax,sq(S.sst),vmin=4,vmax=30)
plt.colorbar(cf)
plt.xlim(0,360)
plt.ylim(-40,40)
plt.title('Original Data on 1/4 degree grid')

plt.figure(2)
xax=grid_out.lonh
yax=grid_out.lath
cf=plt.pcolormesh(xax,yax,R.mean,vmin=4,vmax=30)
plt.colorbar(cf)
plt.xlim(0,360)
plt.ylim(-40,40)
plt.title('Mean on 1 degree grid')

plt.figure(3)
cf=plt.pcolormesh(xax,yax,np.ma.masked_where(R.count==0,R.count))
plt.colorbar(cf)
plt.xlim(0,360)
plt.ylim(-40,40)
plt.title('Number of input points')

plt.figure(4)
cf=plt.pcolormesh(xax,yax,R.std)
plt.colorbar(cf)
plt.xlim(0,360)
plt.ylim(-40,40)
plt.title('Standard deviation')

plt.show()
