from midas.rectgrid import *
import numpy as np
import matplotlib.pyplot as plt


xlon=np.linspace(280.5,289.5,40)
ylat=np.linspace(35.5,44.5,40)
lon,lat = np.meshgrid(xlon,ylat)
grid_out = quadmesh(lon=lon,lat=lat,is_latlon=True,cyclic=True)
grid=quadmesh('http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/smith_sandwell_topo.nc',var='ROSE',simple_grid=True,cyclic=True)
reg=grid.geo_region(x=(280.0,290.),y=(35.,45.))
S=state(grid=grid,geo_region=reg,fields=['ROSE'],path='http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/smith_sandwell_topo.nc')
R=S.subtile(field='ROSE',target=grid_out)

fig=plt.figure(1)
ax1=fig.add_subplot(221)
ax2=fig.add_subplot(222)
ax3=fig.add_subplot(223)
ax4=fig.add_subplot(224)

xax0=S.grid.lonh
yax0=S.grid.lath
cf=ax1.contourf(xax0,yax0,sq(S.ROSE),np.linspace(-3500.,1000.,30),cmap=plt.cm.terrain)
ax1.contour(xax0,yax0,sq(S.ROSE),[0.0,0.0],colors='k',linewidth=2.0)
ax1.set_xlim(280,290)
ax1.set_ylim(35,45)
#plt.colorbar(cf,ax=ax1)
ax1.set_title('Original 2min Bathymetry')

xax=grid_out.lonh
yax=grid_out.lath
cf=ax2.contourf(xax,yax,sq(R.mean),np.linspace(-3500.,1000.,30),cmap=plt.cm.terrain)
ax2.contour(xax0,yax0,sq(S.ROSE),[0.0,0.0],colors='k',linewidth=2.0)
ax2.contour(xax,yax,sq(R.mean),[0.0,0.0],colors='m',linewidth=2.0)
#plt.colorbar(cf,ax=ax2)
ax2.set_xlim(280,290)
ax2.set_ylim(35,45)
ax2.set_title('Un-weighted cell average Bathymetry')

cf=ax3.contourf(xax,yax,sq(R.max),np.linspace(-3500.,1000.,30),cmap=plt.cm.terrain)
#plt.colorbar(cf,ax=ax3)
ax3.contour(xax0,yax0,sq(S.ROSE),[0.0,0.0],colors='k',linewidth=2.0)
ax3.contour(xax,yax,sq(R.max),[0.0,0.0],colors='m',linewidth=2.0)
ax3.set_xlim(280,290)
ax3.set_ylim(35,45)
ax3.set_title('Cell Max Bathymetry')


cf=ax4.contourf(xax,yax,sq(R.std),np.linspace(2.0,200.0,20),cmap=plt.cm.jet)
#plt.colorbar(cf,ax=ax4)
ax4.contour(xax0,yax0,sq(S.ROSE),[0.0,0.0],colors='k',linewidth=2.0)
ax4.contour(xax,yax,sq(R.mean),[0.0,0.0],colors='m',linewidth=2.0)
ax4.set_xlim(280,290)
ax4.set_ylim(35,45)
ax4.set_title('Cell Sub-grid de-trended RMS')

plt.show()
