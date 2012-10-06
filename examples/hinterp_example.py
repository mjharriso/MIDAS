from midas import *
grid=gold_grid('/net3/mjh/models/GOLD_SIS_025/ocean_geometry.nc')

grid_obs=generic_grid(path='/net3/mjh/models/ESM2.5G/input_data/n00an1.nc',var='n00an1')

S=state(path='/net3/mjh/models/ESM2.5G/input_data/n00an1.nc',grid=grid_obs,geo_region=None,fields=['n00an1'])
T=S.horiz_interp('n00an1',target=grid,src_modulo=True,method='bilinear')
xax=T.grid.x_T
yax=T.grid.y_T
sout=np.squeeze(T.n00an1[0,0,:])
plt.pcolormesh(xax,yax,sout)
plt.contour(xax,yax,T.grid.wet,[0.5,0.5],colors='b')
plt.contour(xax,yax,T.grid.D,np.arange(500,6500,500),colors='k',alpha=0.3)
plt.show()

