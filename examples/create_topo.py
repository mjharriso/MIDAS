from midas import *
from mpl_toolkits.basemap import Basemap

use_GEBCO = False
use_CISM = True


lat0=-75.0
lon0=-280.
lenlat=130.
lenlon=360.
#nx=5760
#ny=2240
#ny_scap=480
nx=720
ny=400
ny_scap=48

print 'constructing a mercator grid with (ny,nx) = ',ny,nx
print 'nominal starting lat and starting longitude =',lat0, lon0
print 'and nominal width in latitude = ',lenlat

mercator=supergrid(nx,ny,'mercator','degrees',lat0,lenlat,lon0,lenlon)

print "mercator max/min latitude=", mercator.y.max(),mercator.y.min()
print "mercator starting longitude=",mercator.x[0,0]
print "mercator ending longitude=",mercator.x[0,-1]

lenlat=90.0+mercator.y.min()

antarctic_cap=supergrid(nx,ny_scap,'spherical','degrees',-90.,lenlat,lon0,lenlon,)


r,phi = antarctic_cap.displaced_sp(0.25,130.,0.2)

antarctic_cap.x=phi.copy()
antarctic_cap.y=r.copy()
antarctic_cap.grid_x=antarctic_cap.x[-1,:]
antarctic_cap.grid_y=antarctic_cap.y[:,nx/4]

antarctic_grid = generic_grid(supergrid=antarctic_cap,refine=1)

if use_GEBCO:
    gebco_grid=generic_grid('/tmp/GEBCO_08_v1.nc',var='depth',simple_grid=True)
    gebco_grid.lonq[-1]=360.0-gebco_grid.lonq[-1]
    gebco_grid.latq[0]=-90.0
    s_pole=gebco_grid.geo_region(y=(-90.,antarctic_cap.y.max()+0.1))
    GEBCO=state('/tmp/GEBCO_08_v1.nc',grid=gebco_grid,geo_region=s_pole,fields=['depth'])

    GEBCO.var_dict['depth']['Ztype']='Fixed'

    R=GEBCO.horiz_interp('depth',target=antarctic_grid,src_modulo=True,add_NP=False,method='bilinear')

if use_CISM:
    f=nc.Dataset('/net3/mjh/data/Antarctica_bedmap/Antarctica_5km_withshelves_v0.75.nc')
    lon_in=sq(f.variables['lon'][:])
    lat_in=sq(f.variables['lat'][:])
    grid_in=generic_grid(lon=lon_in,lat=lat_in)
    grid_in.x_T[grid_in.x_T+180.0 < 0.12]=-180.
    BEDMAP=state('/net3/mjh/data/Antarctica_bedmap/Antarctica_5km_withshelves_v0.75.nc',grid=grid_in,fields=['topg'])
    BEDMAP.rename_field('topg','depth')
    BEDMAP.var_dict['depth']['Ztype']='Fixed'

    R=BEDMAP.horiz_interp('depth',target=antarctic_grid,src_modulo=True,add_NP=False,method='bilinear')

    
lon=antarctic_grid.x_T
lat=antarctic_grid.y_T
depth=sq(R.depth)

#lon=grid_in.x_T_bounds
#lat=grid_in.y_T_bounds
#depth=sq(BEDMAP.depth)

print depth.max()
print depth.min()
print depth.mean()

m=Basemap(projection='spstere',lon_0=0.,boundinglat=antarctic_cap.y.max())
x,y=m(lon,lat)
cf=m.contourf(x,y,sq(depth),np.arange(-6500,3750,250),extend='both',cmap=plt.cm.jet)
#cf=m.pcolormesh(x,y,sq(depth),vmin=-6500,vmax=3500,cmap=plt.cm.jet)
m.contour(x,y,sq(depth),np.arange(-6000,6000,1000),colors='k',alpha=0.6)
m.colorbar(cf)
m.contour(x,y,sq(depth),[0.0,0.0],colors='k',linewidth=3.0)

m.plot(sq(x[:,0]),sq(y[:,0]),color='g',linewidth=3.0)
for j in np.arange(0,lat.shape[0],4):
    m.plot(sq(x[j,:]),sq(y[j,:]),color='k')

for i in np.arange(0,lat.shape[1],4):
    m.plot(sq(x[:,i]),sq(y[:,i]),color='k')    
    
#m.drawcoastlines()

if use_GEBCO:
    plt.title('Bathymetry from GEBCO (m) ')
    plt.savefig('antarctica_gebco.png')    
if use_CISM:
    plt.title('Bathymetry from CISM (m) ')        
    plt.savefig('antarctica_cism.png')    


plt.show()



