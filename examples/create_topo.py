from midas import *
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import interp 

PLOT_orig_topo = False
use_GEBCO=False
use_OCCAM=True

lat0=-75.0
lon0=-280.
lenlat=130.
lenlon=360.
nx=2880
ny=1600
ny_scap=192
nfig=0

print 'constructing a mercator grid with (ny,nx) = ',ny,nx
print 'nominal starting lat and starting longitude =',lat0, lon0
print 'and nominal width in latitude = ',lenlat

mercator=supergrid(nx,ny,'mercator','degrees',lat0,lenlat,lon0,lenlon)


print "mercator max/min latitude=", mercator.y.max(),mercator.y.min()
print "mercator starting longitude=",mercator.x[0,0]
print "mercator ending longitude=",mercator.x[0,-1]

lenlat=90.0+mercator.y.min()

antarctic_cap=supergrid(nx,ny_scap,'spherical','degrees',-90.,lenlat,lon0,lenlon,)


r,phi = antarctic_cap.displaced_sp(0.45,100.,0.55)

antarctic_cap.x=phi.copy()
antarctic_cap.y=r.copy()
antarctic_cap.grid_x=antarctic_cap.x[-1,:]
antarctic_cap.grid_y=antarctic_cap.y[:,nx/4]

#antarctic_grid = generic_grid(supergrid=antarctic_cap,refine=1)

f=nc.Dataset('data/Antarctica_5km_withshelves_v0.75.nc')
lon1=sq(f.variables['lon'][:])
lat1=sq(f.variables['lat'][:])
x1=sq(f.variables['x1'][:])
y1=sq(f.variables['y1'][:])
nx=len(x1)
ny=len(y1)
wd=x1[-1]-x1[0]+1
ht=y1[-1]-y1[0]+1
lat0=lat1[ny/2,nx/2]
topg=sq(f.variables['topg'][:])
print 'lat0=',lat0
m = Basemap(projection='stere',width=wd,height=ht,lon_0=0.0,lat_0=lat0,resolution='l')
x2,y2 = m(lon1,lat1,inverse=False)

lon3=antarctic_cap.x.copy()
lon3[lon3>180.]=lon3[lon3>180.]-360.
lon3[lon3<-180.]=lon3[lon3<-180.]+360.
lat3=antarctic_cap.y.copy()
x3,y3 = m(lon3,lat3,inverse=False)

if PLOT_orig_topo:
    nfig=nfig+1
    plt.figure(nfig)
    m.pcolormesh(x2,y2,topg,vmin=-6500.,vmax=3500.)
    
    for i in np.linspace(0,x3.shape[0]-1,8).astype(int):
        m.plot(x3[i,:],y3[i,:],color='k')
            
    for i in np.linspace(0,x3.shape[1]-1,8).astype(int):
        m.plot(x3[:,i],y3[:,i],color='r')

    for i in np.linspace(0,x2.shape[0]-1,8).astype(int):
        m.plot(x2[i,:],y2[i,:],color='g')

    for i in np.linspace(0,x2.shape[1]-1,8).astype(int):
        m.plot(x2[:,i],y2[:,i],color='g')


        plt.savefig('orig_topg_antarctica.png')

data_out=interp(topg,sq(x2[0,:]),sq(y2[:,0]),x3,y3,order=3)            
nfig=nfig+1
plt.figure(nfig)
m.pcolormesh(x3,y3,data_out,vmin=-6500,vmax=3500)
    
m.contour(x3,y3,data_out,[0.0,0.0],colors='b',linewidth=3.0)

m.contour(x3,y3,-data_out,[500.0,500.0],colors='c',linewidth=2.0)

lon4=mercator.x.copy()
lon4[lon4>180.]=lon4[lon4>180.]-360.
lon4[lon4<-180.]=lon4[lon4<-180.]+360.
lat4=mercator.y.copy()

print 'Starting longitude for mercator grid = ',mercator.x[0,0]
print 'Starting longitude at last row of sp grid = ',antarctic_cap.x[-1,0]

x4,y4 = m(lon4,lat4,inverse=False)

data_out2=interp(topg,sq(x2[0,:]),sq(y2[:,0]),x4,y4,order=3)


m.pcolormesh(x4,y4,data_out2,vmin=-6500,vmax=3500)

m.contour(x4,y4,data_out2,[0.0,0.0],colors='b',linewidth=3.0)
m.contour(x4,y4,-data_out2,[500.0,500.0],colors='c',linewidth=2.0)

for i in np.linspace(0,x3.shape[0]-1,8).astype(int):
    m.plot(x3[i,:],y3[i,:],color='k')

m.plot(x3[-1,:],y3[-1,:],color='g',linewidth=4.0)

for i in np.linspace(0,x3.shape[1]-1,8).astype(int):
    m.plot(x3[:,i],y3[:,i],color='r')

for i in np.linspace(0,x4.shape[0]-1,8).astype(int):
    m.plot(x4[i,:],y4[i,:],color='k')

for i in np.linspace(0,x4.shape[1]-1,8).astype(int):
    m.plot(x4[:,i],y4[:,i],color='r')


xp,yp=m(0.,-90,inverse=False)
plt.plot(xp,yp,'+',color='k')
plt.title('Bedrock Bathymetry for Antarctic ice-shelf model')
plt.savefig('topg_interpolated_antarctic.png')


if use_GEBCO:
    ingrid=generic_grid('/tmp/GEBCO_08_v1.nc',var='depth',simple_grid=True)
    merc_reg=ingrid.geo_region(y=(mercator.y.min()-0.1,mercator.y.max()+0.1))
    TOPO=state('/tmp/GEBCO_08_v1.nc',grid=ingrid,geo_region=merc_reg,fields=['depth'])
    TOPO.rename_field('depth','topo')
    TOPO.var_dict['topo']['Ztype']='Fixed'
elif use_OCCAM:
    ingrid=generic_grid('/net3/mjh/data/occam_topo/OCCAM_topo.nc',var='TOPO',simple_grid=True)
    merc_reg=ingrid.geo_region(y=(mercator.y.min()-0.1,mercator.y.max()+0.1))
    TOPO=state('/net3/mjh/data/occam_topo/OCCAM_topo.nc',grid=ingrid,geo_region=merc_reg,fields=['TOPO'])
    TOPO.rename_field('TOPO','topo')    
    TOPO.var_dict['topo']['Ztype']='Fixed'    
    
R=TOPO.grid_overlay_struct('topo',target=mercator,src_modulo=True)

R.write_nc('mercator_topog_struct.nc',['mean','max','min','std','count'])

nfig=nfig+1
plt.figure(nfig)

cf=plt.pcolormesh(mercator.x,mercator.y,sq(R.std))
plt.colorbar(cf)
plt.title('Mercator Grid Fine Grid RMS (m)')
plt.show()




