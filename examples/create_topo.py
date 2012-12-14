#============================================================
#
# Use supergrids (refine = 1 ) to define tiles for the northern/southern
# and Central mercator grids. find overlapping region of fine-grid
# bathymetry to model grid cells.  Use infomation from overlapping
# region to calculate a (un-weighted) bathymetry, including min/max and std.
#
#============================================================



from midas import *
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import interp 
import argparse

def csv(value):
   return map(int, value.split(","))

#### Begin User Input

parser = argparse.ArgumentParser()
parser.add_argument('--layout',type=csv,help='layout(ny,nx)',default=[1,1])
parser.add_argument('--tile',type=str,help='ncap|scap|mercator')
parser.add_argument('--nx',type=int,help='x domain number 0:layout[1]-1',default=0)
parser.add_argument('--ny',type=int,help='y domain number 0:layout[0]-1',default=0)
parser.add_argument('--plot',type=int,help='[0] plot output',default=0)
parser.add_argument('--do_topo',type=int,help='[0] calculate topography',default=0)

args=parser.parse_args()

tile = args.tile

use_GEBCO = True
use_IBCAO = True
use_BEDMAP = True

RECALC_overlay = False

refine=2
if args.do_topo == 1:
   RECALC_overlay = True
   refine=1
   
do_ncap=False
do_mercator=False
do_scap=False

MAKE_figs = False
if args.plot == 1:
   MAKE_figs = True
   RECALC_overlay = False
   args.layout=[1,1]
   args.nx=0
   args.ny=0
   
   
if tile == 'ncap':
   do_ncap = True
if tile == 'mercator':
   do_mercator = True
if tile == 'scap':   
   do_scap = True

lat0=-75.0
lon0=-280.
lenlat=130.
lenlon=360.
nx=2880*refine
ny=1610*refine
ny_scap=210*refine
ny_ncap=275*refine
nfig=0
min_depth = -10.0



#### End User Input

#### Begin Mercator Grid

print 'constructing a mercator supergrid with (ny,nx) = ',ny,nx
print 'nominal starting lat and starting longitude =',lat0, lon0
print 'and nominal width in latitude = ',lenlat

mercator=supergrid(nx,ny,'mercator','degrees',lat0,lenlat,lon0,lenlon,cyclic_x=True)

if refine == 2:
   mercator.grid_metrics()
   mercator.write_nc('mercator_supergrid.nc')

print "mercator max/min latitude=", mercator.y.max(),mercator.y.min()
print "mercator starting longitude=",mercator.x[0,0]
print "mercator ending longitude=",mercator.x[0,-1]

#### Begin Tripolar Cap

lat0_tp=mercator.y.max()
dlat=90.0-lat0_tp

tripolar_n=supergrid(nx,ny_ncap,'spherical','degrees',lat0_tp,dlat,lon0,360.,tripolar_n=True)

if refine == 2:
   tripolar_n.grid_metrics()
   tripolar_n.write_nc('ncap_supergrid.nc')

print "generated a tripolar supergrid of size (ny,nx)= ",ny_ncap,nx
print "tripolar grid starting longitude = ",tripolar_n.x[0,0]
print "tripolar grid starting latitude = ",lat0_tp



if do_ncap:

   xs=np.float(args.nx)/args.layout[1]*nx
   xe=np.float(args.nx+1)/args.layout[1]*nx + 1
   ys=np.float(args.ny)/args.layout[0]*ny_ncap 
   ye=np.float(args.ny+1)/args.layout[0]*ny_ncap  + 1
   xs=np.int(xs); xe=np.int(xe)
   ys=np.int(ys); ye=np.int(ye)
   
   print 'Processing a limited domain of ncap (ydom,xdom)= ',args.ny,args.nx


   if use_IBCAO:

######## Interpolate bathymetry from IBCAO to northern cap 
######## on a np stereo projection
      
      xlen=2904000.0*2.0
      x=np.linspace(0.0,xlen,11617)
      X,Y=np.meshgrid(x,x)
      grid_ibcao = generic_grid(lon=X,lat=Y,is_latlon=False,is_cartesian=True)

      m = Basemap(projection='stere',width=xlen,height=xlen,lon_0=0.0,lat_0=90.0,resolution='l')

      if  RECALC_overlay is True:
         IBCAO=state('IBCAO_V3_500m_RR.grd',grid=grid_ibcao,fields=['z'])
         IBCAO.rename_field('z','topo')

   xx=tripolar_n.x.copy()
   yy=tripolar_n.y.copy()

   xx[xx>180.]=xx[xx>180.]-360.
   xx[xx<-180.]=xx[xx<-180.]+360.
   x2,y2 = m(xx,yy,inverse=False)

   x2=x2[ys:ye,xs:xe]
   y2=y2[ys:ye,xs:xe]
   
   cart_grid_ncap = supergrid(nx,ny_ncap,'cartesian','none',xdat=x2,ydat=y2)

   fnam = 'ncap_topog_'+str(args.ny)+str(args.nx)+'.nc'
   
   
   if RECALC_overlay:
      R=IBCAO.grid_overlay_struct('topo',target=cart_grid_ncap)
      R.write_nc(fnam,['mean','max','min','std','count'])   
   else:
      R=state('ncap_topog.nc',grid=cart_grid_ncap,fields=['mean','std','count'])
   

   if MAKE_figs:

      clevs_std = [5,10,20,50,100,200,250,300,350]
      clevs_dpt = [10,100,200,300,400,500,1000,1500,2000,2500,3000,3500]      
      nfig=nfig+1
      plt.figure(nfig)
      zout=sq(R.mean)
      std=sq(R.std)
      zout=np.ma.masked_where(zout>min_depth,zout)
      std=np.ma.masked_where(zout.mask,std)
      zout=-zout

      x_plt=0.5*(x2+np.roll(x2,shift=-1,axis=1))
      x_plt=0.5*(x_plt+np.roll(x_plt,shift=-1,axis=0))
      x_plt=x_plt[0:-1,0:-1]

      y_plt=0.5*(y2+np.roll(y2,shift=-1,axis=1))
      y_plt=0.5*(y_plt+np.roll(y_plt,shift=-1,axis=0))
      y_plt=y_plt[0:-1,0:-1]   

      cf=plt.contourf(x_plt,y_plt,zout,clevs_dpt,cmap=cm.GMT_relief,extend='max')
      plt.colorbar(cf)

      for i in np.linspace(0,nx/4,5):
         m.plot(x2[:,i],y2[:,i],color='k')
      for i in np.linspace(nx/4,nx/2,5):
         m.plot(x2[:,i],y2[:,i],color='r')
      for i in np.linspace(nx/2,3*nx/4,5):
         m.plot(x2[:,i],y2[:,i],color='g')
      for i in np.linspace(3*nx/4,nx,5):
         m.plot(x2[:,i],y2[:,i],color='b')
      for i in np.linspace(0,ny_ncap,9):
         m.plot(x2[i,:],y2[i,:],color='k')           

      plt.title('Rotated Grid Bathymetry (m)',fontsize=10)
      plt.savefig('arctic_depth.png')
   
      nfig=nfig+1
      plt.figure(nfig)
      cf=plt.contourf(x_plt,y_plt,std,clevs_std,cmap=cm.GMT_relief,extend='max')   
      plt.colorbar(cf)
      
      for i in np.linspace(0,nx/4,5):
         m.plot(x2[:,i],y2[:,i],color='k')
      for i in np.linspace(nx/4,nx/2,5):
         m.plot(x2[:,i],y2[:,i],color='r')
      for i in np.linspace(nx/2,3*nx/4,5):
         m.plot(x2[:,i],y2[:,i],color='g')
      for i in np.linspace(3*nx/4,nx,5):
         m.plot(x2[:,i],y2[:,i],color='b')
      for i in np.linspace(0,ny_ncap,9):
         m.plot(x2[i,:],y2[i,:],color='k')

      
      plt.title('Rotated grid roughness height (m)',fontsize=10)
      plt.savefig('arctic_roughness.png')


#### Begin Mercator
   
if do_mercator:

   xs=np.float(args.nx)/args.layout[1]*nx
   xe=np.float(args.nx+1)/args.layout[1]*nx + 1
   ys=np.float(args.ny)/args.layout[0]*ny
   ye=np.float(args.ny+1)/args.layout[0]*ny + 1

   xs=np.int(xs); xe=np.int(xe)
   ys=np.int(ys); ye=np.int(ye)
   
   print 'Processing a limited domain of mercator (ydom,xdom)= ',args.ny,args.nx
   print 'nx=',nx
   print 'xs,xe= ',xs,xe
   print 'ny=',ny
   print 'ys,ye=',ys,ye

   mercator.x=mercator.x[ys:ye,xs:xe]
   mercator.y=mercator.y[ys:ye,xs:xe]
   mercator.grid_x=mercator.grid_x[xs:xe]
   mercator.grid_y=mercator.grid_y[ys:ye]      
   
   if use_GEBCO:
      ingrid=generic_grid('GEBCO_08_v1.nc',var='depth',simple_grid=True,cyclic=True)
      merc_reg=ingrid.geo_region(y=(mercator.y.min()-1.0,mercator.y.max()+1.0))

      if  RECALC_overlay is True:
         TOPO=state('GEBCO_08_v1.nc',grid=ingrid,geo_region=merc_reg,fields=['depth'])
         TOPO.rename_field('depth','topo')
         TOPO.var_dict['topo']['Ztype']='Fixed'

    
   fnam = 'mercator_topog_'+str(args.ny)+str(args.nx)+'.nc'
   
   if RECALC_overlay:
      R=TOPO.grid_overlay_struct('topo',target=mercator)
      R.write_nc(fnam,['mean','max','min','std','count'])
   else:
      R=state('mercator_topog.nc',grid=mercator,fields=['mean','std','count'])


   if MAKE_figs:

      clevs_std = [5,10,20,50,100,200,250,300,350,400,450,500]
      clevs_dpt = [10,100,200,300,400,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000]
      
      nfig=nfig+1
      plt.figure(nfig)

      zout=sq(R.mean)
      std=sq(R.std)
      zout=np.ma.masked_where(zout>min_depth,zout)
      std=np.ma.masked_where(zout.mask,std)
      zout=-zout
      x_plt=R.var_dict['mean']['xax_data']
      y_plt=R.var_dict['mean']['yax_data']


      cf=plt.contourf(x_plt,y_plt,zout,clevs_dpt,cmap=cm.GMT_relief,extend='max')   
      plt.colorbar(cf)
      plt.title('Mercator Grid Bathymetry (m)',fontsize=8)
      plt.savefig('mercator_depth.png')

      nfig=nfig+1
      plt.figure(nfig)
      cf=plt.contourf(x_plt,y_plt,std,clevs_std,cmap=cm.GMT_relief,extend='max')   
      plt.colorbar(cf)
      plt.title('Mercator grid roughness height (m)',fontsize=10)
      plt.savefig('mercator_roughness.png')
      


#### Begin Antarctic Cap

lenlat=90.0+mercator.y.min()
antarctic_cap=supergrid(nx,ny_scap,'spherical','degrees',-90.,lenlat,lon0,lenlon,displace_pole=True,r0_pole=0.45,lon0_pole=100.,doughnut=0.45)

if refine == 2:
   antarctic_cap.grid_metrics()
   antarctic_cap.write_nc('scap_supergrid.nc')

if do_scap:   

   xs=np.float(args.nx)/args.layout[1]*nx
   xe=np.float(args.nx+1)/args.layout[1]*nx + 1
   ys=np.float(args.ny)/args.layout[0]*ny_scap
   ye=np.float(args.ny+1)/args.layout[0]*ny_scap + 1

   xs=np.int(xs); xe=np.int(xe)
   ys=np.int(ys); ye=np.int(ye)
   
   print 'Processing a limited domain of scap (ydom,xdom)= ',args.ny,args.nx
   print 'nx=',nx
   print 'xs,xe= ',xs,xe
   print 'ny=',ny
   print 'ys,ye=',ys,ye



   if use_BEDMAP:
      f=nc.Dataset('Antarctica_5km_withshelves_v0.75.nc')
      lon1=sq(f.variables['lon'][:])
      lat1=sq(f.variables['lat'][:])
      x1=sq(f.variables['x1'][:])
      y1=sq(f.variables['y1'][:])
      nx1=len(x1)
      ny1=len(y1)
      wd=x1[-1]-x1[0]+1
      ht=y1[-1]-y1[0]+1
      lat0=lat1[ny1/2,nx1/2]
      m = Basemap(projection='stere',width=wd,height=ht,lon_0=0.0,lat_0=lat0,resolution='l')
      X,Y = m(lon1,lat1,inverse=False)
      grid_bedmap = generic_grid(lon=X,lat=Y,is_latlon=False,is_cartesian=True)

      if  RECALC_overlay is True:
         TOPO=state('Antarctica_5km_withshelves_v0.75.nc',grid=grid_bedmap,fields=['topg'])
         TOPO.rename_field('topg','topo')

   xx=antarctic_cap.x.copy()
   yy=antarctic_cap.y.copy()

   xx[xx>180.]=xx[xx>180.]-360.
   xx[xx<-180.]=xx[xx<-180.]+360.
   x2,y2 = m(xx,yy,inverse=False)

   x2=x2[ys:ye,xs:xe]
   y2=y2[ys:ye,xs:xe]
   
   cart_grid_scap = supergrid(nx,ny_scap,'cartesian','none',xdat=x2,ydat=y2)

   fnam = 'scap_topog_'+str(args.ny)+str(args.nx)+'.nc'
   
   if RECALC_overlay:
      R=TOPO.grid_overlay_struct('topo',target=cart_grid_scap)
      R.write_nc(fnam,['mean','max','min','std','count'])   
   else:
      R=state('scap_topog.nc',grid=cart_grid_scap,fields=['mean','std','count'])

   if MAKE_figs:

      clevs_std = [5,10,20,50,100,200,250,300,350]
      clevs_dpt = [10,100,200,300,400,500,1000,1500,2000,2500,3000,3500]
      
      nfig=nfig+1
      plt.figure(nfig)

      zout=sq(R.mean)
      std=sq(R.std)
      zout=np.ma.masked_where(zout>min_depth,zout)
      std=np.ma.masked_where(zout.mask,std)
      zout=-zout
      
      x_plt=0.5*(x2+np.roll(x2,shift=-1,axis=1))
      x_plt=0.5*(x_plt+np.roll(x_plt,shift=-1,axis=0))
      x_plt=x_plt[0:-1,0:-1]

      y_plt=0.5*(y2+np.roll(y2,shift=-1,axis=1))
      y_plt=0.5*(y_plt+np.roll(y_plt,shift=-1,axis=0))
      y_plt=y_plt[0:-1,0:-1]   
   
      cf=plt.contourf(x_plt,y_plt,zout,clevs_dpt,cmap=cm.GMT_relief,extend='max')   
      plt.colorbar(cf)

      for i in np.linspace(0,nx,20):
         plt.plot(x2[:,i],y2[:,i],color='k')

      for i in np.linspace(0,y2.shape[0]-1,5):
         plt.plot(x2[i,:],y2[i,:],color='r')
      
      plt.title('Antarctic Grid Bathymetry (m)',fontsize=8)
      plt.savefig('antarctic_depth.png')

      nfig=nfig+1
      plt.figure(nfig)
      cf=plt.contourf(x_plt,y_plt,std,clevs_std,cmap=cm.GMT_relief,extend='max')   
      plt.colorbar(cf)
      
      for i in np.linspace(0,nx,20):
         plt.plot(x2[:,i],y2[:,i],color='k')
      for i in np.linspace(0,y2.shape[0]-1,5):
         plt.plot(x2[i,:],y2[i,:],color='r')

   
      plt.title('Antarctic grid roughness height (m)',fontsize=10)
      plt.savefig('antarctic_roughness.png')



