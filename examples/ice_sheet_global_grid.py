from midas import *
from midas_grid_gen import *


PI_180 = np.pi / 180.


lat0=-78.0
lon0=-280.
lenlat=140.
lenlon=360.
nx=360
ny=180

print 'constructing a mercator grid with (ny,nx) = ',ny,nx
print 'nominal starting lat and starting longitude =',lat0, lon0
print 'and nominal width in latitude = ',lenlat

mercator=supergrid(nx,ny,'mercator','degrees',lat0,lenlat,lon0,lenlon)

print "mercator max/min latitude=", mercator.y.max(),mercator.y.min()
print "mercator starting longitude=",mercator.x[0,0]
print "mercator ending longitude=",mercator.x[0,-1]

lenlat=90.0+mercator.y.min()

antarctic_cap=supergrid(360,20,'spherical','degrees',-90.,lenlat,lon0,lenlon,)


r,phi = antarctic_cap.displaced_sp(0.6,200.,0.5)

antarctic_cap.x=phi.copy()
antarctic_cap.y=r.copy()

ymax=antarctic_cap.y.max()

phi=phi*PI_180
r=r+90.0
r=r/r.max()

plt.clf()

fig=plt.figure(1)

ax1=fig.add_subplot(111)

for  j in np.arange(0,phi.shape[0]):
    plt.polar(phi[j,:],r[j,:],'r')

for  i in np.arange(0,phi.shape[1]):
    plt.polar(phi[:,i],r[:,i],'k')

plt.polar(phi[:,0],r[:,0],'g',linewidth=3.0)

plt.show()

fig2=plt.figure(2)

ny,nx=mercator.x.shape

print "max/min mercator longitude = ",mercator.x.min(),mercator.x.max()
print "max/min mercator latitude = ",mercator.y.min(),mercator.y.max()

ax2=fig2.add_subplot(111)

ax2.plot(mercator.grid_y)

fig3=plt.figure(3)

ax3=fig3.add_subplot(111)

lat0=mercator.grid_y[-1]

dlat = 90.0-lat0

tripolar_n=supergrid(360,40,'spherical','degrees',lat0,dlat,lon0,360.,tripolar_n=True)

ny,nx = tripolar_n.x.shape

print "generated a tripolar grid of size (ny,nx)= ",ny,nx
print "tripolar grid starting longitude = ",tripolar_n.x[0,0]
print "tripolar grid starting latitude = ",lat0




theta = tripolar_n.x*PI_180
r=(90.0-tripolar_n.y)/dlat

for i in np.arange(0,nx/4):
   plt.polar(theta[:,i],r[:,i],color='k')
for i in np.arange(nx/4,nx/2):
   plt.polar(theta[:,i],r[:,i],color='r')
for i in np.arange(nx/2,3*nx/4):
   plt.polar(theta[:,i],r[:,i],color='g')
for i in np.arange(3*nx/4,nx):
   plt.polar(theta[:,i],r[:,i],color='b')        


for j in np.arange(0,ny):
   plt.polar(theta[j,:],r[j,:],'k')
          
                                            

plt.show()    

    
