from midas.rectgrid import *
import matplotlib.pyplot as plt
import numpy as np

print("""\n
         Produce a Area-weighted histogram
         from AMSRE SST data in the Indian Ocean\n""")

fpath='https://opendap.jpl.nasa.gov:443/opendap/OceanTemperature/amsre/L3/sst_1deg_1mo/tos_AMSRE_L3_v7_200206-201012.nc'
grid=quadmesh(fpath,var='tos')
iop=grid.geo_region(x=(30,160),y=(-10,25))
S=state(fpath,grid=grid,geo_region=iop,fields=['tos'],default_calendar='noleap')

xax=S.grid.lonh
yax=S.grid.lath

fig=plt.figure(1,figsize=(5,4))

cf=plt.contourf(xax,yax,sq(S.tos[0,0,:]),np.arange(299,305,.25),extend='both',cmap=plt.cm.jet)

plt.colorbar(cf)
plt.xlim(30,160)
plt.ylim(-10,25)
plt.grid()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('AMSRE SST')


fig=plt.figure(2)

b = sq(S.tos)
area=S.grid.Ah
area=area[np.newaxis,np.newaxis,:]
vol = np.tile(S.grid.Ah,(S.tos.shape[0],1,1,1))

wt = vol.flatten()
b=b.flatten()
c=b*wt
c=np.ma.masked_where(b.mask,c)
wt=np.ma.masked_where(b.mask,wt)
bm= np.sum(c)/np.sum(wt)


n,bins,patches = plt.hist(b,bins=np.arange(290,303,.25),weights=wt,normed=True,cumulative=False,color='blue')
plt.xlim(299,305)
#plt.plot([bm,bm],[0,5],'r-')
plt.grid()
plt.title('Histogram of Area-weighted SST in the Indian Ocean (b=total;r=sfc)',fontsize=8)
plt.xlabel('SST (degK)',fontsize=8)
plt.ylabel('Frequency',fontsize=8)


plt.show()

