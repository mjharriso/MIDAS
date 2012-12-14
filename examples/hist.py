from midas import *
import matplotlib.pylab as pylab

print """\n
         Produce a volume-weighted histogram
         from WOA09 salinity data in the Indian Ocean\n"""

grid=generic_rectgrid('http://data.nodc.noaa.gov/thredds/dodsC/woa09/salinity_annual_1deg.nc',var='s_an')
iop=grid.geo_region(x=(30,160),y=(-30,25))
S=state('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/salinity_annual_1deg.nc',grid=grid,geo_region=iop,fields=['s_an'],default_calendar='noleap')

xax=S.grid.lonh
yax=S.grid.lath

fig=plt.figure(1,figsize=(5,4))

cf=plt.contourf(xax,yax,sq(S.s_an[0,0,:]),np.arange(32,36.5,.1),extend='both',cmap=plt.cm.jet)
plt.contour(xax,yax,sq(S.s_an[0,0,:]),np.arange(32,36.5,.5),colors='k')

plt.colorbar(cf)
plt.xlim(30,160)
plt.ylim(-30,25)
plt.grid()
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('WOA09 Surface Salinity (psu)')
plt.savefig('sss.png')

b = sq(S.s_an)
area=S.grid.Ah
area=area[np.newaxis,:]
vol = np.tile(S.grid.Ah,(33,1,1))*S.var_dict['s_an']['dz']

wt = vol.flatten()
b=b.flatten()
c=b*wt
c=np.ma.masked_where(b.mask,c)
wt=np.ma.masked_where(b.mask,wt)
bm= np.sum(c)/np.sum(wt)


plt.clf()

n,bins,patches = plt.hist(b,bins=np.arange(32,36.5,.1),weights=wt,normed=True,cumulative=False,color='blue')

#n=n/np.sum(n*np.diff(bins))
#bins_mid = 0.5*(bins + np.roll(bins,shift=1))
#bins_mid = bins_mid[0:-1]
#plt.bar(bins_mid,n)

plt.xlim(34,36.5)
plt.plot([bm,bm],[0,5],'r-')
plt.grid()
plt.title('Histogram of Volume-weighted Salinity in the Indian Ocean',fontsize=8)
plt.xlabel('Salinity (g/kg)',fontsize=8)
plt.ylabel('Frequency',fontsize=8)

S=state('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/salinity_annual_1deg.nc',grid=grid,geo_region=iop,fields=['s_an'],default_calendar='noleap',z_indices=np.arange(0,15))

b = sq(S.s_an)
area=S.grid.Ah
area=area[np.newaxis,:]

vol = np.tile(S.grid.Ah,(15,1,1))*S.var_dict['s_an']['dz']
wt = vol.flatten()
b=b.flatten()
c=b*wt
c=np.ma.masked_where(b.mask,c)
wt=np.ma.masked_where(b.mask,wt)
bm= np.sum(c)/np.sum(wt)

n,bins,patches = plt.hist(b,bins=np.arange(32,36.5,.1),weights=wt,normed=True,cumulative=False,color='red',alpha=0.6)

plt.savefig('hist.png')

