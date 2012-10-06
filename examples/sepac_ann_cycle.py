from midas import *

grid=gold_grid('/net3/mjh/models/CM2G/ocean_geometry.nc')
sepac=grid.geo_region(x=(-240,-69),y=(-40,20))


SCOW=state('/archive/mjh/riga_201012/ClimSCOW-G/gfdl.default-prod/pp/ocean/ts/monthly/25yr/ocean.217501-219912.temp.nc',fields=['temp'],grid=grid,geo_region=sepac,path_interfaces='/archive/mjh/riga_201012/ClimSCOW-G/gfdl.default-prod/pp/ocean/ts/monthly/25yr/ocean.217501-219912.e.nc',interfaces='e')

CORE=state('/archive/mjh/riga_201012/GOLD_SIS_63L_prod/gfdl.default-prod/pp/ocean/ts/monthly/20yr/ocean.206001-207912.temp.nc',fields=['temp'],grid=grid,geo_region=sepac,path_interfaces='/archive/mjh/riga_201012/GOLD_SIS_63L_prod/gfdl.default-prod/pp/ocean/ts/monthly/20yr/ocean.206001-207912.e.nc',interfaces='e')

#CM2G=state(MFpath='/archive/mjh/CM2G63L/CM2G_1860_Control/gfdl.ncrc2-default-prod-openmp/pp/ocean/ts/monthly/5yr/ocean.10??01-10??12.temp.nc',fields=['temp'],grid=grid,geo_region=sepac,MFpath_interfaces='/archive/mjh/CM2G63L/CM2G_1860_Control/gfdl.ncrc2-default-prod-openmp/pp/ocean/ts/monthly/5yr/ocean.10??01-10??12.e.nc',interfaces='e')

tdict=SCOW.var_dict['temp']
tdict['Ztype']='Fixed'
vdict=SCOW.var_dict['temp'].copy()
vdict['z']=[0];vdict['dz']=None;vdict['z_interfaces']=None;vdict['Z']=None
SCOW.create_field('-np.reshape(self.e[:,4,:],(self.e.shape[0],1,self.e.shape[2],self.e.shape[3]))','h_ML',var_dict=vdict)
SCOW.create_field('np.reshape(self.temp[:,0,:],(self.temp.shape[0],1,self.temp.shape[2],self.temp.shape[3]))','sst',var_dict=vdict)
SCOW.monthly_avg('sst',vol_weight=False)
SCOW.monthly_avg('h_ML',vol_weight=False)


tdict=CORE.var_dict['temp']
tdict['Ztype']='Fixed'
vdict=CORE.var_dict['temp'].copy()
vdict['z']=[0];vdict['dz']=None;vdict['z_interfaces']=None;vdict['Z']=None
CORE.create_field('-np.reshape(self.e[:,4,:],(self.e.shape[0],1,self.e.shape[2],self.e.shape[3]))','h_ML',var_dict=vdict)
CORE.create_field('np.reshape(self.temp[:,0,:],(self.temp.shape[0],1,self.temp.shape[2],self.temp.shape[3]))','sst',var_dict=vdict)
CORE.monthly_avg('sst',vol_weight=False)
CORE.monthly_avg('h_ML',vol_weight=False)


vdict=CM2G.var_dict['temp'].copy()
vdict['z']=[0];vdict['dz']=None;vdict['z_interfaces']=None;vdict['Z']=None
CM2G.create_field('-np.reshape(self.e[:,4,:],(self.e.shape[0],1,self.e.shape[2],self.e.shape[3]))','h_ML',var_dict=vdict)
CM2G.create_field('np.reshape(self.temp[:,0,:],(self.temp.shape[0],1,self.temp.shape[2],self.temp.shape[3]))','sst',var_dict=vdict)
CM2G.monthly_avg('sst',vol_weight=False)
CM2G.monthly_avg('h_ML',vol_weight=False)

xax=CORE.grid.lonh
yax=CORE.grid.lath

fig=plt.figure(1,figsize=(8.5,11))

ax1=fig.add_subplot(311)
cf=ax1.contourf(xax,yax,sq(np.min(CORE.h_ML_monthly,axis=0)),np.concatenate((np.arange(2,10,2),np.arange(10,50,10))),extend='max')
plt.colorbar(cf,ax=ax1)
plt.title('CORE-G Monthly Minumum MLD (m)')
plt.grid()

ax2=fig.add_subplot(312)
cf=ax2.contourf(xax,yax,sq(np.min(SCOW.h_ML_monthly,axis=0)),np.concatenate((np.arange(2,10,2),np.arange(10,50,10))),extend='max')
plt.colorbar(cf,ax=ax2)
plt.title('SCOW-G Monthly Minumum MLD (m)')
plt.grid()

ax3=fig.add_subplot(313)
cf=ax3.contourf(xax,yax,sq(np.min(CM2G.h_ML_monthly,axis=0)),np.concatenate((np.arange(2,10,2),np.arange(10,50,10))),extend='max')
plt.colorbar(cf,ax=ax3)
plt.title('CM2G Monthly Minumum MLD (m)')
plt.grid()

plt.show()
