from midas import *
grid=gold_grid('/net3/mjh/models/CM2G/ocean_geometry.nc')
grid.add_mask('tmask','/net3/mjh/models/CM2G/regionmask.nc')
S=state(path='/archive/mjh/CM2G63L/CM2G_1860_Control/gfdl.ncrc2-default-prod-openmp/history/unpacked/06000101.ocean_month.nc',grid=grid,geo_region=None,fields=['temp','salt','SW','PmE','LwLatSens','salt_flux'],z_indices=np.arange(0,1))
S.mask_where('temp','grid.mask!=2.0')
S.mask_where('salt','grid.mask!=2.0')
S.mask_where('SW','grid.mask!=2.0')
S.mask_where('LwLatSens','grid.mask!=2.0')
S.mask_where('PmE','grid.mask!=2.0')
S.mask_where('salt_flux','grid.mask!=2.0')
var_dict=S.var_dict['SW'].copy()
S.create_field('self.SW+self.LwLatSens','heat_flux',var_dict=var_dict)
S.mask_where('heat_flux','grid.mask!=2.0')

f=nc.Dataset('/net3/mjh/models/CM2G/Vertical_coordinate.nc')
R=f.variables['R'][:]
nrho=R.shape[0]
Rb=np.zeros((nrho+1))

Rb[0]=0.0
for i in np.arange(1,nrho):
    Rb[i]=0.5*(R[i]+R[i-1])
    
Rb[nrho]=2.0*Rb[nrho-1]

print Rb

S.sfc_buoyancy_production(sst='temp',sss='salt',heat_flux='heat_flux',fw_flux='PmE',salt_flux='salt_flux',p_ref=2.e7,rho_bounds=Rb)

sout=np.mean(S.thermal_buoyancy_flux_d,axis=0)/1.e6
print sout.max(),sout.min()
plt.plot(sout,R,color='b')

sout=np.mean(S.fw_buoyancy_flux_d,axis=0)/1.e6
print sout.max(),sout.min()

plt.plot(sout,R,color='r')

plt.ylim(1038,1024)
plt.show()


