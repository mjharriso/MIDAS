# IPython log file

from midas import *
grid_obs=generic_grid('/net2/mjh/data/Steele/ptemp_salt.nc',var='PTEMP')
O=state('/net2/mjh/data/Steele/ptemp_salt.nc',fields=['PTEMP'],z_indices=np.arange(0,1),grid=grid_obs,default_calendar='noleap')
grid=gold_grid('/net3/mjh/models/CM2G/ocean_geometry.nc')
O.rename_field('PTEMP','SST')
Obs=O.horiz_interp('SST',target=grid,src_modulo=True)
C=state(MFpath='/archive/mjh/CM2G63L/CM2G_1860_Control/gfdl.ncrc2-default-prod-openmp/pp/ocean/ts/monthly/5yr/ocean.*.SST.nc',grid=grid,fields=['SST'],date_bounds=[datetime(600,1,1),datetime(799,12,31)])
HMIX20=state(MFpath='/archive/mjh/CM2G63L/CM2G_1860_Control_HMIX20/gfdl.ncrc2-default-prod-openmp/pp/ocean/ts/monthly/5yr/ocean.*.SST.nc',grid=grid,fields=['SST'])
HMIX40=state(MFpath='/archive/mjh/CM2G63L/CM2G_1860_Control_HMIX40/gfdl.ncrc2-default-prod-openmp/pp/ocean/ts/monthly/5yr/ocean.*.SST.nc',grid=grid,fields=['SST'])
C.var_dict['SST']['Ztype']='Fixed'
HMIX20.var_dict['SST']['Ztype']='Fixed'
HMIX40.var_dict['SST']['Ztype']='Fixed'
Obs.monthly_avg('SST',vol_weight=False)
C.monthly_avg('SST',vol_weight=False)
HMIX20.monthly_avg('SST',vol_weight=False)
HMIX40.monthly_avg('SST',vol_weight=False)
sst_range_o = np.max(Obs.SST_monthly,axis=0)-np.min(Obs.SST_monthly,axis=0)
sst_range_c = np.max(C.SST_monthly,axis=0)-np.min(C.SST_monthly,axis=0)
sst_range_20 = np.max(HMIX20.SST_monthly,axis=0)-np.min(HMIX20.SST_monthly,axis=0)
sst_range_40 = np.max(HMIX40.SST_monthly,axis=0)-np.min(HMIX40.SST_monthly,axis=0)
xax=grid.x_T
yax=grid.y_T

fig=plt.figure(1)
ax1=fig.add_subplot(211)


cf=ax1.contourf(xax,yax,sq(sst_range_o),np.arange(.25,10.25,.25),extend='both');ax1.contour(xax,yax,sq(sst_range_o),np.arange(0,10,1),colors='k');plt.colorbar(cf,ax=ax1);ax1.set_xlim(-240,-69);ax1.set_ylim(-40,20);ax1.grid();ax1.set_ylabel('Latitude (degN)');ax1.set_title('Observed (WOA) SST max-min seasonal cycle')
ax2=fig.add_subplot(212)
cf=ax2.contourf(xax,yax,sq(sst_range_c),np.arange(.25,10.25,.25),extend='both');ax2.contour(xax,yax,sq(sst_range_c),np.arange(0,10,1),colors='k');plt.colorbar(cf,ax=ax2);ax2.set_xlim(-240,-69);ax2.set_ylim(-40,20);ax2.grid();ax2.set_ylabel('Latitude (degN)');ax2.set_xlabel('Longitude (degE)');ax2.set_title('CM2G HMIN=2m SST max-min seasonal cycle')
plt.savefig('sst_ann_range_cm2g.png')

fig=plt.figure(2)
ax1=fig.add_subplot(211)


cf=ax1.contourf(xax,yax,sq(sst_range_o),np.arange(.25,10.25,.25),extend='both');ax1.contour(xax,yax,sq(sst_range_o),np.arange(0,10,1),colors='k');plt.colorbar(cf,ax=ax1);ax1.set_xlim(-240,-69);ax1.set_ylim(-40,20);ax1.grid();ax1.set_ylabel('Latitude (degN)');ax1.set_title('Observed (WOA) SST max-min seasonal cycle')
ax2=fig.add_subplot(212)
cf=ax2.contourf(xax,yax,sq(sst_range_20),np.arange(.25,10.25,.25),extend='both');ax2.contour(xax,yax,sq(sst_range_20),np.arange(0,10,1),colors='k');plt.colorbar(cf,ax=ax2);ax2.set_xlim(-240,-69);ax2.set_ylim(-40,20);ax2.grid();ax2.set_ylabel('Latitude (degN)');ax2.set_xlabel('Longitude (degE)');ax2.set_title('CM2G HMIN=20m SST max-min seasonal cycle')
plt.savefig('sst_ann_range_hmix20.png')

fig=plt.figure(3)
ax1=fig.add_subplot(211)


cf=ax1.contourf(xax,yax,sq(sst_range_o),np.arange(.25,10.25,.25),extend='both');ax1.contour(xax,yax,sq(sst_range_o),np.arange(0,10,1),colors='k');plt.colorbar(cf,ax=ax1);ax1.set_xlim(-240,-69);ax1.set_ylim(-40,20);ax1.grid();ax1.set_ylabel('Latitude (degN)');ax1.set_title('Observed (WOA) SST max-min seasonal cycle')
ax2=fig.add_subplot(212)
cf=ax2.contourf(xax,yax,sq(sst_range_40),np.arange(.25,10.25,.25),extend='both');ax2.contour(xax,yax,sq(sst_range_40),np.arange(0,10,1),colors='k');plt.colorbar(cf,ax=ax2);ax2.set_xlim(-240,-69);ax2.set_ylim(-40,20);ax2.grid();ax2.set_ylabel('Latitude (degN)');ax2.set_xlabel('Longitude (degE)');ax2.set_title('CM2G HMIN=40m SST max-min seasonal cycle')
plt.savefig('sst_ann_range_hmix40.png')

plt.show()
