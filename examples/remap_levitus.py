from midas import *


grid=gold_grid('/net3/mjh/models/CM2G/ocean_geometry.nc')
grid_obs=generic_grid('/net2/mjh/data/Steele/ptemp_salt.nc',var='PTEMP')
S=state(path='/net2/mjh/data/Steele/ptemp_salt.nc',grid=grid_obs,geo_region=None,fields=['PTEMP','SALT'],date_bounds=[datetime(1900,1,1,0,0,0),datetime(1900,1,30,0,0,0)],default_calendar='noleap')


fvgrid=nc.Dataset('/net3/mjh/models/CM2G/Vertical_coordinate.nc')
R=fvgrid.variables['R'][:]

nkml=2;nkbl=2;min_depth=10.0;p_ref=2.e7;hml=5.0;fit_target=True

T=S.horiz_interp('PTEMP',target=grid,src_modulo=True,method='bilinear')
T=S.horiz_interp('SALT',target=grid,src_modulo=True,method='bilinear',PrevState=T)
T.remap_Z_to_layers('PTEMP','SALT',R,p_ref,grid.wet,nkml,nkbl,hml,fit_target)


print '... Done remapping '

T.add_interface_bounds('temp_remap')


xax=np.squeeze(T.grid.x_T_bounds[0,:])
yax=np.squeeze(T.grid.y_T_bounds[:,0])
zint=T.var_dict['temp_remap']['z_interfaces_ns']

SHOW_FIGS = 0
PICKLE_IT = 0
WRITE_IT = 1

print T.temp_remap.shape
print T.var_dict['temp_remap']['dz'].shape

if SHOW_FIGS > 0:
    fig1=plt.figure(1)
    plt.clf()
    sout=np.squeeze(T.temp_remap[:,:,100])
    hout=np.squeeze(T.var_dict['temp_remap']['dz'][:,:,100])
    sout=np.ma.masked_where(hout<1.e-3,sout)
    cf=plt.pcolormesh(np.arange(0,len(yax)),np.squeeze(zint[:,:,100]),sout,vmin=-2.0,vmax=22.)
    fig1.colorbar(cf)

    for i in np.arange(0,zint.shape[1]):
        plt.plot(np.arange(0,len(yax)),np.squeeze(zint[0,i,:,100]),color='k')

    fig2=plt.figure(2)
    yax=T.var_dict['PTEMP']['yax_data']
    zdat = T.var_dict['PTEMP']['z_interfaces']
    sout=np.squeeze(T.PTEMP[0,:,:,100])
    cf=plt.pcolormesh(np.arange(0,len(yax)),np.squeeze(zdat[:,:,100]),sout,vmin=-2.0,vmax=22.)
    fig2.colorbar(cf)

    plt.show()

if PICKLE_IT == 1:
    print '... Pickling '

    pickle.dump(T,open('./IC.p',"wb"))

if WRITE_IT == 1:
    print '... Write NC '
    T.var_dict['temp_remap']['units']='degC'
    T.var_dict['salt_remap']['units']='psu'
    T.write_nc('IC.nc',['temp_remap','salt_remap'])

