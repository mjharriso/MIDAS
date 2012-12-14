from midas import *


grid_flk=generic_rectgrid('http://thredds.jpl.nasa.gov/thredds/dodsC/OceanWinds/2000_CCMP_MEASURES_ATLAS_L4_OW_L3_0_WIND_VECTORS_FLK_v11l30flk.nc',var='uwnd')
#S=state('http://thredds.jpl.nasa.gov/thredds/dodsC/OceanWinds/2000_CCMP_MEASURES_ATLAS_L4_OW_L3_0_WIND_VECTORS_FLK_v11l30flk.nc',grid=grid_flk,fields=['uwnd','vwnd'],default_calendar='julian',date_bounds=[datetime(2000,1,1),datetime(2000,1,2)])
S=state('http://thredds.jpl.nasa.gov/thredds/dodsC/OceanWinds/2000_CCMP_MEASURES_ATLAS_L4_OW_L3_0_WIND_VECTORS_FLK_v11l30flk.nc',grid=grid_flk,fields=['uwnd','vwnd'],default_calendar='julian',time_indices=np.arange(0,1))


sgrid=supergrid(file='/net3/mjh/models/GIS_0125/grids/supergrid.nc')
grid = ocean_rectgrid(supergrid=sgrid)
f=nc.Dataset('/net3/mjh/models/GIS_0125/ocean/topog_v1.nc')
topog=f.variables['depth'][:]

grid.wet[-1,:]=0.0

S.fill_interior('uwnd')
S.fill_interior('vwnd')

S.var_dict['uwnd']['Ztype']='Fixed'
S.var_dict['vwnd']['Ztype']='Fixed'

T=S.horiz_interp(field_x='uwnd',field_y='vwnd',target=grid,src_modulo=True,method=1)


T.rename_field('uwnd','u10')
T.rename_field('vwnd','v10')

T.write_nc('flk_regrid.nc',['u10','v10'])

