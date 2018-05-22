from midas.rectgrid import *

fpath='https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/nodc.woa98/temperat/seasonal/otemp.gp1deg.nc'
grid=quadmesh(fpath,var='otemp',cyclic=True)
IO = grid.geo_region(x=(30,120.),y=(-30.,25.))
S=state(fpath,grid=grid,geo_region=IO,fields=['otemp'],default_calendar='noleap')
S.volume_integral('otemp','XYZ',normalize=False)
rho=1.e3;cp=3989.
print('Total heat content (J) = ',rho*cp*sq(S.otemp_xyzint))
S.time_avg('otemp')
S.volume_integral('otemp_tav','XYZ',normalize=False)
print('Total heat content (J) = ',rho*cp*sq(S.otemp_tav_xyzint[0]))

