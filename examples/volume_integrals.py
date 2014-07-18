from midas.rectgrid import *


grid=quadmesh('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_annual_1deg.nc',var='t_an',cyclic=True)
IO = grid.geo_region(x=(30,120.),y=(-30.,25.))
S=state('http://data.nodc.noaa.gov/thredds/dodsC/woa/WOA09/NetCDFdata/temperature_monthly_1deg.nc',grid=grid,geo_region=IO,fields=['t_an'],default_calendar='noleap')
S.volume_integral('t_an','XYZ',normalize=False)
rho=1.e3;cp=3989.
print 'Total heat content (J) = ',rho*cp*sq(S.t_an_xyzint)
S.time_avg('t_an')
S.volume_integral('t_an_tav','XYZ',normalize=False)
print 'Total heat content (J) = ',rho*cp*sq(S.t_an_tav_xyzint[0])

