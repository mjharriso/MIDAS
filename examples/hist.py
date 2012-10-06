from midas import *



print """\n
         Produce a volume-weighted histogram
         from GOLD layer output\n"""

grid=Grid('ocean_geometry.nc')
dimsiz=query_size('19080101.ocean_month.nc','temp')
gregion=grid.geo_region(x=(-220,-140),y=(-60,-30))
S=State('19080101.ocean_month.nc',grid,region=gregion)

b = S.temp
vol = S.grid.Ah*S.var_dict['temp']['dz']
wt = vol
result = plt.hist(b,bins=np.arange(-2,33),weights=wt,norm=True)
plt.show()
