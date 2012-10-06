#

import cdms2 as cdms
import numpy as np
import MV2 as mv

#==============================
# Set defaults for Netcdf3 output
#==============================

cdms.setNetcdfDeflateFlag(0)
cdms.setNetcdfDeflateLevelFlag(0)
cdms.setNetcdfShuffleFlag(0)


f=cdms.open('ocean_geometry.nc')
f_hgrid=cdms.open('mosaic/ocean_hgrid.nc')

x_t=np.array(f('geolon'))
[nj,ni]=x_t.shape
y_t=np.array(f('geolat'))
x_c=np.array(f('geolonb'))
ibeg = 3*ni / 4
for i in range (ibeg,ni):
    x_c[nj-1,i]=x_c[nj-1,i] + 360.

y_c=mv.array(f('geolatb'))

area = np.array(f_hgrid('area'))

area_refined=np.zeros((2*nj+1,2*ni+1))

for j in range (0,2*nj):
    area_refined[j,0:2*ni]=area[j,:]
    area_refined[j,2*ni]=area[j,0]

for i in range(0,2*ni):
  area_refined[2*nj,i]=area_refined[2*nj-1,2*ni-i]

area_refined[2*nj,2*ni]=area_refined[2*nj,0]

t_area=np.zeros((nj,ni))
e_area=np.zeros((nj,ni))
c_area=np.zeros((nj,ni))
n_area=np.zeros((nj,ni))

for j in range (0,nj):
  for i in range (0,ni):
    t_area[j,i]=np.sum(area_refined[2*(j-1)+2:2*(j-1)+4,2*(i-1)+2:2*(i-1)+4])
    c_area[j,i]=np.sum(area_refined[2*(j-1)+3:2*(j-1)+5,2*(i-1)+3:2*(i-1)+5])
    e_area[j,i]=np.sum(area_refined[2*(j-1)+2:2*(j-1)+4,2*(i-1)+3:2*(i-1)+4])
    n_area[j,i]=np.sum(area_refined[2*(j-1)+3:2*(j-1)+5,2*(i-1)+2:2*(i-1)+5])

# T and C arrays have ghost points
x_T=np.zeros((nj+2,ni+2))
y_T=np.zeros((nj+2,ni+2))
x_C=np.zeros((nj+2,ni+2))
y_C=np.zeros((nj+2,ni+2))

for j in range (1,nj+1):
  for i in range (1,ni+1):
    x_T[j,i]=x_t[j-1,i-1]  # shift center points by 1 index value to NE
    x_C[j,i]=x_c[j-1,i-1] 
    y_T[j,i]=y_t[j-1,i-1] 
    y_C[j,i]=y_c[j-1,i-1] 


# Fill ghost points

for j in range (1,nj+1):
  x_T[j,0]=x_t[j-1,ni-1]-360.0
  x_T[j,ni+1]=x_t[j-1,0]+360.
  x_C[j,0]=x_c[j-1,ni-1]
  x_C[j,ni+1]=x_c[j-1,0]
  y_T[j,0]=y_t[j-1,ni-1]
  y_T[j,ni+1]=y_t[j-1,0]
  y_C[j,0]=y_c[j-1,ni-1]
  y_C[j,ni+1]=y_c[j-1,0]

x_T[0,:]=x_T[1,:]
x_C[0,:]=x_C[1,:]
y_T[0,:]=y_T[1,:]-0.5*(y_T[2,:]-y_T[1,:])
y_C[0,:]=y_C[1,:]-0.5*(y_C[2,:]-y_C[1,:])

for j in range (1,nj):
  x_C[j,0]=x_C[j,0]-360.
  x_C[j,ni+1]=x_C[j,ni+1]+360.

x_C[0,0]=x_C[1,0]

for i in range(0,ni+2):
  x_T[nj+1,i]=x_T[nj,ni+1-i]
  x_C[nj+1,i]=x_C[nj-1,ni+1-i]
  y_T[nj+1,i]=y_T[nj,ni+1-i]
  y_C[nj+1,i]=y_C[nj-1,ni+1-i]

for i in range (1,ni+1):
  y_C[0,i]=y_t[0,i-1]-(y_c[0,i-1]-y_t[0,i-1])

for j in range (1,nj+1):
  y_C[j,0]=y_c[j-1,ni-1]

y_C[0,0]=y_C[0,1]

# Face-centered points (no ghost points)

x_e = np.zeros((nj,ni))
x_n = np.zeros((nj,ni))
y_e = np.zeros((nj,ni))
y_n = np.zeros((nj,ni))

for j in range (0,nj):
  for i in range (0,ni):
    x_e[j,i]=0.5*(x_C[j,i+1]+x_C[j+1,i+1])
    y_e[j,i]=0.5*(y_C[j,i+1]+y_C[j+1,i+1])
    x_n[j,i]=0.5*(x_C[j+1,i]+x_C[j+1,i+1])
    y_n[j,i]=0.5*(y_C[j+1,i]+y_C[j+1,i+1])

x_t_bounds=np.zeros((4,nj,ni))
y_t_bounds=np.zeros((4,nj,ni))
x_e_bounds=np.zeros((4,nj,ni))
y_e_bounds=np.zeros((4,nj,ni))
x_c_bounds=np.zeros((4,nj,ni))
y_c_bounds=np.zeros((4,nj,ni))
x_n_bounds=np.zeros((4,nj,ni))
y_n_bounds=np.zeros((4,nj,ni))

for j in range (0,nj):
  for i in range (0,ni):
    x_t_bounds[0,j,i]=x_C[j,i]
    x_t_bounds[1,j,i]=x_C[j,i+1]
    x_t_bounds[2,j,i]=x_C[j+1,i+1]
    x_t_bounds[3,j,i]=x_C[j+1,i]
    y_t_bounds[0,j,i]=y_C[j,i]
    y_t_bounds[1,j,i]=y_C[j,i+1]
    y_t_bounds[2,j,i]=y_C[j+1,i+1]
    y_t_bounds[3,j,i]=y_C[j+1,i]
    x_c_bounds[0,j,i]=x_T[j+1,i+1]
    x_c_bounds[1,j,i]=x_T[j+1,i+2]
    x_c_bounds[2,j,i]=x_T[j+2,i+2]
    x_c_bounds[3,j,i]=x_T[j+2,i+1]
    y_c_bounds[0,j,i]=y_T[j+1,i+1]
    y_c_bounds[1,j,i]=y_T[j+1,i+2]
    y_c_bounds[2,j,i]=y_T[j+2,i+2]
    y_c_bounds[3,j,i]=y_T[j+2,i+1]
    x_e_bounds[0,j,i]=0.5*(x_C[j,i]+x_C[j,i+1])
    x_e_bounds[1,j,i]=0.5*(x_C[j,i+1]+x_C[j,i+2])
    x_e_bounds[2,j,i]=0.5*(x_C[j+1,i+1]+x_C[j+1,i+2])
    x_e_bounds[3,j,i]=0.5*(x_C[j+1,i]+x_C[j+1,i+1])
    y_e_bounds[0,j,i]=0.5*(y_C[j,i]+y_C[j,i+1])
    y_e_bounds[1,j,i]=0.5*(y_C[j,i+1]+y_C[j,i+2])
    y_e_bounds[2,j,i]=0.5*(y_C[j+1,i+1]+y_C[j+1,i+2])
    y_e_bounds[3,j,i]=0.5*(y_C[j+1,i]+y_C[j+1,i+1])
    x_n_bounds[0,j,i]=0.5*(x_T[j+1,i+1]+x_T[j+1,i])
    x_n_bounds[1,j,i]=0.5*(x_T[j+1,i+1]+x_T[j+1,i+2])
    x_n_bounds[2,j,i]=0.5*(x_T[j+2,i+1]+x_T[j+2,i+2])
    x_n_bounds[3,j,i]=0.5*(x_T[j+2,i+1]+x_T[j+2,i])
    y_n_bounds[0,j,i]=0.5*(y_T[j+1,i+1]+y_T[j+1,i])
    y_n_bounds[1,j,i]=0.5*(y_T[j+1,i+1]+y_T[j+1,i+2])
    y_n_bounds[2,j,i]=0.5*(y_T[j+2,i+1]+y_T[j+2,i+2])
    y_n_bounds[3,j,i]=0.5*(y_T[j+2,i+1]+y_T[j+2,i])


##### Kludge Fix for tripolar (mjh)
x_e_bounds[3,nj-1,0]=x_e_bounds[0,nj-1,0]
x_e_bounds[2,nj-1,ni-1]=x_e_bounds[3,nj-1,ni-1]

x_n_bounds[2,nj-1,0:ni/4]=x_n_bounds[2,nj-1,0:ni/4]-360.
x_n_bounds[3,nj-1,0:ni/4]=x_n_bounds[3,nj-1,0:ni/4]-360.

x_n_bounds[2,nj-1,3*ni/4:]=x_n_bounds[2,nj-1,3*ni/4:]+360.
x_n_bounds[3,nj-1,3*ni/4:]=x_n_bounds[3,nj-1,3*ni/4:]+360.

x_t_bounds[3,nj-1,0]=x_t_bounds[0,nj-1,0]


x_T_ = mv.array(f('geolon'))
y_T_ = mv.array(f('geolat'))
x_C_ = mv.array(f('geolonb'))
y_C_ = mv.array(f('geolatb'))

x_T=x_T[1:nj+1,1:ni+1]
y_T=y_T[1:nj+1,1:ni+1]
x_C=x_C[1:nj+1,1:ni+1]
y_C=y_C[1:nj+1,1:ni+1]


xh=x_T_.getAxis(1)
yh=x_T_.getAxis(0)
xq=x_C_.getAxis(1)
yq=x_C_.getAxis(0)

x_T=mv.array(x_T)
y_T=mv.array(y_T)
x_C=mv.array(x_C)
y_C=mv.array(y_C)

setattr(x_T,'id','x_T')
setattr(x_T,'name','x_T')
setattr(y_T,'id','y_T')
setattr(y_T,'name','y_T')
x_T.setAxisList((yh,xh))
y_T.setAxisList((yh,xh))

T_area=mv.array(t_area)
setattr(T_area,'id','area_T')
setattr(T_area,'name','area_T')
T_area.setAxisList((yh,xh))
T_area.setAxisList((yh,xh))

setattr(x_C,'id','x_C')
setattr(x_C,'name','x_C')
setattr(y_C,'id','y_C')
setattr(y_C,'name','y_C')
x_C.setAxisList((yq,xq))
y_C.setAxisList((yq,xq))

C_area=mv.array(c_area)
setattr(C_area,'id','area_C')
setattr(C_area,'name','area_C')
C_area.setAxisList((yq,xq))
C_area.setAxisList((yq,xq))

x_E=mv.array(x_e)
y_E=mv.array(y_e)
setattr(x_E,'id','x_E')
setattr(x_E,'name','x_E')
setattr(y_E,'id','y_E')
setattr(y_E,'name','y_E')
x_E.setAxisList((yh,xq))
y_E.setAxisList((yh,xq))

E_area=mv.array(e_area)
setattr(E_area,'id','area_E')
setattr(E_area,'name','area_E')
E_area.setAxisList((yh,xq))
E_area.setAxisList((yh,xq))

x_N=mv.array(x_n)
y_N=mv.array(y_n)
setattr(x_N,'id','x_N')
setattr(x_N,'name','x_N')
setattr(y_N,'id','y_N')
setattr(y_N,'name','y_N')
x_N.setAxisList((yq,xh))
y_N.setAxisList((yq,xh))

N_area=mv.array(n_area)
setattr(N_area,'id','area_N')
setattr(N_area,'name','area_N')
N_area.setAxisList((yq,xh))
N_area.setAxisList((yq,xh))

vertex=cdms.createAxis([0,1,2,3])
vertex.id='vertex'

x_T_bounds=mv.array(x_t_bounds)    
y_T_bounds=mv.array(y_t_bounds)    
setattr(x_T_bounds,'id','x_T_bounds')
setattr(x_T_bounds,'name','x_T_bounds')
setattr(y_T_bounds,'id','y_T_bounds')
setattr(y_T_bounds,'name','y_T_bounds')
x_T_bounds.setAxisList((vertex,yh,xh))
y_T_bounds.setAxisList((vertex,yh,xh))


x_E_bounds=mv.array(x_e_bounds)    
y_E_bounds=mv.array(y_e_bounds)    
setattr(x_E_bounds,'id','x_E_bounds')
setattr(x_E_bounds,'name','x_E_bounds')
setattr(y_E_bounds,'id','y_E_bounds')
setattr(y_E_bounds,'name','y_E_bounds')
x_E_bounds.setAxisList((vertex,yh,xq))
y_E_bounds.setAxisList((vertex,yh,xq))

x_C_bounds=mv.array(x_c_bounds)    
y_C_bounds=mv.array(y_c_bounds)    
setattr(x_C_bounds,'id','x_C_bounds')
setattr(x_C_bounds,'name','x_C_bounds')
setattr(y_C_bounds,'id','y_C_bounds')
setattr(y_C_bounds,'name','y_C_bounds')
x_C_bounds.setAxisList((vertex,yq,xq))
y_C_bounds.setAxisList((vertex,yq,xq))

x_N_bounds=mv.array(x_n_bounds)    
y_N_bounds=mv.array(y_n_bounds)    
setattr(x_N_bounds,'id','x_N_bounds')
setattr(x_N_bounds,'name','x_N_bounds')
setattr(y_N_bounds,'id','y_N_bounds')
setattr(y_N_bounds,'name','y_N_bounds')
x_N_bounds.setAxisList((vertex,yq,xh))
y_N_bounds.setAxisList((vertex,yq,xh))

g=cdms.open('gold_cell_bounds.nc','w')

g.write(x_T)
g.write(y_T)
g.write(x_C)
g.write(y_C)
g.write(x_E)
g.write(y_E)
g.write(x_N)
g.write(y_N)
g.write(T_area)
g.write(E_area)
g.write(C_area)
g.write(N_area)

g.write(x_T_bounds)
g.write(x_E_bounds)
g.write(x_C_bounds)
g.write(x_N_bounds)
g.write(y_T_bounds)
g.write(y_E_bounds)
g.write(y_C_bounds)
g.write(y_N_bounds)

g.close()

