######################################################################
#
#  This file is copied from ROMSTOOLS
#
#  Build ROMS configuration for the Baroclinic Jet test case
#
#  Create grid file and initial file
# 
#  Further Information:  
#  http://www.brest.ird.fr/Roms_tools/
#  
# 
#
#  ROMSTOOLS is free software you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation either version 2 of the License,
#  or (at your option) any later version.
#
#  ROMSTOOLS is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
#  Xavier Capet and Patrick Marchesiello - Nov 2012
#
#  Yves Soufflet for the translation to Python -Jan 2013
#
######################################################################
#clear all
#close all
##################### USERS DEFINED VARIABLES ########################
#
#  Title
#
title_case='JET'
#
#  Names
#
parent_grd ='jet_grd.nc'
parent_ini ='jet_ini.nc'
parent_clm ='jet_clm.nc'
#
# Domain size
#
xmax = 250e3      # Domain half length LAT
ymax = 1000e3     # Domain half length LON
H0   = 4000       # Total Depth
#
# Horizontal grid resolution
#
dx=5.e3
#
# Vertical grid parameters
#
N=30
theta_s=5
theta_b=0
hc=100
#
# Jet parameters
#
jetwidth=100e3    # jet width
lat=45            # latitude of jet center
bplane=1          # bplane=1 or 0 (beta-plane or f-plane)
perturb=1         # add perturbation to zonal jet
#
#   subsurface Phillips mode properties
rhomaxn=27.7573   # northern max density
rhomaxs=27.75     # southern max density
vbmaxn=9.8e-6     # northern background stratification
vbmaxs=9.8e-6     # southern background stratification
zs1=-1000;dzs=700 # northern profile parameter
zn1=-400 ;dzn=300 # southern profile parameter
drhos=1.4         # phillips tends to dominate
#
#   surface Charney mode properties
drhosfs=1.5       # anomaly in the south
drhosfn=0.0       # anomaly in the north
z0=-300           # vertical penetration
drhosf=0.00       # weakens Charney mode (canceled for 1.05, max for 0)
z0p=-110
#
rho0=1024.4       # Mean ocean density
g=9.81            # Gravity acceleration
R=6367442.76      # Earth radius
#
################### END USERS DEFINED VARIABLES #######################
#
#
# Horizontal Grid
#

import numpy as np
from Preprocessing_tools import zlevs

dy=dx
x=np.arange(-xmax-dx/2,xmax+dx,dx)
y=np.arange(-ymax-dy/2,ymax+dy,dy)
X,Y=np.meshgrid(x,y)

(Mp,Lp)=X.shape;

M=Mp-1 ;L=Lp-1
Mm=Mp-2; Lm=Lp-2
pm=np.ones((Mp,Lp))/dx

r=np.ones(N)

pn=np.ones((Mp,Lp))/dy

xr=np.reshape(np.outer(r,X),(N,Mp,Lp))
yr=np.reshape(np.outer(r,Y),(N,Mp,Lp))

#
# Vertical grid
#
zeta=np.zeros((Mp,Lp))

h=H0+0*X

zw = zlevs.zlevs(h,zeta,theta_s,theta_b,hc,N,'w')
zr = zlevs.zlevs(h,zeta,theta_s,theta_b,hc,N,'r')
zr1d  = zlevs.zlevs(H0,0,theta_s,theta_b,hc,N,'r')
zr1d=zr1d.squeeze()
#
# Coriolis term (beta plane)
#
pi=np.pi
deg2rad=pi/180
omega=2*pi/(24*3600)
f0=2*omega*np.sin(deg2rad*lat)
if bplane:
    beta=2*omega/R*np.cos(deg2rad*lat)
else:
    beta=0

f=f0+beta*Y
#
#--------------------------------------------------------
# Make northern and southern density profiles
#--------------------------------------------------------
#
# First, get a gentle stratification that does not depend 
# on jet side. It is there to ensure static stability. 
#

rhon=np.zeros(N)
rhos=rhon
for k in range(N):
    zv = zr1d[k]
    rhon[k]=rhomaxn-vbmaxn*(zv+4000)
    rhos[k]=rhomaxs-vbmaxs*(zv+4000)

#
# Second, get main north/south contrast with a distorded 
# tanh shape of variable amplitude. Distorsion of the tanh 
# is done with the asym function that increases 
# stratification in the upper ocean 
#
dzs2=1.3*dzs
#zrm=asym_strat(zr1d,zs1,dzs2)
for k in range(N):
    zv = (zr1d[k]-zs1)/np.sqrt((1+0.5*(zr1d[k]-zs1-abs(zr1d[k]-zs1))**2)/dzs2**2)
    rhos[k]=rhos[k]-drhos*(0.5+0.5*np.tanh((zv-zs1)/dzs))

dzn2=1.3*dzn

#zrm=asym_strat(zr1d,zn1,dzn2)
zrm=np.zeros(N)
for k in range(N):
    zrm[k]= (zr1d[k]-zn1)/np.sqrt((1+0.5*(zr1d[k]-zn1-abs(zr1d[k]-zn1))**2)/dzn2**2)

drhon=(rhon[-1]-rhos[-1])/(0.5+0.5*np.tanh((zrm[-1]-zn1)/dzn))
for k in range(N):
    zv = zrm[k]
    rhon[k]=rhon[k]-drhon*(0.5+0.5*np.tanh((zv-zn1)/dzn))

for k in range(N):
    zv = zr1d[k]
    rhon[k]=rhon[k]-drhosfn*0.5*(1+np.tanh((zv-z0)/abs(z0)))/np.tanh(1)
    rhos[k]=rhos[k]-drhosfs*0.5*(1+np.tanh((zv-z0)/abs(z0)))/np.tanh(1)
    rhon[k]=rhon[k]-drhosf*(np.exp((zv-z0p)/abs(z0p)))/(np.exp(+1)) ## 0 for setup1
    rhos[k]=rhos[k]-drhosf*(np.exp((zv-z0p)/abs(z0p)))/(np.exp(+1)) ## 0 for setup1
#
#--------------------------------------------------------
# Fit southern and northern profiles using linear and
# sin^2 functions. fthy is the shape integral in y
#--------------------------------------------------------
#
ywidth = jetwidth/2
ydist  = ymax-ywidth
j0 = 0
j1 = j0+int(ydist/dy)
j2 = j1+int(ywidth/dy)
j3 = j2
j4 = j3+int(ywidth/dy)
j5 = Mp
jmidjet=j2
tmpy=np.zeros(Mp)

yfac1=pi*0.5/(j2-j1)
yfac3=pi*0.5/(j4-j3)

for j in range(j0,j1):
      tmpy[j]=0.

for j in range(j4,j5):
    tmpy[j]=0.

for j in range(j1+1,j2):
    tmpy[j]=np.sin((j-j1)*yfac1)**2

for j in range(j2+1,j3):
    tmpy[j]=1.

for j in range(j3+1,j4):
    tmpy[j]=np.cos((j-j3)*yfac3)**2


# Now integrate in y to find fthy.
# tmpy is at rho (u) points, fthy at psi (v) points.
fthy=np.zeros(Mp)

rho=np.zeros((Mp,N))
for j in range(1,Mp):
  fthy[j]=fthy[j-1]+0.5*(tmpy[j]+tmpy[j-1])

yfac3=fthy[-1]
fthy=fthy/yfac3

# fill in rho
for k in range(N):
    for j in range(Mp):
        rho[j,k]=rhos[k]*(1-fthy[j])+rhon[k]*fthy[j]
#        fsouth=1.-fthy[j]
#        fnorth=fthy[j]
        


t=np.zeros((N,Mp,Lp))

for k in range(N):
    for j in range(Mp):
        for i in range(Lp):
            t[k,j,i]=rho[j,k]

rho=t
#
#--------------------------------------------------------
# Add rho perturbation to zonal fields
#--------------------------------------------------------
#
# Perturbation function:
# find wavelength perturbation 2*pi/Kper near deformation 
# radius wavelength so that there is an integer number 
# of wave numbers in x direction
#
bvf=rho
if perturb:
    cff=g/rho0
    bvf[1:,:,:]=-cff*(rho[1:,:,:]-rho[:-1,:,:])/(zr[1:,:,:]-zr[:-1,:,:])
    mbvf=np.sqrt(bvf.max())
    Lrad=2*pi*mbvf*abs(2*z0)/f0
    Lper=2*xmax/np.floor(2*xmax/Lrad)
    Kper=2*pi/Lper      # perturbation wavenumber
    Cper=0.01           # perturbation magnitude
    t=t*(1+Cper*np.cos(Kper*xr)*np.exp(-(yr/jetwidth)**2)*np.exp(zr/abs(z0)))


#
#--------------------------------------------------------
# Compute dynamical fields
#--------------------------------------------------------
#
rho=1000+t
rhosurf=np.squeeze(rho[N-1,:,:])
rho_w=.5*(rho[:-1,:,:]+rho[1:,:,:])
dz_w=zr[1:,:,:]-zr[:-1,:,:]
#
#  Compute pressure field
#
pres=0*rho
pres[1,:,:]=-zr[1,:,:]*1.e4 #  initialize bottom pressure in Pa
for k in range(N-1):
  pres[k+1,:,:]=pres[k,:,:]-rho_w[k,:,:]*g*dz_w[k,:,:]

#
#  compute SSH and remove area averaged
#
ssh=np.squeeze(pres[N-1,:,:])/(rhosurf*g)
avgssh=ssh.mean()
ssh=ssh-avgssh
zeta=ssh
#
#  Compute geostrophic baroclinic velocities
#
p_u=0.5*(pres[:,:,:-1]+pres[:,:,1:])#rho2u_3d[pres]
p_v=0.5*(pres[:,:-1,:]+pres[:,1:,:])#rho2v_3d[pres]
px=np.zeros(pres.shape)
px[:,:,1:-1]=p_u[:,:,1:]-p_u[:,:,:-1]
px[:,:,0]=2*px[:,:,1]-px[:,:,2]
px[:,:,-1]=2.*px[:,:,-1]-px[:,:,-2]
py=np.zeros(pres.shape)
py[:,1:-2,:]=p_v[:,1:-1,:]-p_v[:,:-2,:]
py[:,0,:]=2.*py[:,1,:]-py[:,2,:]
py[:,-1,:]=2.*py[:,-2,:]-py[:,-3,:]


pn3d=np.reshape(np.outer(r,pn),(N,Mp,Lp))#tridim(pn,N) 
pm3d=np.reshape(np.outer(r,pm),(N,Mp,Lp))#tridim(pm,N) 
f3d=np.reshape(np.outer(r,f),(N,Mp,Lp))#tridim(pm,N)np.outer(f,r)#tridim(f,N) 
u_r=-pn3d*py/(rho0*f3d)
v_r= pm3d*px/(rho0*f3d)
u=0.5*(u_r[:,:,:-1]+u_r[:,:,1:])#rho2u_3d(u_r)
v=0.5*(v_r[:,:-1,:]+v_r[:,1:,:])#rho2v_3d(v_r)
#
# Compute barotropic velocities
#
zw=zlevs.zlevs(h,zeta,theta_s,theta_b,hc,N,'w')
zr=zlevs.zlevs(h,zeta,theta_s,theta_b,hc,N,'r')
dz=zw[1:,:,:]-zw[:1,:,:]
dzu=0.5*(dz[:,:,:-1]+dz[:,:,1:])
dzv=0.5*(dz[:,:-1,:]+dz[:,1:,:])
hu=sum(dzu*u)
hv=sum(dzv*v)
D_u=sum(dzu)
D_v=sum(dzv)
ubar=hu/D_u
vbar=hv/D_v
#
#--------------------------------------------------------
# Create and fill netcdf files
#--------------------------------------------------------
#
#  Create the grid file
#
(Mp,Lp)=zeta.shape
M=Mp-1
L=Lp-1
print(' ')
print('GRID DIMENSIONS Lm Mm : %d %d ',L-1,M-1)

#nc=netcdf(parent_grd, 'clobber')
#redef(nc)
#nc('xi_rho') = Lp
#nc('eta_rho') = Mp
#nc('xi_psi') = L
#nc('eta_psi') = M
#nc('one') = 1
import netCDF4
nc=netCDF4.Dataset('test.nc','w',clobber=True,format='NETCDF3_CLASSIC')
#nc._redef()
nc.createDimension('xi_rho',Lp)
nc.createDimension('eta_rho',Mp)
nc.createDimension('xi_psi',L)
nc.createDimension('eta_psi',M)
nc.createDimension('one',1)

el=nc.createVariable('el','f4',dimensions=('one',))
xl=nc.createVariable('xl','f4',dimensions=('one',))
yl=nc.createVariable('yl','f4',dimensions=('one',))
spherical=nc.createVariable('spherical','c',dimensions=('one',))
nc_h=nc.createVariable('h','f4',dimensions=('eta_rho','xi_rho'))
nc_f=nc.createVariable('f','f4',dimensions=('eta_rho','xi_rho'))
nc_pm=nc.createVariable('pm','f4',dimensions=('eta_rho','xi_rho'))
nc_pn=nc.createVariable('pn','f4',dimensions=('eta_rho','xi_rho'))

x_rho=nc.createVariable('x_rho','f4',dimensions=('eta_rho','xi_rho'))
y_rho=nc.createVariable('y_rho','f4',dimensions=('eta_rho','xi_rho'))
mask_rho=nc.createVariable('mask_rho','f4',dimensions=('eta_rho','xi_rho'))
#nc._endef()

xl[:]=x[-1]-x[0]
yl[:]=y[-1]-y[0]
spherical[:]='F'
nc_h[:]=H0
nc_f[:]=f
nc_pm[:]=1/dx
nc_pn[:]=1/dy
x_rho[:]=X
y_rho[:]=Y

mask_rho[:]=1+0*Y
nc.close()


#nc{'el'} = ncdouble('one')
#nc{'xl'} = ncdouble('one')
#nc{'spherical'} = ncchar('one')
#nc{'h'} = ncdouble('eta_rho', 'xi_rho')
#nc{'f'} = ncdouble('eta_rho', 'xi_rho')
#nc{'pm'} = ncdouble('eta_rho', 'xi_rho')
#nc{'pn'} = ncdouble('eta_rho', 'xi_rho')
#nc{'x_rho'} = ncdouble('eta_rho', 'xi_rho')
#nc{'y_rho'} = ncdouble('eta_rho', 'xi_rho')
#nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho')
#endef(nc)
#nc.title = ncchar(title_case)
#nc.title = title_case
#nc.date = ncchar(date)
#nc.date = date
#nc.type = ncchar('ROMS grid file')
#nc.type = 'ROMS grid file'
##
##  fill the grid file
##
#nc{'xl'}(:)=x(end)-x(1)
#nc{'el'}(:)=y(end)-y(1)
#nc{'spherical'}(:)='F'
#nc{'h'}(:)=H0
#nc{'f'}(:)=f
#nc{'pm'}(:)=1/dx
#nc{'pn'}(:)=1/dy
#nc{'x_rho'}(:)=X
#nc{'y_rho'}(:)=Y
#nc{'mask_rho'}(:)=1+0.*Y
#close(nc)
##
##  Create and fill the initial file
##
#create_inifile(parent_ini,parent_grd,title_case,...
#               theta_s,theta_b,hc,N,0,'clobber')
#nc=netcdf(parent_ini,'write')
#nc{'u'}(:) =  u 
#nc{'v'}(:) =  v 
#nc{'zeta'}(:) =  zeta 
#nc{'ubar'}(:) =  ubar 
#nc{'vbar'}(:) =  vbar 
#nc{'temp'}(:) =  t 
#close(nc)
##
##  Create and fill the climatology file
##
#create_climfile(parent_clm,parent_grd,title_case,...
#                theta_s,theta_b,hc,N,[25 75],100,'clobber')
#nc=netcdf(parent_clm,'write')
#nc{'u'}(1,:,:,:) =  u
#nc{'v'}(1,:,:,:) =  v
#nc{'zeta'}(1,:,:,:) =  zeta
#nc{'ubar'}(1,:,:,:) =  ubar
#nc{'vbar'}(1,:,:,:) =  vbar
#nc{'temp'}(1,:,:,:) =  t
#nc{'u'}(2,:,:,:) =  u
#nc{'v'}(2,:,:,:) =  v
#nc{'zeta'}(2,:,:,:) =  zeta
#nc{'ubar'}(2,:,:,:) =  ubar
#nc{'vbar'}(2,:,:,:) =  vbar
#nc{'temp'}(2,:,:,:) =  t
#close(nc)

from pylab import *
#
#--------------------------------------------------------
# Make a few plots
#--------------------------------------------------------
#
figure()
plot(rhos,zr1d,'r') #hold on 
plot(rhon,zr1d,'b')
plot(0.5*(rhos+rhon),zr1d,'k') #hold off
title('Density profiles')
#
X=X/1000 ;Y=Y/1000
xr=xr/1000; yr=yr/1000
zu=0.5*(zr[:,:,:-1]+zr[:,:,1:])
xu=0.5*(xr[:,:,:-1]+xr[:,:,1:])
yu=0.5*(yr[:,:,:-1]+yr[:,:,1:])
#
figure()

temp1=np.squeeze(u[N-1,:,:])
ur=np.zeros((Mp,Lp))
ur[:,1:-1]=0.5*(temp1[:,:-1]+temp1[:,1:])
ur[:,0]=ur[:,1]
ur[:,-1]=ur[:,-2]

temp2=np.squeeze(v[N-1,:,:])
vr=np.zeros((Mp,Lp))
vr[1:-1,:]=0.5*(temp2[:-1,:]+temp2[1:,:])
vr[0,:]=vr[1,:]
vr[-1,:]=vr[-2,:]


spd=np.sqrt(ur**2+vr**2)
pcolor(X,Y,spd)
#shading flat
#axis image
colorbar
#hold on
quiver(X,Y,ur,vr)
#hold off
title('Velocity')
#
figure()
(M,L)=Y.shape
imid=round(L/2)
contourf(np.squeeze(yr[:,:,imid]),np.squeeze(zr[:,:,imid]),np.squeeze(t[:,:,imid]),20)
colorbar
title('Temperature')
#
figure()
(M,L)=Y.shape
imid=round(L/2)
contourf(np.squeeze(yu[:,:,imid]),np.squeeze(zu[:,:,imid]),np.squeeze(u[:,:,imid]),10)
colorbar
title('Zonal velocity')
#
figure()
contourf(X,Y,zeta,10)
colorbar
#axis image
title('Surface height')
#
figure()
Xbar=0.5*(X[:-1,:]+X[1:,:])
Ybar=0.5*(Y[:-1,:]+Y[1:,:])
contourf(Xbar,Ybar,np.squeeze(v[N-1,:,:]),10)
colorbar
#axis image
title('Meridional Velocity (m/s)')
##
#us=squeeze(u[N-1,:,:])
#vs=squeeze(v[N-1,:,:])
#nc=netcdf(parent_grd)
#pm=nc{'pm'}(:)
#pn=nc{'pn'}(:)
#xr=1e-3*nc{'x_rho'}(:)
#yr=1e-3*nc{'y_rho'}(:)
#[xu,xv,xp]=rho2uvp(xr)
#[yu,yv,yp]=rho2uvp(yr)
#[fu,fv,fp]=rho2uvp(f)
#close(nc)
#[vort]=vorticity(us,vs,pm,pn)
#vort=vort./fp
#figure
#map=colormap(jet(20))
#map(10:11,:)=1*[1 1 1  1 1 1]
#colormap(map)
#contourf(xp,yp,vort,20) shading flat
#colorbar
#caxis([-0.8 0.8])
#axis image
#title('Vorticity/f')
#
