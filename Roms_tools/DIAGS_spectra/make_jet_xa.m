%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build ROMS configuration for the Baroclinic Jet test case
%
%  Create grid file and initial file
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Xavier Capet and Patrick Marchesiello - Nov 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
%  Title
%
title_case='JET';
%
%  Names
%
parent_grd ='jet_grd.nc';
parent_ini ='jet_ini.nc';
parent_clm ='jet_clm.nc';
%
% Domain size
%
xmax = 250e3;      % Domain half length LAT
ymax = 1000e3;     % Domain half length LON
H0   = 4000;       % Total Depth
%
% Horizontal grid resolution
%
dx=5.e3;
%
% Vertical grid parameters
%
N=30;
theta_s=5;
theta_b=0;
hc=100;
%
% Jet parameters
%
jetwidth=100e3;    % jet width
lat=45;            % latitude of jet center
bplane=1;          % bplane=1 or 0 (beta-plane or f-plane)
perturb=1;         % add perturbation to zonal jet
%
%   subsurface Phillips mode properties
rhomaxn=27.7573;   % northern max density
rhomaxs=27.75;     % southern max density
vbmaxn=9.8e-6;     % northern background stratification
vbmaxs=9.8e-6;     % southern background stratification
zs1=-1000;dzs=700; % northern profile parameter
zn1=-400; dzn=300; % southern profile parameter
drhos=1.4;         % phillips tends to dominate
%
%   surface Charney mode properties
drhosfs=1.5;       % anomaly in the south
drhosfn=0.0;       % anomaly in the north
z0=-300;           % vertical penetration
drhosf=0.00;       % weakens Charney mode (canceled for 1.05, max for 0)
z0p=-110;
%
rho0=1024.4;       % Mean ocean density
g=9.81;            % Gravity acceleration
R=6367442.76;      % Earth radius
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
%
% Horizontal Grid
%
dy=dx;
x=[-xmax-dx/2:dx:xmax+dx/2];
y=[-ymax-dy/2:dy:ymax+dy/2];
[X,Y]=meshgrid(x,y);
[Mp,Lp]=size(X);
M=Mp-1; L=Lp-1;
Mm=Mp-2; Lm=Lp-2;
pm=ones(Mp,Lp)./dx;
pn=ones(Mp,Lp)./dy;
xr=tridim(X,N);
yr=tridim(Y,N);

%
% Vertical grid
%
zeta=zeros(Mp,Lp);
h=H0+0*X;
zw = zlevs(h,zeta,theta_s,theta_b,hc,N,'w');
zr = zlevs(h,zeta,theta_s,theta_b,hc,N,'r');
zr1d  = zlevs(H0,0,theta_s,theta_b,hc,N,'r');

%
% Coriolis term (beta plane)
%
deg2rad=pi/180;
omega=2*pi/(24*3600);
f0=2*omega*sin(deg2rad*lat);
if bplane,
 beta=2*omega/R*cos(deg2rad*lat);
else
 beta=0;
end
f=f0+beta*Y;
%
%--------------------------------------------------------
% Make northern and southern density profiles
%--------------------------------------------------------
%
% First, get a gentle stratification that does not depend 
% on jet side. It is there to ensure static stability. 
%
for k=1:N
  zv = zr1d(k);
  rhon(k)=rhomaxn-vbmaxn*(zv+4000);
  rhos(k)=rhomaxs-vbmaxs*(zv+4000);
end
%
% Second, get main north/south contrast with a distorded 
% tanh shape of variable amplitude. Distorsion of the tanh 
% is done with the asym function that increases 
% stratification in the upper ocean 
%
dzs2=1.3*dzs;
zrm=(zr1d-zs1)/sqrt((1+0.5*power(zr1d-zs1-abs(zr1d-zs1),2))/power(dzs2,2))
%zeros(N,1);%asym_strat(zr1d,zs1,dzs2);
for k=1:N
    
  zv = zrm(k);
  rhos(k)=rhos(k)-drhos*(0.5+0.5*tanh((zv-zs1)/dzs));
end
dzn2=1.3*dzn;
zrm=(zr1d-zn1)/sqrt((1+0.5*power(zr1d-zn1-abs(zr1d-zn1),2))/power(dzn2,2))
%zeros(N,1);%asym_strat(zr1d,zn1,dzn2);
drhon=(rhon(N)-rhos(N))/(0.5+0.5*tanh((zrm(N)-zn1)/dzn));
for k=1:N
  zv = zrm(k);
  rhon(k)=rhon(k)-drhon*(0.5+0.5*tanh((zv-zn1)/dzn));
end
for k=1:N
  %zv = (zr1d(k)-zn1)/sqrt((1+0.5*power(zr1d(k)-zn1-abs(zr1d(k)-zn1),2))/power(dzn,2));
  zv = zr1d(k);
  rhon(k)=rhon(k)-drhosfn*0.5*(1+tanh((zv-z0)/abs(z0)))./tanh(1);
  rhos(k)=rhos(k)-drhosfs*0.5*(1+tanh((zv-z0)/abs(z0)))./tanh(1);
  rhon(k)=rhon(k)-drhosf*(exp((zv-z0p)/abs(z0p)))/(exp(+1)); %% 0 for setup1
  rhos(k)=rhos(k)-drhosf*(exp((zv-z0p)/abs(z0p)))/(exp(+1)); %% 0 for setup1
end
%
%--------------------------------------------------------
% Fit southern and northern profiles using linear and
% sin^2 functions. fthy is the shape integral in y
%--------------------------------------------------------
%
ywidth = jetwidth/2;
ydist  = ymax-ywidth;
j0 = 1;
j1 = j0+round(ydist/dy);
j2 = j1+round(ywidth/dy);
j3 = j2;
j4 = j3+round(ywidth/dy);
j5 = Mp;
jmidjet=j2;

yfac1=pi*0.5/(j2-j1);
yfac3=pi*0.5/(j4-j3);

for j=j0:j1
  tmpy(j)=0.;
end
for j=j4:j5
  tmpy(j)=0.;
end
for j=j1+1:j2
  tmpy(j)=sin((j-j1)*yfac1)^2;
end
for j=j2+1:j3
  tmpy(j)=1.;
end
for j=j3+1:j4
  tmpy(j)=cos((j-j3)*yfac3)^2;
end

% Now integrate in y to find fthy.
% tmpy is at rho (u) points, fthy at psi (v) points.
fthy=zeros(Mp,1);
fthy(1)=0.;
for j=2:Mp;
  fthy(j)=fthy(j-1)+0.5*(tmpy(j)+tmpy(j-1));
end
yfac3=fthy(Mp);
fthy=fthy./yfac3;

% fill in rho
for k=1:N
  for j=1:Mp
    fsouth=1.-fthy(j);
    fnorth=fthy(j);
    rho(j,k)=rhos(k)*fsouth+rhon(k)*fnorth;
  end
end

for k=1:N
  for j=1:Mp
    for i=1:Lp
      t(k,j,i)=rho(j,k);
    end
  end
end
rho=t;
%
%--------------------------------------------------------
% Add rho perturbation to zonal fields
%--------------------------------------------------------
%
% Perturbation function:
% find wavelength perturbation 2*pi/Kper near deformation 
% radius wavelength so that there is an integer number 
% of wave numbers in x direction
%
if perturb,

 cff=g/rho0;
 bvf(2:N,:,:)=-cff*(rho(2:N,:,:)-rho(1:N-1,:,:))./ ...
                   (zr(2:N,:,:)-zr(1:N-1,:,:));
 mbvf=sqrt(max(max(max(bvf))));
 Lrad=2*pi*mbvf*abs(2*z0)/f0;
 Lper=2*xmax/floor(2*xmax/Lrad);
 Kper=2*pi/Lper;      % perturbation wavenumber
 Cper=0.01;           % perturbation magnitude

 t=t.*(1+ ...
       Cper*cos(Kper*xr).*exp(-(yr./jetwidth).^2).* ...
       exp(zr./abs(z0)));

end
%
%--------------------------------------------------------
% Compute dynamical fields
%--------------------------------------------------------
%
rho=1000+t;
rhosurf=squeeze(rho(N,:,:));
rho_w=.5*(rho(1:N-1,:,:)+rho(2:N,:,:));
dz_w=zr(2:N,:,:)-zr(1:N-1,:,:);
%
%  Compute pressure field
%
pres=0*rho;
pres(1,:,:)=-zr(1,:,:)*1.e4; %  initialize bottom pressure in Pa
for k=1:N-1;
  pres(k+1,:,:)=pres(k,:,:)-rho_w(k,:,:).*g.*dz_w(k,:,:);
end
%
%  compute SSH and remove area averaged
%
ssh=squeeze(pres(N,:,:))./(rhosurf*g);
avgssh=mean(mean(ssh));
ssh=ssh-avgssh;
zeta=ssh;
%
%  Compute geostrophic baroclinic velocities
%
p_u=rho2u_3d(pres);
p_v=rho2v_3d(pres);
px(:,:,2:Lp-1)=p_u(:,:,2:Lp-1)-p_u(:,:,1:Lp-2);
px(:,:,1)=2.*px(:,:,2)-px(:,:,3);
px(:,:,Lp)=2.*px(:,:,Lp-1)-px(:,:,Lp-2);
py(:,2:Mp-1,:)=p_v(:,2:Mp-1,:)-p_v(:,1:Mp-2,:);
py(:,1,:)=2.*py(:,2,:)-py(:,3,:);
py(:,Mp,:)=2.*py(:,Mp-1,:)-py(:,Mp-2,:);
pn3d=tridim(pn,N); 
pm3d=tridim(pm,N); 
f3d=tridim(f,N); 
u_r=-pn3d.*py./(rho0*f3d);
v_r= pm3d.*px./(rho0*f3d);
u=rho2u_3d(u_r);
v=rho2v_3d(v_r);
%
% Compute barotropic velocities
%
zw=zlevs(h,zeta,theta_s,theta_b,hc,N,'w');
zr=zlevs(h,zeta,theta_s,theta_b,hc,N,'r');
dz=zw(2:end,:,:)-zw(1:end-1,:,:);
dzu=0.5*(dz(:,:,1:end-1)+dz(:,:,2:end));
dzv=0.5*(dz(:,1:end-1,:)+dz(:,2:end,:));
hu(:,:)=sum(dzu.*u);
hv(:,:)=sum(dzv.*v);
D_u(:,:)=sum(dzu);
D_v(:,:)=sum(dzv);
ubar(:,:)=hu./D_u;
vbar(:,:)=hv./D_v;
%
%--------------------------------------------------------
% Create and fill netcdf files
%--------------------------------------------------------
%
%  Create the grid file
%
[Mp,Lp]=size(zeta);
M=Mp-1;
L=Lp-1;
disp(' ')
disp(['GRID DIMENSIONS Lm Mm : ',num2str(L-1),'  ',num2str(M-1)])
nc=netcdf(parent_grd, 'clobber');
redef(nc);
nc('xi_rho') = Lp;
nc('eta_rho') = Mp;
nc('xi_psi') = L;
nc('eta_psi') = M;
nc('one') = 1;
nc{'el'} = ncdouble('one');
nc{'xl'} = ncdouble('one');
nc{'spherical'} = ncchar('one');
nc{'h'} = ncdouble('eta_rho', 'xi_rho');
nc{'f_data'} = ncdouble('eta_rho', 'xi_rho');
nc{'pm'} = ncdouble('eta_rho', 'xi_rho');
nc{'pn'} = ncdouble('eta_rho', 'xi_rho');
nc{'x_rho'} = ncdouble('eta_rho', 'xi_rho');
nc{'y_rho'} = ncdouble('eta_rho', 'xi_rho');
nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho');
endef(nc);
nc.title = ncchar(title_case);
nc.title = title_case;
nc.date = ncchar(date);
nc.date = date;
nc.type = ncchar('ROMS grid file');
nc.type = 'ROMS grid file';
%
%  fill the grid file
%
nc{'xl'}(:)=x(end)-x(1);
nc{'el'}(:)=y(end)-y(1);
nc{'spherical'}(:)='F';
nc{'h'}(:)=H0;

nc{'pm'}(:)=1/dx;
nc{'pn'}(:)=1/dy;

% nc{'f_data'}(:)=f;
% nc{'x_rho'}(:)=X;
% nc{'y_rho'}(:)=Y;
% nc{'mask_rho'}(:)=1+0.*Y;
close(nc);
% 
% %
% %  Create and fill the initial file
% %
% create_inifile(parent_ini,parent_grd,title_case,...
%                theta_s,theta_b,hc,N,0,'clobber')
% nc=netcdf(parent_ini,'write');
% nc{'u'}(:) =  u; 
% nc{'v'}(:) =  v; 
% nc{'zeta'}(:) =  zeta; 
% nc{'ubar'}(:) =  ubar; 
% nc{'vbar'}(:) =  vbar; 
% nc{'temp'}(:) =  t; 
% close(nc)
% %
% %  Create and fill the climatology file
% %
% create_climfile(parent_clm,parent_grd,title_case,...
%                 theta_s,theta_b,hc,N,[25 75],100,'clobber')
% nc=netcdf(parent_clm,'write');
% nc{'u'}(1,:,:,:) =  u;
% nc{'v'}(1,:,:,:) =  v;
% nc{'zeta'}(1,:,:,:) =  zeta;
% nc{'ubar'}(1,:,:,:) =  ubar;
% nc{'vbar'}(1,:,:,:) =  vbar;
% nc{'temp'}(1,:,:,:) =  t;
% nc{'u'}(2,:,:,:) =  u;
% nc{'v'}(2,:,:,:) =  v;
% nc{'zeta'}(2,:,:,:) =  zeta;
% nc{'ubar'}(2,:,:,:) =  ubar;
% nc{'vbar'}(2,:,:,:) =  vbar;
% nc{'temp'}(2,:,:,:) =  t;
% close(nc)
%
%--------------------------------------------------------
% Make a few plots
%--------------------------------------------------------
%
figure;
plot(rhos,zr1d,'r'); hold on; 
plot(rhon,zr1d,'b');
plot(0.5*(rhos+rhon),zr1d,'k'); hold off
title('Density profiles')
%
X=X./1000; Y=Y./1000;
xr=xr./1000; yr=yr./1000;
zu=0.5*(zr(:,:,1:end-1)+zr(:,:,2:end));
xu=0.5*(xr(:,:,1:end-1)+xr(:,:,2:end));
yu=0.5*(yr(:,:,1:end-1)+yr(:,:,2:end));
%
figure
ur=u2rho_2d(squeeze(u(N,:,:)));
vr=v2rho_2d(squeeze(v(N,:,:)));
spd=sqrt(ur.^2+vr.^2);
pcolor(X,Y,spd)
shading flat
axis image
colorbar
hold on
quiver(X,Y,ur,vr,'k')
hold off
title('Velocity')
%
figure
[M,L]=size(Y);
imid=round(L/2);
contourf(squeeze(yr(:,:,imid)),squeeze(zr(:,:,imid)),squeeze(t(:,:,imid)),20);
colorbar
title('Temperature')
%
figure
[M,L]=size(Y);
imid=round(L/2);
contourf(squeeze(yu(:,:,imid)),squeeze(zu(:,:,imid)),squeeze(u(:,:,imid)),10);
colorbar
title('Zonal velocity')
%
figure
contourf(X,Y,zeta,10)
colorbar
axis image
title('Surface height')
%
figure
contourf(rho2v_2d(X),rho2v_2d(Y),squeeze(v(N,:,:)),10)
colorbar
axis image
title('Meridional Velocity (m/s)')
%
us=squeeze(u(N,:,:));
vs=squeeze(v(N,:,:));
nc=netcdf(parent_grd);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
xr=1e-3*nc{'x_rho'}(:);
yr=1e-3*nc{'y_rho'}(:);
[xu,xv,xp]=rho2uvp(xr);
[yu,yv,yp]=rho2uvp(yr);
[fu,fv,fp]=rho2uvp(f);
close(nc);
[vort]=vorticity(us,vs,pm,pn);
vort=vort./fp;
figure
map=colormap(jet(20));
map(10:11,:)=1*[1 1 1 ; 1 1 1];
colormap(map)
contourf(xp,yp,vort,20); shading flat
colorbar
caxis([-0.8 0.8])
axis image
title('Vorticity/f')

