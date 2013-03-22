clear all;
close all;
%=================================================================
% Barotropic Vorticity Budget
%
% Compute all terms of the equation for the rate of change of 
% the vorticity of the depth-integrated flow
%
% Procedure: depth integrate the mean momentum equations computed
% in ROMS and then cross differentiate. 
%
% The routine reads time-averaged momentum budget terms as well as 
% time-averaged velocities and free surface.
%
% See Wells (1995, JPO, 25, 2569-2582) for a reference and 
% Mertz and Wright (1992, JPO, 22, 301-305) and Cane et al.
% (1998, JPO 28, 519-526) for a discussion
%
% Patrick Marchesiello - IRD 2007
%===================================================================

DIR       = 'ROMS_FILES/';

diafile   = [DIR,'roms_diaM_avg.nc'];
mfile     = [DIR,'roms_avg.nc'];

grdfile   = [DIR,'roms_grd.nc'];
clmfile   = [DIR,'roms_clm.nc'];
frcfile   = [DIR,'roms_frc.nc'];
coastfile = [DIR,'coastline_l.mat'];

ISEA    =  5;   % 5 for annual mean

%
%lonw=-85; lone=-75;   lats=-15; latn=-7;     % PERU TEST
lonw=8; lone=22;   lats=-38; latn=-25.5;      % BENGUELA_LR TEST

compute    = 1;      % 1 : compute vorticity terms and store in *.mat file
                     % 0 : read vorticity terms in *.mat file 

local_budget = 0;
zonal_budget = 0;
maps         = 0;    % 1 : make a few plots

plot_scale   = 5.e7;
subsample    = 0;
plot_res     = 1/4;
%====================================================================
%
% Get grid parameters
%
nc=netcdf(grdfile);
mask=nc{'mask_psi'}(:);
mask_u=nc{'mask_u'}(:);
mask_v=nc{'mask_v'}(:);
mask_r=nc{'mask_rho'}(:);
lat_psi=nc{'lat_psi'}(:);
lon_psi=nc{'lon_psi'}(:);
lat_v=nc{'lat_v'}(:);
lon_v=nc{'lon_v'}(:);
lat_r=nc{'lat_rho'}(:);
lon_r=nc{'lon_rho'}(:);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
h=nc{'h'}(:);
f=nc{'f'}(:);
close(nc);
[Mp Lp]=size(pm);
L=Lp-1; M=Mp-1;
mask(find(mask==0))=NaN;
mask_u(find(mask_u==0))=NaN;
mask_v(find(mask_v==0))=NaN;
mask_r(find(mask_r==0))=NaN;
h_psi=0.25*(h(1:end-1,1:end-1)+h(1:end-1,2:end)+ ...
            h(2:end,1:end-1)+h(2:end,2:end));
%
% Frame computational domain inside of a 
% smaller domain (roms domain minus ZONE [degrees])
%
ZONE=1.0; 
lon1=min(min(lon_psi)); lon2=max(max(lon_psi));
lat1=min(min(lat_psi)); lat2=max(max(lat_psi));
dlon=lon_psi(1,2)-lon_psi(1,1);
dlat=lat_psi(2,1)-lat_psi(1,1);
lonw=max(lonw,lon1+ZONE);
lone=min(lone,lon2-ZONE);
lats=max(lats,lat1+ZONE);
latn=min(latn,lat2-ZONE);
Istr=min(find(lon_psi(1,:)>lonw-2*dlon & lon_psi(1,:)<lonw+2*dlon));
Iend=max(find(lon_psi(1,:)>lone-2*dlon & lon_psi(1,:)<lone+2*dlon));
Jstr=min(find(lat_psi(:,1)>lats-2*dlat & lat_psi(:,1)<lats+2*dlat));
Jend=max(find(lat_psi(:,1)>latn-2*dlat & lat_psi(:,1)<latn+2*dlat));

if compute,
%===========================================================================
%
% Starts computation of vorticity balance terms and store them in
% *.mat file (compute=1); or read them in *.mat file from a 
% previous computation (compute=0).
%
%===========================================================================
%
% get mean dynamical variables
%
nc=netcdf(mfile);
zeta=squeeze(nc{'zeta'}(ISEA,:,:));
ubar=squeeze(nc{'ubar'}(ISEA,:,:));
vbar=squeeze(nc{'vbar'}(ISEA,:,:));
u=squeeze(nc{'u'}(ISEA,:,:));
v=squeeze(nc{'v'}(ISEA,:,:));
close(nc);

%
% Vertical coordinates
%
nc=netcdf(clmfile);
N       = length(nc('s_rho'));
theta_s = nc{'theta_s'}(:);
theta_b = nc{'theta_b'}(:);
hc      = nc{'hc'}(:);
close(nc)
zw= zlevs(h,zeta,theta_s,theta_b,hc,N,1,'w');
Hz(1:N,:,:)=zw(2:N+1,:,:)-zw(1:N,:,:);
zu=rho2u_3d(zw); zv=rho2v_3d(zw);
%clear zw

%
% Get wind stress
% for later wind stress curl computation
% (unit is already  m^2/s^2 
%
nc=netcdf(frcfile);
rho0=1025;
u_taus=squeeze(mean(nc{'sustr'}(:)))./rho0;
v_taus=squeeze(mean(nc{'svstr'}(:)))./rho0;
close(nc);

%
% Get transport by the mean flow
%

[Mu_xadv,Mu_yadv,Mv_xadv,Mv_yadv]=roms_hadv_uv(u,v,Hz,pm,pn);

%
% Get U and V Momentum balance terms
%
nc=netcdf(diafile);

u_xadv   = squeeze(nc{'u_xadv'}(ISEA,:,:,:));
u_yadv   = squeeze(nc{'u_yadv'}(ISEA,:,:,:));
u_cor    = squeeze(nc{'u_cor'}(ISEA,:,:,:));
u_vmix   = squeeze(nc{'u_vmix'}(ISEA,:,:,:));
u_hmix   = squeeze(nc{'u_hmix'}(ISEA,:,:,:));
u_Prsgrd = squeeze(nc{'u_Prsgrd'}(ISEA,:,:,:));
u_rate   = squeeze(nc{'u_rate'}(ISEA,:,:,:));

v_xadv   = squeeze(nc{'v_xadv'}(ISEA,:,:,:));
v_yadv   = squeeze(nc{'v_yadv'}(ISEA,:,:,:));
v_cor    = squeeze(nc{'v_cor'}(ISEA,:,:,:));
v_vmix   = squeeze(nc{'v_vmix'}(ISEA,:,:,:));
v_hmix   = squeeze(nc{'v_hmix'}(ISEA,:,:,:));
v_Prsgrd = squeeze(nc{'v_Prsgrd'}(ISEA,:,:,:));
v_rate   = squeeze(nc{'v_rate'}(ISEA,:,:,:));

close(nc);

debug=0;
if debug,
%  utot=-u_rate+u_xadv+u_yadv+u_zadv+u_cor+u_vmix+u_Prsgrd;
  figure; pcolor(squeeze(u_cor(N,:,:)));
  shading flat; colorbar;
  stop

%  vtot=-v_rate+v_xadv+v_yadv+v_zadv+v_cor+v_vmix+v_Prsgrd;
  figure; pcolor(squeeze(vtot(N,:,:)*1.e8).*mask_v);
  shading flat;caxis([-1 1]); colorbar;
  stop
end

%
% Get depth integrated Momentum terms from ROMS diag file
%
dzu=zu(2:N+1,:,:)-zu(1:N,:,:);
dzv=zv(2:N+1,:,:)-zv(1:N,:,:); 
%dzu=rho2u_3d(Hz); dzv=rho2v_3d(Hz);

u_xadv_zint   = squeeze(sum(u_xadv.*dzu,1))   .*mask_u;
u_yadv_zint   = squeeze(sum(u_yadv.*dzu,1))   .*mask_u;
Mu_xadv_zint  = squeeze(sum(Mu_xadv.*dzu,1))  .*mask_u;
Mu_yadv_zint  = squeeze(sum(Mu_yadv.*dzu,1))  .*mask_u;
u_cor_zint    = squeeze(sum(u_cor.*dzu,1))    .*mask_u;
u_vmix_zint   = squeeze(sum(u_vmix.*dzu,1))   .*mask_u;
u_hmix_zint   = squeeze(sum(u_hmix.*dzu,1))   .*mask_u;
u_Prsgrd_zint = squeeze(sum(u_Prsgrd.*dzu,1)) .*mask_u;
u_rate_zint   = squeeze(sum(u_rate.*dzu,1))   .*mask_u;

v_xadv_zint   = squeeze(sum(v_xadv.*dzv,1))   .*mask_v;
v_yadv_zint   = squeeze(sum(v_yadv.*dzv,1))   .*mask_v;
Mv_xadv_zint  = squeeze(sum(Mv_xadv.*dzv,1))  .*mask_v;
Mv_yadv_zint  = squeeze(sum(Mv_yadv.*dzv,1))  .*mask_v;
v_cor_zint    = squeeze(sum(v_cor.*dzv,1))    .*mask_v;
v_vmix_zint   = squeeze(sum(v_vmix.*dzv,1))   .*mask_v;
v_hmix_zint   = squeeze(sum(v_hmix.*dzv,1))   .*mask_v;
v_Prsgrd_zint = squeeze(sum(v_Prsgrd.*dzv,1)) .*mask_v;
v_rate_zint   = squeeze(sum(v_rate.*dzv,1))   .*mask_v;

u_taus_zint   = u_taus  .*mask_u;
v_taus_zint   = v_taus  .*mask_v;

%    compute avg terms for STRETCHING and JEBAR
hu=-squeeze(zu(1,:,:));
hv=-squeeze(zv(1,:,:));
u_cor_zavg    = squeeze(sum(u_cor.*dzu,1))    ./hu.*mask_u;
v_cor_zavg    = squeeze(sum(v_cor.*dzv,1))    ./hv.*mask_v;
u_Prsgrd_zavg = squeeze(sum(u_Prsgrd.*dzu,1)) ./hu.*mask_u;
v_Prsgrd_zavg = squeeze(sum(v_Prsgrd.*dzv,1)) ./hv.*mask_v;

clear dzu dzv 

%-----------------------------------------------------
% Compute vorticity terms:
% Cross differentiate Momemtum U and V equations
%-----------------------------------------------------
odx1 = .5*(pm(2:end,:)+pm(1:end-1,:));
odx  = .5*(odx1(:,1:end-1)+odx1(:,2:end));
ody1 = .5*(pn(2:end,:)+pn(1:end-1,:));
ody  = .5*(ody1(:,1:end-1)+ody1(:,2:end));

% Vadv
u_xterm_y=(u_xadv_zint(2:end,:)-u_xadv_zint(1:end-1,:)).*ody;
u_yterm_y=(u_yadv_zint(2:end,:)-u_yadv_zint(1:end-1,:)).*ody;
u_term_y=u_xterm_y+u_yterm_y;

v_xterm_x=(v_xadv_zint(:,2:end)-v_xadv_zint(:,1:end-1)).*odx;
v_yterm_x=(v_yadv_zint(:,2:end)-v_yadv_zint(:,1:end-1)).*odx;
v_term_x=v_xterm_x+v_yterm_x;

Vadv=(v_term_x-u_term_y).*mask;

% MVadv
u_xterm_y=(Mu_xadv_zint(2:end,:)-Mu_xadv_zint(1:end-1,:)).*ody;
u_yterm_y=(Mu_yadv_zint(2:end,:)-Mu_yadv_zint(1:end-1,:)).*ody;
u_term_y=u_xterm_y+u_yterm_y;

v_xterm_x=(Mv_xadv_zint(:,2:end)-Mv_xadv_zint(:,1:end-1)).*odx;
v_yterm_x=(Mv_yadv_zint(:,2:end)-Mv_yadv_zint(:,1:end-1)).*odx;
v_term_x=v_xterm_x+v_yterm_x;

MVadv=(v_term_x-u_term_y).*mask;

% Vcor
u_term_y=(u_cor_zint(2:end,:)-u_cor_zint(1:end-1,:)).*ody;
v_term_x=(v_cor_zint(:,2:end)-v_cor_zint(:,1:end-1)).*odx;

Vcor=(v_term_x-u_term_y).*mask;

% Vpg
u_term_y=(u_Prsgrd_zint(2:end,:)-u_Prsgrd_zint(1:end-1,:)).*ody;
v_term_x=(v_Prsgrd_zint(:,2:end)-v_Prsgrd_zint(:,1:end-1)).*odx;

Vpg=(v_term_x-u_term_y).*mask;

% Vbpg: Barotropic Pressure term
dedx=0.5*(zeta(2:end,2:end)+zeta(1:end-1,2:end)- ...
          zeta(2:end,1:end-1)-zeta(1:end-1,1:end-1)).*odx;
dhdx=0.5*(h(2:end,2:end)+h(1:end-1,2:end)- ...
          h(2:end,1:end-1)+h(1:end-1,1:end-1)).*odx;
dedy=0.5*(zeta(2:end,2:end)-zeta(2:end,1:end-1)- ...
          zeta(1:end-1,2:end)+zeta(1:end-1,1:end-1)).*ody;
dhdy=0.5*(h(2:end,2:end)+h(2:end,1:end-1)- ...
          h(1:end-1,2:end)-h(1:end-1,1:end-1)).*ody;

Vbpg=-9.8*(dedx.*dhdy+dedy.*dhdx).*mask;

% Vjebar
u_term_y=(u_Prsgrd_zavg(2:end,:)-u_Prsgrd_zavg(1:end-1,:)).*ody;
v_term_x=(v_Prsgrd_zavg(:,2:end)-v_Prsgrd_zavg(:,1:end-1)).*odx;

Vjebar=h_psi.*(v_term_x-u_term_y).*mask;

% Vstr
u_term_y=(u_cor_zavg(2:end,:)-u_cor_zavg(1:end-1,:)).*ody;
v_term_x=(v_cor_zavg(:,2:end)-v_cor_zavg(:,1:end-1)).*odx;

Vstr=(h_psi.*(v_term_x-u_term_y)-Vcor).*mask;

% Vtau
u_term_y=(u_vmix_zint(2:end,:)-u_vmix_zint(1:end-1,:)).*ody;
v_term_x=(v_vmix_zint(:,2:end)-v_vmix_zint(:,1:end-1)).*odx;

Vtau=(v_term_x-u_term_y).*mask;

% Vtaus
u_term_y=(u_taus_zint(2:end,:)-u_taus_zint(1:end-1,:)).*ody;
v_term_x=(v_taus_zint(:,2:end)-v_taus_zint(:,1:end-1)).*odx;

Vtaus=(v_term_x-u_term_y).*mask;

% Lateral Friction
u_term_y=(u_hmix_zint(2:end,:)-u_hmix_zint(1:end-1,:)).*ody;
v_term_x=(v_hmix_zint(:,2:end)-v_hmix_zint(:,1:end-1)).*odx;

Vfric=(v_term_x-u_term_y).*mask;

% Vrate
u_term_y=(u_rate_zint(2:end,:)-u_rate_zint(1:end-1,:)).*ody;
v_term_x=(v_rate_zint(:,2:end)-v_rate_zint(:,1:end-1)).*odx;

Vrate=(v_term_x-u_term_y).*mask;

% Vtot
Vtot=(-Vrate+Vcor+Vtau+Vpg+Vadv+Vfric).*mask;

save 'vorticity_budget.mat' Vtot Vrate Vadv MVadv Vcor Vpg Vbpg Vjebar Vstr Vtau Vtaus Vfric


else % compute

  load vorticity_budget.mat

end; % compute

%=====================================================================================
%
% COMPUTE MEAN BUDGETS
%
%=====================================================================================

%--------------------------------
% Area mean budget
%--------------------------------

ds=1./(0.0625*(pm(1:end-1,1:end-1)+pm(2:end,1:end-1)+ ...
               pm(1:end-1,2:end)+pm(2:end,2:end)).* ...
              (pn(1:end-1,1:end-1)+pn(2:end,1:end-1)+ ...
               pn(1:end-1,2:end)+pn(2:end,2:end)) );
DS=sum(sum( ds(Jstr:Jend,Istr:Iend) ));

VTMP=Vrate; VTMP(isnan(VTMP))=0;
VRATE=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

VTMP=Vtau; VTMP(isnan(VTMP))=0;
VTAU=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

VTMP=Vtaus; VTMP(isnan(VTMP))=0;
VTAUS=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

VTAUB=VTAU-VTAUS;

VTMP=Vcor; VTMP(isnan(VTMP))=0;
VCOR=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

VTMP=Vadv; VTMP(isnan(VTMP))=0;
VADV=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

VTMP=MVadv; VTMP(isnan(VTMP))=0;
MVADV=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

EVADV=VADV-MVADV;

VTMP=Vfric; VTMP(isnan(VTMP))=0;
VFRIC=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

VTMP=Vpg; VTMP(isnan(VTMP))=0;
VBPT=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

VTMP=Vbpg; VTMP(isnan(VTMP))=0;
VBBPT=sum(sum( VTMP(Jstr:Jend,Istr:Iend).*ds(Jstr:Jend,Istr:Iend) ))./DS;

scale=1.e12;
VRATE=scale*VRATE;
VCOR=scale*VCOR;
VTAU=scale*VTAU;
VTAUS=scale*VTAUS;
VTAUB=scale*VTAUB;
VBPT=scale*VBPT;
VBBPT=scale*VBBPT;
VADV=scale*VADV;
MVADV=scale*MVADV;
EVADV=scale*EVADV;
VFRIC=scale*VFRIC;

VTOT=-VRATE+VTAU+VCOR+VBPT+VADV+VFRIC;

disp(' ')
disp([' ============================================================== '])
disp([' =========     INTEGRATED VORTICITY BUDGET (m/s2)    ========= '])
disp([' ============================================================== '])
disp(' ')
disp([' TOTAL                         = ',num2str(VTOT)]);
disp([' RATE OF CHANGE                = ',num2str(VRATE)]);
disp(' ........')
disp(' Topographic Sverdrup Relation ')
disp(' ............................. ')
disp([' BETA                          = ',num2str(VCOR)]);
disp([' WIND STRESS                   = ',num2str(VTAUS)]);
disp([' BOTTOM FRICTION               = ',num2str(VTAUB)]);
disp([' BOTTOM PRESSURE TORQUE        = ',num2str(VBPT)]);
disp([' ADVECTION                     = ',num2str(VADV)]);
disp([' LATERAL FRICTION              = ',num2str(VFRIC)]);
disp(' ........')
disp([' MEAN ADVECTION                = ',num2str(MVADV)]);
disp([' EDDY ADVECTION                = ',num2str(EVADV)]);
disp(' ........')
disp([' BAROTROPIC BPT                = ',num2str(VBBPT)]);
disp([' BAROCLINIC BPT                = ',num2str(VBPT-VBBPT)]);


if zonal_budget,
%------------------------------------
% Meridional mean  Budget
%------------------------------------

dy=4./(pn(1:end-1,1:end-1)+pn(2:end,1:end-1)+ ...
       pn(1:end-1,2:end)+pn(2:end,2:end));
DY=sum(dy(Jstr:Jend,:));
x=squeeze(lon_psi(1,:));

VTMP=Vrate; VTMP(isnan(VTMP))=0;
ZVRATE=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=Vtau; VTMP(isnan(VTMP))=0;
ZVTAU=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=Vtaus; VTMP(isnan(VTMP))=0;
ZVTAUS=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=Vcor; VTMP(isnan(VTMP))=0;
ZVCOR=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=Vpg; VTMP(isnan(VTMP))=0;
ZVPG=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=Vbpg; VTMP(isnan(VTMP))=0;
ZVBPG=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=Vadv; VTMP(isnan(VTMP))=0;
ZVADV=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=MVadv; VTMP(isnan(VTMP))=0;
ZMVADV=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;

VTMP=Vfric; VTMP(isnan(VTMP))=0;
ZVFRIC=sum( VTMP(Jstr:Jend,:).*dy(Jstr:Jend,:) )./DY;


scale=1e11;
fscale=0.8;
ZVCOR=scale*filter_1D(x,ZVCOR,fscale);
ZVPG=scale*filter_1D(x,ZVPG,fscale);
ZVBPG=scale*filter_1D(x,ZVBPG,fscale);
ZVADV=scale*filter_1D(x,ZVADV,fscale);
ZVMADV=scale*filter_1D(x,ZMVADV,fscale);
ZVTAU=scale*filter_1D(x,ZVTAU,fscale);
ZVTMP=ZVTAU+ZVCOR+ZVPG+ZVADV;
ZVEADV=ZVADV-ZVMADV;
ZVTAUB=ZVTAU-ZVTAUS;
ZVFRIC=scale*filter_1D(x,ZVFRIC,fscale);

figure
hh=plot(x,ZVTMP);
set(hh,'Color','k','LineWidth',2);
legend('TOTAL');
axis([155 188 -100 100])
title(['Mean meridional vorticity terms - Scale factor = 1.10^{11}'])
grid

figure
hh=plot(x,ZVCOR,'-',x,ZVPG+ZVADV,'--');
set(hh,'Color','k','LineWidth',2);
legend('BETA','BPT+ADVECTION');
axis([155 188 -100 100])
title(['Mean meridional vorticity terms - Scale factor = 1.10^{11}'])
grid

figure
hh=plot(x,ZVADV,'-',x,ZVMADV,'--',x,ZVEADV,':');
set(hh,'Color','k','LineWidth',2);
legend('Total ADVECTION','Mean ADVECTION','Eddy ADVECTION');
axis([155 188 -100 100])
title(['Mean meridional vorticity terms - Scale factor = 1.10^{11}'])
grid

figure
%hh=plot(x,ZVCOR,'-',x,ZVPG,'--',x,ZVADV,':',x,ZVFRIC,'--');
hh=plot(x,ZVCOR,'-',x,ZVPG,'--',x,ZVADV,':');
set(hh,'Color','k','LineWidth',2);
%legend('BETA','BPT','ADVECTION','FRICTION');
legend('BETA','BPT','ADVECTION');
axis([155 188 -100 100])
xlabel('Longitude','FontSize',11);
title(['Mean meridional vorticity terms - Scale factor = 1.10^{11}'])
grid

figure
hh=plot(x,ZVCOR,':',x,0.1*ZVBPG,'-',x,0.1*(ZVPG-ZVBPG),'--');
set(hh,'Color','k','LineWidth',2);
legend('BETA','Barotropic BPT /10','Baroclinic BPT /10');
axis([155 188 -100 100])
xlabel('Longitude','FontSize',11);
title(['Mean meridional vorticity terms - Scale factor = 1.10^{11}'])
grid

%figure
%hh=plot(x,ZVCOR,':',x,abs(ZVBPG)-abs(ZVPG-ZVBPG),'-');
%set(hh,'Color','k','LineWidth',2);
%legend('BETA','abs diff in BPT');
%axis([155 188 -100 100])
%xlabel('Longitude','FontSize',11);
%title(['Mean meridional vorticity terms - Scale factor = 1.10^{11}'])
%grid

end % Meridonal mean budget

if local_budget,
%----------------------------
% Local budget
%----------------------------
lon0=170; lat0=-22.5;
[I,J]=find(lon_psi<lon0+1/12 & lon_psi>lon0-1/12 & lat_psi<lat0+1/12 & lat_psi>lat0-1/12);
I=min(I); J=min(J);
disp(' ')
disp([' ============================================================== '])
disp([' =========       LOCAL VORTICITY BUDGET (m3/s2)       ========= '])
disp([' ============================================================== '])
disp(' ')
disp([' TOTAL                         = ',num2str(Vtot(I,J))]);
disp([' RATE OF CHANGE                = ',num2str(Vrate(I,J))]);
disp(' ........')
disp(' Topographic Sverdrup Relation ')
disp(' ............................. ')
disp([' BETA                          = ',num2str(Vcor(I,J))]);
disp([' WIND & BOTTOM STRESS          = ',num2str(Vtau(I,J))]);
disp([' SLOPE CONTROL                 = ',num2str(Vpg(I,J))]);
disp([' ADVECTION                     = ',num2str(Vadv(I,J))]);
disp(' ........')
disp([' MEAN ADVECTION                = ',num2str(MVadv(I,J))]);
disp([' EDDY ADVECTION                = ',num2str(Vadv(I,J)-MVadv(I,J))]);
end

%
%=============================================================================
%
% Make a few plots
%
%=============================================================================
if (maps)

scale=plot_scale; scale2=scale/10;
SCALE=num2str(1/scale); SCALE2=num2str(1/scale2);
foh=1.e7*f./h.*mask_r; axis_foh=[-.2:.02:0];
%foh=h.*mask_r; axis_foh=[0:500:10000];

if subsample,
  [LON,LAT]=meshgrid([lonw:plot_res:lone],[lats:plot_res:latn]);
  Vcor=interp2(lon_psi,lat_psi,Vcor,LON,LAT);
  Vtau=interp2(lon_psi,lat_psi,Vtau,LON,LAT);
  Vtaus=interp2(lon_psi,lat_psi,Vtaus,LON,LAT);
  Vpg=interp2(lon_psi,lat_psi,Vpg,LON,LAT);
  Vbpg=interp2(lon_psi,lat_psi,Vbpg,LON,LAT);
  Vadv=interp2(lon_psi,lat_psi,Vadv,LON,LAT);
  MVadv=interp2(lon_psi,lat_psi,MVadv,LON,LAT);
  Vfric=interp2(lon_psi,lat_psi,Vfric,LON,LAT);
  lon_psi=LON; lat_psi=LAT;
end;

%
% Main Budget
%
minval=0.1;
m_proj('mercator','lon',[lonw lone],'lat',[lats latn]);

figure;
subplot(2,2,1)
var=Vcor.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on; 
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['BETA - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,2)
var=Vtau.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Wind and Bottom Stress - ',SCALE,' m/s^2'])
%title(['Lateral Friction - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,3)
var=Vpg.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Bottom Pressure Torque - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,4)
var=Vadv.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Advection - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

return

% Advection decompositions
%
figure
subplot(2,2,1)
var=Vadv.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['TOTAL ADVECTION - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,2)
var=MVadv.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['MEAN ADVECTION - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,3)
var=(Vadv-MVadv).*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['EDDY ADVECTION - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

% Wind and Bottom stress curl
%
figure
subplot(2,2,1)
var=Vtau.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Wind and Bottom Stress - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,2)
var=Vtaus.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Wind Stress - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,3)
var=(Vtau-Vtaus).*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Bottom Stress - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

% BPT decomposition
%
figure
subplot(2,2,1)
var=Vpg.*scale; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Bottom Pressure Torque - ',SCALE,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,2)
var=Vbpg.*scale2; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Barotropic Bottom Pressure Torque - ',SCALE2,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

subplot(2,2,3)
var=(Vpg-Vbpg).*scale2; var(abs(var)<minval)=NaN;
m_pcolor(lon_psi,lat_psi,var); shading flat
caxis([-1 1]); axis([lonw lone lats latn]); hold on;
[C1,h1]=m_contour(lon_r,lat_r,foh,axis_foh,'k'); set(h1,'LineStyle','--'); hold off;
title(['Baroclinic Bottom Pressure Torque - ',SCALE2,' m/s^2'])
m_usercoast(coastfile,'patch',[.0 .0 .0]);
m_grid('box','fancy');

end  % plot
return

