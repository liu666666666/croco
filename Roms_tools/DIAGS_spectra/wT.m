%==============================================================
% COMPUTE T FLUX DIAGNOSTICS OF EMBEDDED GRID SOLUTIONS
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================

clear all
close all

diags_params

coastfile = ['/data2/coastline_l.mat'];
print_and_keep=1;

filter_smart = 0;
navg         = 1;      % time filter: avg period =navg*2+1 days
R            = 50.e3; % space filter: decorrelation radius [m]

snap=1;
if snap,
  lit=[31-navg:31+navg];
else
  lit=[3:96];
end

%=============================================================

for ilevel=1:3;

  nhan=3.^(ilevel-1);
  dx=dx_lev(ilevel);
  model=[model0 num2str(ilevel)];
  ctl=ctlload(root(2:end),his,model,0,lastfile_indx,Nrec);
  lims=squeeze(lims_lev(:,ilevel+1));
  if mod(lims(2)-lims(1)+1,2)==1
     lims(2)=lims(2)+1;
  end
  if mod(lims(4)-lims(3)+1,2)==1
     lims(4)=lims(4)+1;
  end

  grd=rnt_gridload(model);
  theta_s=grd.thetas;theta_b=grd.thetab;hc=grd.hc;Method=grd.Method;N=grd.N;h=grd.h;
  disp(['hc is set to ' num2str(hc) ' remember bug mc2'])

  ctlavg=ctlload(root(2:end),avg,model,0,64,32);

  imin=lims(1);imax=lims(2);jmin=lims(3);jmax=lims(4);
  ikm=2*min(lims(2)-lims(1),lims(4)-lims(3));
  L=imax-imin+1;M=jmax-jmin+1;
  kmin=kchoice-1;kmax=kchoice+1;
  pm=grd.pm;    pm=pm(imin:imax,jmin:jmax);
  pn=grd.pn;    pn=pn(imin:imax,jmin:jmax);
  lon=grd.lonr; lon=lon(imin:imax,jmin:jmax);
  lat=grd.latr; lat=lat(imin:imax,jmin:jmax);
  h=h(imin:imax,jmin:jmax);
  zw=perm(zlevs(perm(h),0,theta_s,theta_b,hc,N,'w'));
  dz=zeros(size(zw));
  dz(:,:,2:N)=zw(:,:,3:N+1)-zw(:,:,2:N);
  dz(:,:,1)=dz(:,:,2);
  dz=dz(:,:,kchoice:N);

  wt=zeros(L,M);
  hbl=zeros(L,M);
  cffb=-g./rho0;

  for it=lit(1)+navg:lit(end)-navg

    disp(['processing rec ' num2str(it) ]);

    % Total T
    tmp=rnt_loadvar_partialz(ctlavg,it,'temp',kchoice,kchoice);
    ttot=tmp(imin:imax,jmin:jmax);
    % Filtered and subfiltered T
    if filter_smart,
      tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'temp',kchoice,kchoice)),3);
      tmp=tmp(imin:imax,jmin:jmax);
      tfilt=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
    else
      tmp=low_filter1(ctlavg,'temp',it,kchoice,nhan,navg);
      tfilt=tmp(imin:imax,jmin:jmax);
    end
    tsub=ttot-tfilt;

    % Total, filtered and subfiltered  w
    tmp=rnt_loadvar_partialz(ctlavg,it,'w',kchoice,kchoice);
    tmp(tmp>1e30)=NaN;
    wtot=tmp(imin:imax,jmin:jmax);
    if filter_smart,
      tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'w',kchoice,kchoice)),3);
      tmp=tmp(imin:imax,jmin:jmax);
      wfilt=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
    else
      tmp=low_filter1(ctlavg,'w',it,kchoice,nhan,navg);
      wfilt=tmp(imin:imax,jmin:jmax);
    end
    wsub=wtot-wfilt;

    % T flux associated with submesoscale
    tmp=wsub.*tsub;
    wt=wt+tmp;

    % HBL
    tmp=rnt_loadvar(ctlavg,it,'hbl');
    hbl=hbl+tmp(imin:imax,jmin:jmax);

  end

  tmp=length([1+navg:length(lit)-navg]);
  wt=wt./tmp;
  hbl=hbl./tmp;
  dwtdz=wt./sum(dz,3);

  if ilevel==1,
   dwtdz_1=mean(dwtdz,1);
   hbl_1=mean(hbl,1);
   lat_1=lat(1,:);
  elseif ilevel==2,
   dwtdz_2=mean(dwtdz,1);
   hbl_2=mean(hbl,1);
   lat_2=lat(1,:);
  elseif ilevel==3,
   dwtdz_3=mean(dwtdz,1);
   hbl_3=mean(hbl,1);
   lat_3=lat(1,:);
  end

  %======== PLOT =============================

  %x1=dx.*[0:L-1]/1000;y1=dx.*[0:M-1]/1000;
  %[X,Y]=ndgrid(x1,y1);
  %figure;pcolor(X,Y,wbsub);shading flat;colorbar

  %lonw=min(min(lon)); lone=max(max(lon));
  %lats=min(min(lat)); latn=max(max(lat));
  lonw=228; lone=239.9; lats=-3.8; latn=8.4;

  cff_scale=86400;
  if snap,
    cmin=-2.; cint=0.1; cmax=2.;
  else
    cmin=-10.; cint=1; cmax=10.;
  end

  lx=20; ly=20;
  figure('units','centimeters','position', ...
           [0 0 lx ly],'paperpositionmode','auto')
 
  map=colormap(cool(40));
%  map(10:11,:)=0.95*[1 0.9 0.9 ; 1 0.9 0.9];
  map(20:21,:)=1*[1 1 1 ; 1 1 1];
  colormap(map)

  subplot('position',[0 0.15 1 .75]);
  m_proj('mercator','lon',[lonw lone],'lat',[lats latn]);

  m_contourf(lon,lat,cff_scale*wt,[cmin:cint:cmax],'k:'); shading flat;
  %m_pcolor(lon,lat,dwtdz); shading flat;

  caxis([cmin cmax]);
  m_usercoast(coastfile,'patch',[.0 .0 .0]);
  m_grid('box','fancy','fontsize',16);
  title('ROMS submesoscale T flux [m K day^{-1}]','Fontsize',15)
  my_colorbar([cmin cmax cint],'h','',20);
  set(gcf, 'PaperPositionMode', 'auto');
  set(gca,'fontsize',16)

  if print_and_keep,
    if snap,
     outname=['wt_TIW_l',num2str(ilevel),'_s',num2str(season)];
    else
     outname=['mwt_TIW_l',num2str(ilevel),'_s',num2str(season)];
    end
    warning off
    eval(['print -painter -depsc2 ',outname,'.eps;'])
    warning on
    eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])

    %eval(['! mv ',outname,'.eps ',dirout_EPS]);
    %eval(['! mv ',outname,'.jpg ',dirout_JPG]);
    eval(['! scp ',outname,'.eps ',dirout_EPS]);
    eval(['! scp ',outname,'.jpg ',dirout_JPG]);
    eval(['! rm -f ',outname,'.eps']);
    eval(['! rm -f ',outname,'.jpg']);
  end

end

return

% WT

lx=20; ly=20;
figure('units','centimeters','position', ...
           [0 0 lx ly],'paperpositionmode','auto')
pt1=plot(cff_scale*dwtdz_1,lat_1); hold on;
pt2=plot(cff_scale*dwtdz_2,lat_2); hold on;
pt3=plot(cff_scale*dwtdz_3,lat_3); hold on;
set(gca,'fontsize',16)

hl=line([0 0],[-4 10]); hold off
set(hl,'Linewidth',2,'Linestyle','-','Color','r')

set(pt1,'Linewidth',2,'Linestyle',':','Color','k')
set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
set(pt3,'Linewidth',2,'Linestyle','-','Color','k')
legend('COARSE 36km','MEDIUM 12km','FINE 4km','Location','NorthEast')

set(gca,'Yscale','linear','Xscale','linear');
set(gca,'Ylim',[-4 10],'Xlim',[-2 2]);
ylabel('Latitude')
xlabel('Submesoscale vertical transport [K/day]')
set(gca,'fontsize',16)

% HBL

lx=20; ly=20;
figure('units','centimeters','position', ...
           [0 0 lx ly],'paperpositionmode','auto')
hbl_1(hbl_1==0)=NaN;hbl_2(hbl_2==0)=NaN;hbl_3(hbl_3==0)=NaN;
pt1=plot(hbl_1,lat_1); hold on;
pt2=plot(hbl_2,lat_2); hold on;
pt3=plot(hbl_3,lat_3); hold on;
set(gca,'fontsize',16)

hl=line([0 0],[-4 10]); hold off
set(hl,'Linewidth',2,'Linestyle','-','Color','r')

set(pt1,'Linewidth',2,'Linestyle',':','Color','k')
set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
set(pt3,'Linewidth',2,'Linestyle','-','Color','k')
legend('COARSE 36km','MEDIUM 12km','FINE 4km','Location','NorthEast')

set(gca,'Yscale','linear','Xscale','linear');
set(gca,'Ylim',[-4 10],'Xlim',[0 50]);
ylabel('Latitude')
xlabel('Hbl [m]')
set(gca,'fontsize',16)


