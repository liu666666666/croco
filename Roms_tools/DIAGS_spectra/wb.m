%==============================================================
% COMPUTE KE DIAGNOSTICS OF EMBEDDED GRID SOLUTIONS
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================

clear all
close all

diags_params
snap=1;
if snap,
  lit=[29:33];
else
  lit=[3:96];
end
coastfile = ['/data2/coastline_l.mat'];
print_and_keep=1;

filter_smart = 1;
navg         = 1;      % time filter: avg period =navg*2+1 days
R            = 50.e3; % space filter: decorrelation radius [m]

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
  pm=grd.pm;    pm=pm(imin:imax,jmin:jmax);
  pn=grd.pn;    pn=pn(imin:imax,jmin:jmax);
  lon=grd.lonr; lon=lon(imin:imax,jmin:jmax);
  lat=grd.latr; lat=lat(imin:imax,jmin:jmax);
  xr=grd.xr;    xr=xr(imin:imax,jmin:jmax);
  yr=grd.yr;    yr=yr(imin:imax,jmin:jmax);
  h=h(imin:imax,jmin:jmax);
  zrw=perm(zlevs(perm(h),0,theta_s,theta_b,hc,N,'r'));
  zrw=squeeze(zrw(:,:,kchoice-1:kchoice));

  wbsub=zeros(L,M);
  cffb=-g./rho0;
  irec=0;

  for it=lit(1)+navg:lit(end)-navg

    irec=irec+1;
    disp(['processing rec ' num2str(it) ]);

    % Total T,S,rho
    tmp=rnt_loadvar_partialz(ctlavg,it,'temp',kchoice,kchoice);
    ttot=tmp(imin:imax,jmin:jmax);
    tmp=rnt_loadvar_partialz(ctlavg,it,'salt',kchoice,kchoice);
    stot=tmp(imin:imax,jmin:jmax);
    rtot=rnt_rho_potential_nomex(ttot,stot);

    % Filtered and subfiltered T,S
    if filter_smart,
      tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'temp',kchoice,kchoice)),3);
      tmp=tmp(imin:imax,jmin:jmax);
      tfilt=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
    else
      tmp=low_filter1(ctlavg,'temp',it,kchoice,nhan,navg);
      tfilt=tmp(imin:imax,jmin:jmax);
    end
    tsub=ttot-tfilt;
    %
    if filter_smart,
      tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'salt',kchoice,kchoice)),3);
      tmp=tmp(imin:imax,jmin:jmax);
      sfilt=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
    else
      tmp=low_filter1(ctlavg,'salt',it,kchoice,nhan,navg);
      sfilt=tmp(imin:imax,jmin:jmax);
    end
    ssub=stot-sfilt;

    % Filtered and subfiltered rho
    rfilt=rnt_rho_potential_nomex(tfilt,sfilt);
    rsub=rtot-rfilt;

    % Total, filtered and subfiltered  w
    tmp=rnt_loadvar_partialz(ctlavg,it,'w',kchoice,kchoice);
    tmp(tmp>1e30)=NaN;
    wtot=tmp(imin:imax,jmin:jmax);
    %
    if filter_smart,
      tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'w',kchoice,kchoice)),3);
      tmp=tmp(imin:imax,jmin:jmax);
      wfilt=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
    else
      tmp=low_filter1(ctlavg,'w',it,kchoice,nhan,navg);
      wfilt=tmp(imin:imax,jmin:jmax);
    end
    wsub=wtot-wfilt;

    % Buoyancy flux associated with submesoscale
    % b=-g/rho0*rho
    tmp=cffb*wsub.*rsub;
    wbsub=wbsub+tmp;

  end

  tmp=length([1+navg:length(lit)-navg]);
  wbsub=squeeze(wbsub./tmp);

  %x1=dx.*[0:L-1]/1000;y1=dx.*[0:M-1]/1000;
  %[X,Y]=ndgrid(x1,y1);
  %figure;pcolor(X,Y,wbsub);shading flat;colorbar

  %lonw=min(min(lon)); lone=max(max(lon));
  %lats=min(min(lat)); latn=max(max(lat));
  lonw=228; lone=239.9; lats=-3.8; latn=8.4;

  if snap,
    cmin0=-4.; cint0=0.4; cmax0=4.;
  else
    cmin0=-1.; cint0=0.1; cmax0=1.;
  end
  cff_scale=1.e7;
  wbsub=cff_scale*wbsub;

  lx=20; ly=20;
  figure('units','centimeters','position', ...
           [0 0 lx ly],'paperpositionmode','auto')
  subplot('position',[0 0.15 1 .75]);
  m_proj('mercator','lon',[lonw lone],'lat',[lats latn]);

  m_contourf(lon,lat,wbsub,[cmin0:cint0:cmax0]); shading flat;
  %m_pcolor(lon,lat,wbsub); shading flat;

  caxis([cmin0 cmax0]);
  m_usercoast(coastfile,'patch',[.0 .0 .0]);
  m_grid('box','fancy','fontsize',16);
  title('ROMS submesoscale KE injection [10^{-7} m^3s^{-3}]','Fontsize',15)
  my_colorbar([cmin0 cmax0 cint0],'h','',20);
  set(gcf, 'PaperPositionMode', 'auto');
  set(gca,'fontsize',16)

  if print_and_keep & snap,
    outname=['wb_TIW_l',num2str(ilevel),'_s',num2str(season)];
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




