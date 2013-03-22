%==============================================================
% COMPUTE T FLUX DIAGNOSTICS OF EMBEDDED GRID SOLUTIONS
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================

clear all
%close all

diags_params
snap=0;
lit=[3:93];
coastfile = ['/data2/coastline_l.mat'];
print_and_keep=1;

filter_smart = 1;
navg         = 1;      % time filter: avg period =navg*2+1 days
R            = 50.e3; % space filter: decorrelation radius [m]

loaddata=1;
%=============================================================

for ilevel=1:3;

 if loaddata,

  if ilevel==1,
   load eke_l1
   eke_1=hanning_1d(eke_1,3);
  elseif ilevel==2,
   load eke_l2
   eke_2=hanning_1d(eke_2,10);
  elseif ilevel==3,
   load eke_l3
   eke_3=hanning_1d(eke_3,30);
  end

 else  % compute and save data

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

  imin=lims(1);imax=lims(2);jmin=lims(3);jmax=lims(4);
  ikm=2*min(lims(2)-lims(1),lims(4)-lims(3));
  L=imax-imin+1;M=jmax-jmin+1;
  pm=grd.pm;    pm=pm(imin:imax,jmin:jmax);
  pn=grd.pn;    pn=pn(imin:imax,jmin:jmax);
  lon=grd.lonr; lon=lon(imin:imax,jmin:jmax);
  lat=grd.latr; lat=lat(imin:imax,jmin:jmax);
  xr=grd.xr;    xr=xr(imin:imax,jmin:jmax);
  yr=grd.yr;    yr=yr(imin:imax,jmin:jmax);

  ctlavg=ctlload(root(2:end),avg,model,0,64,32);

  eke=zeros(L,M);

  for it=lit(1)+navg:lit(end)-navg

    disp(['processing rec ' num2str(it) ]);

    % Total u,v
    tmp=rnt_loadvar_partialz(ctlavg,it,'u',kchoice,kchoice);
    tmp=perm(u2rho(perm(tmp)));
    utot=tmp(imin:imax,jmin:jmax);
    %
    tmp=rnt_loadvar_partialz(ctlavg,it,'v',kchoice,kchoice);
    tmp=perm(v2rho(perm(tmp)));
    vtot=tmp(imin:imax,jmin:jmax);
    %
    tmp=rnt_loadvar_partialz(ctlavg,it,'w',kchoice,kchoice);
    wtot=tmp(imin:imax,jmin:jmax);

    % Filtered and subfiltered u,v,w
    if filter_smart,
      tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'u',kchoice,kchoice)),3);
      tmp=perm(u2rho(perm(tmp)));
      tmp=tmp(imin:imax,jmin:jmax);
      ufilt=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
    else
      tmp=low_filter1(ctlavg,'u',it,kchoice,nhan,navg);
      tmp=perm(u2rho(perm(tmp)));
      ufilt=tmp(imin:imax,jmin:jmax);
    end
    usub=utot-ufilt;
    %
    if filter_smart,
      tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'v',kchoice,kchoice)),3);
      tmp=perm(v2rho(perm(tmp)));
      tmp=tmp(imin:imax,jmin:jmax);
      vfilt=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
    else
      tmp=low_filter1(ctlavg,'v',it,kchoice,nhan,navg);
      tmp=perm(v2rho(perm(tmp)));
      vfilt=tmp(imin:imax,jmin:jmax);
    end
    vsub=vtot-vfilt;
    %
    %tmp=low_filter1(ctlavg,'w',it,kchoice,nhan,navg);
    %wfilt=tmp(imin:imax,jmin:jmax);
    %wsub=wtot-wfilt;

    eke=0.5*(usub.^2+vsub.^2); %+wsub.^2);

  end

  tmp=length([1+navg:length(lit)-navg]);
  eke=eke./tmp;
  eke=mean(eke,1);

  if ilevel==1,
   eke_1=eke;
   lat_1=lat(1,:);
   save eke_l1 lat_1 eke_1
  elseif ilevel==2,
   eke_2=eke;
   lat_2=lat(1,:);
   save eke_l2 lat_2 eke_2
  elseif ilevel==3,
   eke_3=eke;
   lat_3=lat(1,:);
   save eke_l3 lat_3 eke_3
  end

 end % load or save data

end

% PLOT EKE

lx=20; ly=20;
figure('units','centimeters','position', ...
           [0 0 lx ly],'paperpositionmode','auto')
pt1=plot(lat_1,eke_1); hold on;
pt2=plot(lat_2,eke_2); hold on;
pt3=plot(lat_3,eke_3); hold on;
set(gca,'fontsize',16)

set(pt1,'Linewidth',2,'Linestyle',':','Color','k')
set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
set(pt3,'Linewidth',2,'Linestyle','-','Color','k')
legend('COARSE 36km','MEDIUM 12km','FINE 4km','Location','NorthWest')

set(gca,'Yscale','linear','Xscale','linear');
set(gca,'Xlim',[-4 10],'Ylim',[0 1.2e-4]);
xlabel('Latitude')
ylabel('Submesoscale EKE [m^2/s^2]')
set(gca,'fontsize',16)

if print_and_keep,
   outname=['eke_prof_s',num2str(season)];
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

