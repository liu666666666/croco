%==============================================================
% COMPUTE T FLUX DIAGNOSTICS OF EMBEDDED GRID SOLUTIONS
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================

clear all
close all

diags_params
snap=0;
if snap,
  lit=[48:52];
else
  lit=[3:93];
end
coastfile = ['/data2/coastline_l.mat'];
print_and_keep=1;
navg=2;
%=============================================================

for ilevel=1:3;

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

  hbl=zeros(L,M);
  sst=zeros(L,M);

  for it=1+navg:length(lit)-navg

    disp(['processing rec ' num2str(it) ]);

    % HBL
    tmp=rnt_loadvar(ctl,it,'hbl');
    hbl=hbl+tmp(imin:imax,jmin:jmax);

    % SST
    tmp=squeeze(rnt_loadvar_partialz(ctl,it,'temp',N,N));
    sst=sst+tmp(imin:imax,jmin:jmax);

  end

  tmp=length([1+navg:length(lit)-navg]);
  hbl=hbl./tmp;
  sst=sst./tmp;

  if ilevel==1,
   hbl_1=mean(hbl,1);
   sst_1=mean(sst,1);
   lat_1=lat(1,:);
  elseif ilevel==2,
   hbl_2=mean(hbl,1);
   sst_2=mean(sst,1);
   lat_2=lat(1,:);
  elseif ilevel==3,
   hbl_3=mean(hbl,1);
   sst_3=mean(sst,1);
   lat_3=lat(1,:);
  end

end

% HBL

lx=20; ly=20;
figure('units','centimeters','position', ...
           [0 0 lx ly],'paperpositionmode','auto')
hbl_1(hbl_1==0)=NaN;hbl_2(hbl_2==0)=NaN;hbl_3(hbl_3==0)=NaN;
pt1=plot(lat_1,-hbl_1); hold on;
pt2=plot(lat_2,-hbl_2); hold on;
pt3=plot(lat_3,-hbl_3); hold on;
set(gca,'fontsize',16)

%hl=line([0 0],[-4 10]); hold off
%set(hl,'Linewidth',2,'Linestyle','-','Color','r')

set(pt1,'Linewidth',2,'Linestyle',':','Color','k')
set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
set(pt3,'Linewidth',2,'Linestyle','-','Color','k')
legend('COARSE 36km','MEDIUM 12km','FINE 4km','Location','NorthEast')

set(gca,'Yscale','linear','Xscale','linear');
set(gca,'Xlim',[-4 10],'Ylim',[-48 -12]);
xlabel('Latitude')
ylabel('-HBL [m]')
set(gca,'fontsize',16)

if print_and_keep,
   outname=['hbl_prof_s',num2str(season)];
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

% SST

lx=20; ly=20;
figure('units','centimeters','position', ...
           [0 0 lx ly],'paperpositionmode','auto')
pt1=plot(lat_1,sst_1); hold on;
pt2=plot(lat_2,sst_2); hold on;
pt3=plot(lat_3,sst_3); hold on;
set(gca,'fontsize',16)

%hl=line([0 0],[-4 10]); hold off
%set(hl,'Linewidth',2,'Linestyle','-','Color','r')

set(pt1,'Linewidth',2,'Linestyle',':','Color','k')
set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
set(pt3,'Linewidth',2,'Linestyle','-','Color','k')
legend('COARSE 36km','MEDIUM 12km','FINE 4km','Location','SouthEast')

set(gca,'Yscale','linear','Xscale','linear');
set(gca,'Xlim',[-4 10],'Ylim',[22.5 27.5]);
xlabel('Latitude')
ylabel('SST [^oC]')
set(gca,'fontsize',16)

if print_and_keep,
   outname=['sst_prof_s',num2str(season)];
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

