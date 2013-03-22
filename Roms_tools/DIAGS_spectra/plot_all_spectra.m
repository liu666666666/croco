%==============================================================================
% MAKE ALL SPECTRA PLOTS
%==============================================================================
clear all
close all
%==============================================================================

diags_params

print_and_keep=0;
ilevel=1;
%==============================================================================
%
% KE SPECTRUM .............
%
lx=30; ly=15;
figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')
%
for ilevel=Nlev:-1:1
%
  eval(['!cp ',root,'spectrum_KE_l',num2str(ilevel),'.mat spectrum_KE.mat'])
  load spectrum_KE.mat
%
  if window==3,
    istr=2; iend=length(amp);
  else
    istr=1; iend=length(amp);
  end
  amp(2:end-1)=0.5*amp(2:end-1)+0.25*(amp(1:end-2)+amp(3:end));
  if ilevel==1,
    pt1=plot_spectrum(ktmp(istr:iend),amp(istr:iend)./dk(istr:iend),nan,nan,ilevel); hold on;
  elseif ilevel==2,
    pt2=plot_spectrum(ktmp(istr:iend),amp(istr:iend)./dk(istr:iend),nan,nan,ilevel); hold on;
  else
    pt3=plot_spectrum(ktmp(istr:iend),amp(istr:iend)./dk(istr:iend),nan,nan,ilevel); hold on;
  end
  set(gca,'fontsize',16)
%
end

%legend([pt1,pt2,pt3],'COARSE 36km','MEDIUM 12km','FINE 4km','Location','NorthEast')

set(gca,'Ylim',[5e-5 5e3],'Xlim',[xmin xmax]);
xlabel('Wavenumber [rad/m]')
ylabel('KE (10m) [m^3/s^2]')
set(gca,'fontsize',16)
%
x=[3000,1500,1000,700,500,400:-100:100,80:-10:10,5];
LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
if print_and _keep,
 outname=['KE_spectra_s',num2str(season)];
 warning off
 eval(['print -painter -depsc2 ',outname,'.eps;'])
 orient landscape
 eval(['print -dpdf ',outname,'.pdf;'])
 warning on
 eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
%
 eval(['! mv ',outname,'.eps ',dirout_EPS]);
 eval(['! mv ',outname,'.jpg ',dirout_JPG]);
 eval(['! mv ',outname,'.pdf ',dirout_PDF]);
 eval(['! rm -f ',outname,'.eps']);
 eval(['! rm -f ',outname,'.jpg']);
 eval(['! rm -f ',outname,'.pdf']);
end
%=========================================================================

cff_scale=1/6.e-9;
lx=30; ly=15;
figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')
%
for ilevel=Nlev:-1:1
  eval(['!cp ',root,specname,'_KE_l',num2str(ilevel),'.mat spectrum_KE.mat'])
  load spectrum_KE.mat
  if window==3,
    istr=2; iend=length(amp);
  else
    istr=1; iend=length(amp);
  end
  A=amp(istr:iend)./dk(istr:iend);
  K=ktmp(istr:iend);
  A=hanning_1d(A,3);
  %A=interp(A,4);
  %K=interp(K,4);
  B=cff_scale*A.*K.^2;
  if ilevel==1,
    pt1=plot(K,B,'linewidth',2,'color','b'); hold on;
  elseif ilevel==2,
    pt2=plot(K,B,'linewidth',2,'color','r'); hold on;
  else
    pt3=plot(K,B,'linewidth',2,'color','g'); hold on;
  end
  set(gca,'fontsize',16)
end
%legend([pt1,pt2,pt3],'COARSE 36km','MEDIUM 12km','FINE 4km','Location','NorthEast')
hl=line([3e-6 2e-3],[1 1]); hold off
set(hl,'Linewidth',1,'Linestyle','-','Color','k')
set(gca,'Yscale','linear','Xscale','log');
set(gca,'Ylim',[-Inf Inf],'Xlim',[xmin xmax]);
xlabel('Wavenumber [rad/m]')
ylabel('k^2KE (10m) [m/s^2]')
x=[5000,3000,1500,1000,700,500,400:-100:100,80:-10:10,5];
LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
outname=['KE_spectra_compensated'];
warning off
eval(['print -painter -depsc2 ',outname,'.eps;'])
orient landscape
eval(['print -dpdf ',outname,'.pdf;'])
warning on
eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
%
eval(['! mv ',outname,'.eps ',dirout_EPS]);
eval(['! mv ',outname,'.jpg ',dirout_JPG]);
eval(['! mv ',outname,'.pdf ',dirout_PDF]);
eval(['! rm -f ',outname,'.eps']);
eval(['! rm -f ',outname,'.jpg']);
eval(['! rm -f ',outname,'.pdf']);


%=========================================================================

return

%
% TRACER SPECTRUM .............
%
lx=30; ly=15;
figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')
%
for ilevel=Nlev:-1:1
%
  eval(['!cp ',root,specname,'_T_l',num2str(ilevel),'.mat spectrum_T.mat'])
  load spectrum_T.mat
%
  [figt,pt]=plot_spectrum(ktmp(3:end),amp(3:end)/dk,0,nan,nan,ilevel);
  set(gca,'fontsize',16)
%
end 
set(gca,'Ylim',[1e-3 5e4],'Xlim',[xmin xmax]);
xlabel('k [rad/m]')
ylabel('T (10m) [^oC^2m]')
set(gca,'fontsize',16)
%
x=[1200,1000,700,500,400:-100:100,80:-10:10,5];
LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
if print_and_keep,
  outname=['KE_spectra_s',num2str(season)];
  warning off
  eval(['print -painter -depsc2 ',outname,'.eps;'])
  eval(['print -dpdf ',outname,'.pdf;'])
  warning on
  eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])

  eval(['! mv ',outname,'.eps ',dirout_EPS]);
  eval(['! mv ',outname,'.jpg ',dirout_JPG]);
  eval(['! mv ',outname,'.pdf ',dirout_PDF]);
  eval(['! rm -f ',outname,'.eps']);
  eval(['! rm -f ',outname,'.jpg']);
  eval(['! rm -f ',outname,'.pdf']);
end


