%==============================================================================
% MAKE ALL SPECTRA PLOTS
%==============================================================================
clear all
%close all
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
%
eval(['!cp ',root,specname,'_KE_l',num2str(ilevel),'_Zint.mat spectrum_KE.mat'])
load spectrum_KE.mat
%
if window==3,
  istr=2; iend=length(amp);
else
  istr=1; iend=length(amp);
end
pt3=plot_spectrum(ktmp(istr:iend),amp(istr:iend)./dk(istr:iend),nan,nan,2); hold on;
set(gca,'fontsize',16)
set(pt3,'Linewidth',3,'Color','r');
%
set(gca,'Ylim',[-Inf Inf],'Xlim',[xmin xmax]);
xlabel('k [rad/m]')
ylabel('KE(10m) [m^3/s^2]')
set(gca,'fontsize',16)
%
x=[5000,3000,2000,1500,1000,700,500,400:-100:100,80:-10:10,5];
LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
if print_and_keep,
 outname=['KE_spectrum'];
 warning off
 eval(['print -painter -depsc2 ',outname,'.eps;'])
 warning on
 eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
 eval(['! mv ',outname,'.eps ',dirout_EPS]);
 eval(['! mv ',outname,'.jpg ',dirout_JPG]);
end

hold off

return
