%==============================================================
% SPECTRAL KE BUDGET
% convention manu with first dim being xsi and second eta
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================
close all
clear all

diags_params

print_and_keep=0;
windwork=0;
ilevel=1;
%================================================================

  eval(['!cp ',root,specname,'_KE_ALL_l',num2str(ilevel),'.mat spectrum_KE.mat'])
  load spectrum_KE.mat

%===============        PLOT        ================================
lx=30; ly=15;
 figure('units','centimeters','position', ...
        [0 0 lx ly],'paperpositionmode','auto')

% I
K=ktmp(1:end);
A=cff_scale*ampI(1:end)./dk(1:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
K=interp(K,4);
A=interp(A,4);
pt1=plot(K,A); hold on;
set(gca,'fontsize',16)
% A
A=cff_scale*ampA(1:end)./dk(1:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
A=interp(A,4); 
pt2=plot(K,A); hold on;
set(gca,'fontsize',16)
% DH
A=10*cff_scale*ampDH(1:end)./dk(1:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
A=interp(A,4);
pt3=plot(K,A); hold on;
set(gca,'fontsize',16) 
% DV
A=cff_scale*(ampDV(1:end)+ampP(1:end))./dk(1:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
A=interp(A,4); 
pt4=plot(K,A); hold on;
set(gca,'fontsize',16)
% P
A=cff_scale*ampP(1:end)./dk(1:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
A=interp(A,4); 
% pt5=plot(K,A); hold on;
% set(gca,'fontsize',16)
% W
if windwork & kmax==N,
 A=cff_scale*ampW(1:end)./dk(1:end);
 if filtamp
  A=hanning_1d(A,filtcff);
 end
 A=interp(A,4); 
 pt22=plot(K,A,'g'); hold on;
 %set(gca,'fontsize',16)
end
% TOTAL
A=cff_scale*(ampI(1:end)+ampA(1:end)+ampDH(1:end)+ampDV(1:end)... 
            +ampP(1:end)-ampT(1:end) )./dk(1:end);
          %   +ampP(1:end)-ampT(1:end))./dk(1:end);
if filtamp
 A=hanning_1d(A,filtcff);
end
A=interp(A,4); 
%pt6=plot(K,A); hold on;
%set(gca,'fontsize',16)
%
% TENDENCY
A=cff_scale*(ampT(1:end))./dk(1:end);
if filtamp
 A=hanning_1d(A,filtcff);
end
A=interp(A,4); 
%pt7=plot(K,A,'r'); hold on;
%set(gca,'fontsize',16)
%

%
% FINALIZE PLOT
%
line([3e-6 2e-3],[0 0],'Color','r'); hold off;
set(pt1,'Linewidth',2,'Linestyle','-' ,'Color','k')
set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
set(pt3,'Linewidth',2,'Linestyle','--' ,'Color',.7*[1 1 1])
set(pt4,'Linewidth',2,'Linestyle','-' ,'Color',.7*[1 1 1])

legend('I','A','D_H x10','D_V+P','Location','SouthEast')

set(gca,'Yscale','linear','Xscale','log');
set(gca,'Ylim',[ymin ymax],'Xlim',[xmin xmax]);
xlabel('k [rad/m]')
[coef,expon] =strread(strrep(sprintf('%E',cff_scale),'E','#'),'%f#%f');
ylabel(['KE tendencies [10^{-',num2str(expon),'} m^2/s^3]'])
set(gca,'fontsize',16)
%
x=[2000,1500,1000,700,500,400:-100:100,80:-10:10,5];
LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
if print_and_keep,
 outname=['KE_budget_Zint'];
 warning off
 eval(['print -painter -depsc2 ',outname,'.eps;'])
 warning on
 eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])

 eval(['! mv ',outname,'.eps ',dirout_EPS]);
 eval(['! mv ',outname,'.jpg ',dirout_JPG]);
 eval(['! rm -f ',outname,'.eps']);
 eval(['! rm -f ',outname,'.jpg']);
end


