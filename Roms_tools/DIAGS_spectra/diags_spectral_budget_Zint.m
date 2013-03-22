%==============================================================
% SPECTRAL KE BUDGET
% convention manu with first dim being xsi and second eta
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================
%close all
clear all

diags_params

loadmode = input('loadmode? (0/no 1/yes) ')
print_and_keep=0;
windwork=0;
checkplot=1;

ilevel=1;
%================================================================

dx=dx_lev(ilevel);
if Nlev==1,
  lims=lims_lev;
  model=model0;
else
  lims=squeeze(lims_lev(:,ilevel+1));
  model=[model0 num2str(ilevel)];
end
if mod(lims(2)-lims(1)+1,2)==1
   lims(2)=lims(2)+1;
end
if mod(lims(4)-lims(3)+1,2)==1
   lims(4)=lims(4)+1;
end
ctl=ctlload(root(2:end),his,model,0,lastfile_indx,Nrec);
ctlhis=ctlload(root(2:end),his,model,0,lastfile_indx,Nrec);

Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
grd=rnt_gridload(model);
pm=grd.pm;pm=squeeze(pm(Lmin:Lmax,Mmin:Mmax));
pn=grd.pn;pn=squeeze(pn(Lmin:Lmax,Mmin:Mmax));
%lon=grd.lonr;lon=squeeze(lon(Lmin:Lmax,Mmin:Mmax));
%lat=grd.latr;lat=squeeze(lat(Lmin:Lmax,Mmin:Mmax));
f=grd.f;f=squeeze(f(Lmin:Lmax,Mmin:Mmax));

kmin=klist(1);kmax=klist(end);
hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;h=grd.h;h0=h(1,1);method=grd.Method;
zr=zlevs(h,0,thetas,thetab,hc,N,'r');
zw=zlevs(h,0,thetas,thetab,hc,N,'w');
zr=permute(zr,[2,3,1]);
zw=permute(zw,[2,3,1]);
dzw=zeros(size(zw));
dzw(:,:,2:N)=zr(:,:,2:N)-zr(:,:,1:N-1);
dzw(:,:,1)=dzw(:,:,2);
dzw=dzw(Lmin:Lmax,Mmin:Mmax,kmin:kmax);
dzr=zeros(size(zr));
dzr(:,:,2:N)=zw(:,:,3:N+1)-zw(:,:,2:N);
dzr(:,:,1)=dzr(:,:,2);
dzr=dzr(Lmin:Lmax,Mmin:Mmax,kmin:kmax);

if loadmode,
  eval(['!cp ',root,specname,'_KE_ALL_l',num2str(ilevel),'.mat spectrum_KE.mat'])
  load spectrum_KE.mat
else

 [L,M]=size(pm);
 if window==3,
  amp=zeros(1,min(L,M));
 else
  amp=zeros(1,min(L,M)/2);
 end
 ampI=amp;
 ampA=amp;
 ampP=amp;
 ampDH=amp;
 ampDV=amp;
 ampT=amp;
 ampCO=amp;
 kk=0;

 for kchoice=klist,  %%%%%%%%%%%%%%%

  kk=kk+1;
  cffz=mean(mean(dzr(:,:,kk)./sum(dzr,3)));
  disp([' %%%%%% KCHOICE : ',num2str(kchoice),'   %%%%%%%%%'])

  disp(' BAROCLINIC CONVERSION ')
  [amp,count,ktmp,dk]=compute_spectral_flux_wb(ctl,model,dx,kchoice,window,lims,lit);
  ampI=ampI+cffz*amp; 
  ampIk=amp;%save for 3D pressure work computation
%
  disp('  HORIZONTAL ADVECTION ')
  [amp,count,ktmp,dk]=compute_spectral_flux_advh(ctl,model,dx,kchoice,window,lims,lit);
  ampAH=amp;
%
  disp(' VERTICAL ADVECTION ')
  [amp,count,ktmp,dk]=compute_spectral_flux_advw(ctl,model,dx,kchoice,window,lims,lit);
  ampA=ampA+cffz*(amp+ampAH);
%
  disp(' HORIZONTAL DISSIPATION ')
  [amp,count,ktmp,dk]=compute_spectral_flux_diffh(ctl,model,dx,kchoice,window,lims,lit);
  ampDH=ampDH+cffz*amp;
%
  disp(' VERTICAL DISSIPATION ')
  [amp,count,ktmp,dk]=compute_spectral_flux_diffv(ctl,model,dx,kchoice,window,lims,lit);
  ampDV=ampDV+cffz*amp;
%
  disp(' 3D PRESSURE WORK ')
  [amp,count,ktmp,dk]=compute_spectral_flux_P(ctl,model,dx,kchoice,window,lims,lit);
  amp=amp-ampIk;
  ampP=ampP+cffz*amp;
%
  disp(' CORIOLIS ')
  [amp,count,ktmp,dk]=compute_spectral_flux_cor(ctl,model,dx,kchoice,window,lims,lit);
  ampCO=ampCO+cffz*amp;
%
  disp(' TIME TENDENCY ')
  compT=2;
  if compT==1,
   [amp,count,ktmp,dk]=compute_spectral_flux_T(ctl,model,dx,kchoice,window,lims,lit);
  else
   it1=lit(1); it2=lit(end);
   for it=[it1-1 it1 it2 it2+1]
    disp(['Treating rec ' num2str(it)])
    u=rnt_loadvar_partialz(ctlhis,it,'u',kchoice,kchoice); % u is 3d
    v=rnt_loadvar_partialz(ctlhis,it,'v',kchoice,kchoice); % v is 3d
    u=perm(u2rho(perm(u)));
    u=squeeze(u(Lmin:Lmax,Mmin:Mmax));
    v=perm(v2rho(perm(v)));
    v=squeeze(v(Lmin:Lmax,Mmin:Mmax));
    [u,v,u,v]=windowing(u,v,u,v,window);
    ke=0.5*(u.^2+v.^2); 
    [M,L]=size(u);
    fcoefu=fft2(u); fcoefv=fft2(v);
    tmpamp=0.5*(real(conj(fcoefu).*fcoefu+conj(fcoefv).*fcoefv));
    if it==it1-1,
      amp_tim=tmpamp;
      ke_tim=ke;
    elseif it==it1,
      amp_ti=tmpamp;
      ke_ti=ke;
      if it1==it2; 
       amp_tf=tmpamp;
       ke_tf=ke;
      end
    elseif it==it2,
      amp_tf=tmpamp;
      ke_tf=ke;
    elseif it==it2+1,
      amp_tfp=tmpamp;
      ke_tfp=ke;
    end
   end % temporal loop
   amp_tf=0.5*(amp_tf+amp_tfp);
   amp_ti=0.5*(amp_ti+amp_tim);
   tmpamp=(amp_tf-amp_ti)/((it2-it1+1)*86400);
   tmpamp=1/(L*M)^2*tmpamp;     
   tmpamp=reorganize_fft2d(tmpamp);  
   [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,2);
   ket=0.5*((ke_tf+ke_tfp)-(ke_ti+ke_tim))./((it2-it1+1)*86400);
  end
%
   ampT=ampT+cffz*amp;
%
 end % kchoice

 if windwork & kmax==N,
  disp(' WIND WORK  ')
  for it=lit
   disp(['Treating rec ' num2str(it)])
   sustr=rnt_loadvar(ctlhis,it,'sustr')./(rho0*5);
   svstr=rnt_loadvar(ctlhis,it,'svstr')./(rho0*5);
   us=rnt_loadvar_partialz(ctlhis,it,'u',N,N);
   vs=rnt_loadvar_partialz(ctlhis,it,'v',N,N);
   sustr=perm(u2rho_2d(perm(sustr)));sustr=cut_var(sustr,lims);
   svstr=perm(v2rho_2d(perm(svstr)));svstr=cut_var(svstr,lims);
   us=perm(u2rho_2d(perm(us)));us=cut_var(us,lims);
   vs=perm(u2rho_2d(perm(vs)));vs=cut_var(vs,lims);
   us=us-mean(mean(us));vs=vs-mean(mean(vs));
   sustr=sustr-mean(mean(sustr)); svstr=svstr-mean(mean(svstr));
   [us,vs,sustr,svstr]=windowing(us,vs,sustr,svstr,window);
   fcoefus=fft2(us); fcoefvs=fft2(vs);
   fcoeftmpu=fft2(sustr); fcoeftmpv=fft2(svstr);
   tmpamp=tmpamp+real(conj(fcoefus).*fcoeftmpu+conj(fcoefvs).*fcoeftmpv);
  end
  tmpamp=1/((it2-it1+1)*(L*M)^2)*tmpamp;
  tmpamp=reorganize_fft2d(tmpamp);       % reorganize to be centered
  [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,2);
%
  cffz=mean(mean(1./sum(dzr,3)));
  ampW=cffz*amp;
 end % windwork

 save spectrum_KE.mat ampI ampA ampDH ampDV ampP ampT ampCO ktmp dk
 eval(['!mv spectrum_KE.mat ',root,specname,'_KE_ALL_l',num2str(ilevel),'.mat'])

end %loadmode

%===============        PLOT        ================================
lx=30; ly=15;
figure('units','centimeters','position', ...
        [0 0 lx ly],'paperpositionmode','auto')
%
% TOTAL
K=ktmp(istr:end);
A=cff_scale*(ampI(istr:end)+ampA(istr:end)+ampDH(istr:end)+ampDV(istr:end)...
            +ampCO(istr:end)+ampP(istr:end)-ampT(istr:end) )./dk(istr:end);
ATOT=A;
if filtamp
 A=hanning_1d(A,filtcff);
end
if supersamp,
 K=interp(K,4);
 A=interp(A,4);
 K=K(4:end);
 A=A(4:end);
end
if checkplot
 pt6=plot(K,A,'g'); hold on;
 set(gca,'fontsize',16)
end
% I
A=cff_scale*ampI(istr:end)./dk(istr:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
if supersamp,
 A=interp(A,4);
 A=A(4:end);
end
pt1=plot(K,A); hold on;
set(gca,'fontsize',16)
% A
A=cff_scale*ampA(istr:end)./dk(istr:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
if supersamp,
 A=interp(A,4);
 A=A(4:end);
end
pt2=plot(K,A); hold on;
set(gca,'fontsize',16)
% DH
A=10*cff_scale*ampDH(istr:end)./dk(istr:end);
if filtamp
  A=hanning_1d(A,filtcff);
end
if supersamp,
 A=interp(A,4);
 A=A(4:end);
end
pt3=plot(K,A); hold on;
set(gca,'fontsize',16) 
% CO
A=cff_scale*(ampCO(istr:end))./dk(istr:end);
ACOR=A;
if filtamp
  A=hanning_1d(A,filtcff);
end
if supersamp,
 A=interp(A,4);
 A=A(4:end);
end
if checkplot
 pt33=plot(K,A,'b'); hold on;
 set(gca,'fontsize',16)
end
% P
A=cff_scale*ampP(istr:end)./dk(istr:end);
P=A;
A=A+ACOR;
if filtamp
  A=hanning_1d(A,filtcff);
end
if supersamp,
 A=interp(A,4);
 A=A(4:end);
end
pt5=plot(K,A); hold on;
set(gca,'fontsize',16)
% DV
A=cff_scale*(ampDV(istr:end))./dk(istr:end);
%A=A-ATOT;
if filtamp
  A=hanning_1d(A,filtcff);
end
if supersamp,
 A=interp(A,4);
 A=A(4:end);
end
pt4=plot(K,A); hold on;
set(gca,'fontsize',16)
% W
if windwork & kmax==N,
 A=cff_scale*ampW(istr:end)./dk(istr:end);
 if filtamp
  A=hanning_1d(A,filtcff);
 end
 if supersamp,
  A=interp(A,4);
  A=A(4:end);
 end
 pt22=plot(K,A,'g'); hold on;
 %set(gca,'fontsize',16)
end
% TENDENCY
A=cff_scale*(ampT(istr:end))./dk(istr:end);
if filtamp
 A=hanning_1d(A,filtcff);
end
if supersamp,
 A=interp(A,4);
 A=A(4:end);
end
if checkplot
 pt7=plot(K,A,'r'); hold on;
 set(gca,'fontsize',16)
end
%
% FINALIZE PLOT
%
line([3e-6 2e-3],[0 0],'Color','r'); hold off;
set(pt1,'Linewidth',2,'Linestyle','-' ,'Color','k')
set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
set(pt3,'Linewidth',2,'Linestyle','--' ,'Color',.7*[1 1 1])
set(pt4,'Linewidth',2,'Linestyle','-' ,'Color',.7*[1 1 1])
set(pt5,'Linewidth',2,'Linestyle',':' ,'Color',.7*[1 1 1])
if checkplot
 legend('TOT','I','A','D_H (x10)','Cor','P','D_V','T','Location','SouthEast')
else
 legend('I','A','D_H (x10)','D_V','P','Location','SouthEast')
end
if windwork & kmax==N,
  set(pt6,'Linewidth',2,'Linestyle','--' ,'Color','b')
  legend('I','A','D_H (x10)','D_V','P','W','Location','SouthEast')
end

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
 outname=['KE_budget_spectra_Zint'];
 warning off
 eval(['print -painter -depsc2 ',outname,'.eps;'])
 warning on
 eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
 eval(['! mv ',outname,'.eps ',dirout_EPS]);
 eval(['! mv ',outname,'.jpg ',dirout_JPG]);
 eval(['! rm -f ',outname,'.eps']);
 eval(['! rm -f ',outname,'.jpg']);
end



