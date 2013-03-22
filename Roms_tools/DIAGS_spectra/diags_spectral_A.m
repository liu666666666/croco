%==============================================================
% SPECTRAL KE ADVECTION
% convention manu with first dim being xsi and second eta
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================
close all
clear all

diags_params

loadmode = input('loadmode? (0/no 1/yes) ');
print_and_keep=1;
cff_scale=1.e4;
%===============================================================

lx=30; ly=15;
figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')
%

for ilevel=1:Nlev

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
  kmin=klist(1);kmax=klist(end);
  grd=rnt_gridload(model);
  pm=grd.pm;pm=squeeze(pm(Lmin:Lmax,Mmin:Mmax));
  pn=grd.pn;pn=squeeze(pn(Lmin:Lmax,Mmin:Mmax));
  %lon=grd.lonr;lon=squeeze(lon(Lmin:Lmax,Mmin:Mmax));
  %lat=grd.latr;lat=squeeze(lat(Lmin:Lmax,Mmin:Mmax));
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

  kk=0;
  [L,M]=size(pm);
  amp=zeros(1,min(L,M)/2);
  ampA=amp;

%
% HORIZONTAL ADVECTION
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_Ah_l',num2str(ilevel),'.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    for kchoice=klist,  %%%%%%%%%%%%%%%
      kk=kk+1;
      cffz=mean(mean(dzr(:,:,kk)./sum(dzr,3)));
      disp([' %%%%%% KCHOICE : ',num2str(kchoice),'   %%%%%%%%%'])
      [amp,count,ktmp,dk]=compute_spectral_flux_advh(ctl,model,dx,kchoice,window,lims,lit);
      ampA=ampA+cffz*amp;
    end %kchoice
    amp=ampA;
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_Ah_l',num2str(ilevel),'.mat'])
  end
  amptmp=amp;
%
% VERTICAL ADVECTION
%
  if loadmode==1,
    eval(['!cp ',root,'spectrum_KE_Av_l',num2str(ilevel),'.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    kk=0;
    for kchoice=klist,  %%%%%%%%%%%%%%%
      kk=kk+1;
      cffz=mean(mean(dzr(:,:,kk)./sum(dzr,3)));
      disp([' %%%%%% KCHOICE : ',num2str(kchoice),'   %%%%%%%%%'])
      [amp,count,ktmp,dk]=compute_spectral_flux_advw(ctl,model,dx,kchoice,window,lims,lit);
      ampA=ampA+cffz*amp;
    end %kchoice
    amp=ampA;
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,'spectrum_KE_Av_l',num2str(ilevel),'.mat'])
  end
  amp=amp+amptmp;

  K=ktmp(istr:end);
  A=cff_scale*amp(istr:end)./dk(istr:end);
  if filtamp
    A(2:end-1)=0.7*A(2:end-1)+0.15*(A(1:end-2)+A(3:end));
  end
  K=interp(K,4);
  A=interp(A,4);

  if ilevel==1,
    pt1=plot(K,A); hold on;
    set(pt1,'Linewidth',2,'Linestyle',':','Color','k')
  elseif ilevel==2,
    pt2=plot(K,A); hold on;
    set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
  else
    pt3=plot(K,A); hold on;
    set(pt3,'Linewidth',2,'Linestyle','-','Color','k')
  end
  set(gca,'fontsize',16)

end

hl=line([3e-6 2e-3],[0 0]); hold off
set(hl,'Linewidth',2,'Linestyle','-','Color','r')

%legend('COARSE 36km','MEDIUM 12km','FINE 4km','Location','NorthEast')

set(gca,'Yscale','linear','Xscale','log');
set(gca,'Ylim',[ymin ymax],'Xlim',[xmin xmax]);
xlabel('k [rad/m]')
[coef,expon] =strread(strrep(sprintf('%E',cff_scale),'E','#'),'%f#%f');
ylabel(['spectral flux [10^{-',num2str(expon),'} m^2/s^3]'])
set(gca,'fontsize',16)
%
x=[2000:-500:500,400:-100:100,80:-10:10,5];
LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
if print_and_keep,
  outname=['KE_A_spectra'];
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



