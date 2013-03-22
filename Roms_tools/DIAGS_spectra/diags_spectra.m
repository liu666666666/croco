%==============================================================
% COMPUTE KE SPECTRA OF EMBEDDED GRID SOLUTIONS
% convention manu with first dim being xsi and second eta
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================
clear all
close all

loadmode = input('loadmode? (0/no 1/yes) ');

diags_params
nestrange=[Nlev:-1:1];
npts=500; % spectrum
print_and_keep=0;
%==============================================================

%---------------- LOOP ON NESTED LEVELS ------------------
for ilevel=nestrange

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

  kk=0;
  [L,M]=size(pm);
  if npts>0 & npts>L & npts>M,
   M=npts; L=npts;
  end

  if window==3,
   amp=zeros(1,min(L,M));
  else
   amp=zeros(1,min(L,M)/2);
  end
  ampM=amp;

%
% Compute KE SPECTRUM ................
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_l',num2str(ilevel),'_Zint.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    for kchoice=klist
      kk=kk+1;
      cffz=mean(mean(dzr(:,:,kk)./sum(dzr,3)));
      disp([' %%%%%% KCHOICE : ',num2str(kchoice),'   %%%%%%%%%'])
      [amp,count,ktmp,dk]=compute_power_spectrum(ctl,dx,kchoice,window,npts,lims,lit);
      ampM=ampM+cffz*amp;
    end
    amp=ampM;
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_l',num2str(ilevel),'_Zint.mat'])
  end

%
% PLOT spectrum
%
  if window==3,
    istr=2; iend=length(amp);
  else
    istr=1; iend=length(amp);
  end
  K=ktmp(istr:iend);
  A=amp(istr:iend)./dk(istr:iend);

  lx=30; ly=15;
  figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')

  h1=plot_spectrum(K,A,NaN,NaN,2); hold on;
  set(gca,'fontsize',16)
  set(h1,'Linewidth',3,'Color','r');
%
  set(gca,'Ylim',[-Inf Inf],'Xlim',[xmin xmax]);
  xlabel('k [rad/m]')
  ylabel('KE(10m) [m^3/s^2]')
  set(gca,'fontsize',16)
%
  x=[5000,3000,2000,1500,1000,700,500,400:-100:100,80:-10:10,5];
  LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
  hold off
  if print_and_keep,
   outname=['KE_spectrum'];
   warning off
   eval(['print -painter -depsc2 ',outname,'.eps;'])
   warning on
   eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
   eval(['! mv ',outname,'.eps ',dirout_EPS]);
   eval(['! mv ',outname,'.jpg ',dirout_JPG]);
  end

% Compensated spectrum
  cff_scale=1.e10;
  B=cff_scale*A.*K.^2;
  figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')

  h1=plot(K,B); hold on;
  set(gca,'fontsize',16)
  set(h1,'Linewidth',3,'Color','r');
  set(gca,'Yscale','linear','Xscale','log');
%
  hl=line([xmin xmax],[0.55 0.55]); hold off
  set(hl,'Linewidth',1,'Linestyle','-','Color','k')

  set(gca,'Ylim',[-Inf Inf],'Xlim',[xmin xmax]);
  xlabel('k [rad/m]')
  [coef,expon] =strread(strrep(sprintf('%E',cff_scale),'E','#'),'%f#%f');
  ylabel(['k^2KE (10m) [10^{-',num2str(expon),'} m/s^2]'])
  set(gca,'fontsize',16)
%
  x=[5000,3000,2000,1500,1000,700,500,400:-100:100,80:-10:10,5];
  LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
  if print_and_keep,
   outname=['KE_spectrum_compensated'];
   warning off
   eval(['print -painter -depsc2 ',outname,'.eps;'])
   warning on
   eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
   eval(['! mv ',outname,'.eps ',dirout_EPS]);
   eval(['! mv ',outname,'.jpg ',dirout_JPG]);
  end
  hold off

%
% TRACER SPECTRUM .............
%
%  [amp,count,ktmp]=compute_power_spectrum_tracer(ctl,dx,kchoice,window,lims,lit,'temp');
%  dk=ktmp(10)-ktmp(9);
%  save spectrum_T amp count ktmp dk
%  eval(['!mv spectrum_T.mat ',root,specname,'_T_l',num2str(ilevel),'.mat'])

end % ilevel ----------------------

