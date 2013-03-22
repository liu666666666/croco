%==============================================================
% SPECTRAL KE ADVECTION
% convention manu with first dim being xsi and second eta
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================
%close all
clear all

diags_params

print_and_keep=0;
cff_scale=1.e4;
%===============================================================

lx=30; ly=15;
figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')

for ilevel=1:Nlev;

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
  it1=lit(1); it2=lit(end);
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
  ampT=amp;

  for kchoice=klist,  %%%%%%%%%%%%%%%
    kk=kk+1;
    cffz=mean(mean(dzr(:,:,kk)./sum(dzr,3)));
    disp([' %%%%%% KCHOICE : ',num2str(kchoice),'   %%%%%%%%%'])

    for it=[it1 it2]

     disp(['Treating rec ' num2str(it)])
     u=rnt_loadvar_partialz(ctl,it,'u',kchoice,kchoice); % u is 3d
     v=rnt_loadvar_partialz(ctl,it,'v',kchoice,kchoice); % v is 3d
     u=perm(u2rho(perm(u)));
     u=squeeze(u(Lmin:Lmax,Mmin:Mmax));
     v=perm(v2rho(perm(v)));
     v=squeeze(v(Lmin:Lmax,Mmin:Mmax));
     u=u-mean(mean(u));
     v=v-mean(mean(v));
  
     [u,v,u,v]=windowing(u,v,u,v,window);

     [M,L]=size(u);

     fcoefu=fft2(u); fcoefv=fft2(v);
     tmpamp=0.5*(real(conj(fcoefu).*fcoefu+conj(fcoefv).*fcoefv));

     if it==it1,
       amp_ti=tmpamp;
     else
       amp_tf=tmpamp;
     end

    end % temporal loop

    tmpamp=(amp_tf-amp_ti)/((it2-it1)*86400);

    tmpamp=1/(L*M)^2*tmpamp;               % parseval equality
    tmpamp=reorganize_fft2d(tmpamp);       % reorganize to be centered
    method=2;
    [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,method); % 1D spectra
    ampT=ampT+cffz*amp;
  end % kchoice
%
  amp=ampT;
%  

  K=ktmp(1:end);
  A=cff_scale*amp(1:end)./dk(1:end);
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

end %ilevel

hl=line([3e-6 2e-3],[0 0]); hold off
set(hl,'Linewidth',2,'Linestyle','-','Color','r')

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
  outname=['KE_T_spectra'];
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


return

