function [amp,count,ktmp,dk]=compute_spectral_flux_advh(ctl,model,dx,kchoice,window,lims,lit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [amp,count,ktmp]=compute_spectral_flux_adv(ctl,dx,kchoice,window,lims,lit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  this function computes velocity power spectrum for a netcdf file series. 
%%%% inputs: 
%           ctl    : structure with netcdf files loaded
%           dx     : meshgrid size (!!! supposed to be uniform !!!)
%           kchoice: sigma level to be processed (z level not implemented)
%           window : windowing option. There are 2 options. 
%                       apply a filter (rectangular 1, triangular 2 or hanning 3)
%                       apply a 2d dct    (0) or 1d dct (xsi dir, -1)
%           lims   : limits of the domain to be processed (4 elements vector with
%                       lims(1) and lims(2) mins and max in xsi direction. 
%           t1,t2  : time limits in the series. Default: 1 and length(ctl.time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
pi=3.1415926535;

grd=rnt_gridload(model);
pm=grd.pm;pm=squeeze(pm(Lmin:Lmax,Mmin:Mmax));
pn=grd.pn;pn=squeeze(pn(Lmin:Lmax,Mmin:Mmax));

hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;h=grd.h;h0=h(1,1);method=grd.Method;
kmin=kchoice-1;kmax=min(N,kchoice+1);
zr=zlevs(h0,0,thetas,thetab,hc,N,'r');
zw=zlevs(h0,0,thetas,thetab,hc,N,'w');
pz=zeros(size(zr));
pz(2:N-1)=0.5.*(zr(3:N)-zr(1:N-2));pz(1)=pz(2);pz(N)=pz(N-1);
pz=1./pz; pz=pz(kmin:kmax);zraux=zr(kmin:kmax);zwaux=zw(kmin:kmax+1);

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


if nargin==5
%  it1=1; it2=length(ctl.time);
  lit=[1:length(ctl.time)];
end

for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])
  w=rnt_loadvar_partialz(ctl,it,'w',kchoice,kchoice);     % w is 2d
  if isempty(w)
  w=mean(rnt_loadvar_partialz(ctl,it,'omega',kchoice,kchoice+1),3);     % w is 2d
  end
  w=squeeze(w(Lmin:Lmax,Mmin:Mmax));

  u=rnt_loadvar_partialz(ctl,it,'u',kmin,kmax); % u is 3d
  u=perm(u2rho(perm(u)));
  u=squeeze(u(Lmin:Lmax,Mmin:Mmax,:));
  k=2;
  if kchoice==N
    advu=(u(:,:,k)-u(:,:,k-1))./dzw(:,:,k);
  else
    advu=(u(:,:,k+1)-u(:,:,k  ))./dzw(:,:,k+1);
  end
  %advu=partialz(u,pz);
  u=squeeze(u(:,:,k));
  advu=-w.*advu;

  v=rnt_loadvar_partialz(ctl,it,'v',kmin,kmax); % v is 3d
  v=perm(v2rho(perm(v)));
  v=squeeze(v(Lmin:Lmax,Mmin:Mmax,:));
  if kchoice==N
    advv=(v(:,:,k)-v(:,:,k-1))./dzw(:,:,k);
  else
    advv=(v(:,:,k+1)-v(:,:,k  ))./dzw(:,:,k+1);
  end
  %advv=partialz(v,pz);
  v=squeeze(v(:,:,k));
  advv=-w.*advv;

  [u,v,advu,advv]=windowing(u,v,advu,advv,window);

  [M,L]=size(u);
  if it==lit(1);
     tmpamp=zeros(M,L); 
  end
  fcoefadvu=fft2(advu); fcoefadvv=fft2(advv);
  fcoefu=fft2(u); fcoefv=fft2(v);
  tmpamp=tmpamp+real(conj(fcoefu).*fcoefadvu+conj(fcoefv).*fcoefadvv);

end % temporal loop

%----------------------------------------------------------------
% if check parseval's relationship
% lhswt=mean(mean(w.*T));
% rhswt=1/(Lu*Mu)^2*sum(sum(tmpampwT));
%
% do the reorganization and processing of results.
tmpamp=1/(length(lit)*(L*M)^2)*tmpamp; % parseval equality
tmpamp=reorganize_fft2d(tmpamp);       % reorganize to be centered

%
% Changes to get 1D spectra
%
new=1;
if new==1
  method=2;
  [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,method);
else
  Lh=L/2;Mh=M/2;i2L=1/L^2;i2M=1/M^2;
  kk=repmat(perm([1:L])-Lh,[1 M]); ll=repmat([1:M]-Mh,[L 1]);
  Kh=i2L*kk.^2+i2M*ll.^2; Kh=6.28/(dx)*sqrt(Kh);Kh=perm(Kh);
  leng=L;
  aux=leng/Kh(M,L);ktmp=[1:leng]/aux;ktmp=perm(ktmp);
  amp=zeros(leng,1);count=zeros(leng,1);

  ik = min(round(aux*Kh)+1,leng);
  for ind=1:leng
    II=find(ik==ind);
    amp(ind)=sum(tmpamp(II));
    count(ind)=length(II);
  end
end

return
