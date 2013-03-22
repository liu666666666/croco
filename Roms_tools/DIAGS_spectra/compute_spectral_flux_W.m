function [amp,count,ktmp,dk]=compute_spectral_flux_diffv(ctl,model,dx,kchoice,window,lims,lit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [amp,count,ktmp]=compute_spectral_flux_diffv(ctl,dx,kchoice,window,lims,lit);
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
rho0=1025;

grd=rnt_gridload(model);
pm=grd.pm;pm=squeeze(pm(Lmin:Lmax,Mmin:Mmax));
pn=grd.pn;pn=squeeze(pn(Lmin:Lmax,Mmin:Mmax));

hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;h=grd.h;h0=h(1,1);method=grd.Method;
kmin=kchoice-1;kmax=min(N,kchoice+1);
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
  lit=[1:length(ctl.time)];
end

for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])
%
  u=rnt_loadvar_partialz(ctl,it,'u',kchoice,kchoice); % u is 3d
  u=perm(u2rho(perm(u)));
  u=squeeze(u(Lmin:Lmax,Mmin:Mmax,:));
  v=rnt_loadvar_partialz(ctl,it,'v',kchoice,kchoice); % v is 3d
  v=perm(v2rho(perm(v)));
  v=squeeze(v(Lmin:Lmax,Mmin:Mmax,:));
%
  dissu=zeros(Lmax-Lmin+1,Mmax-Mmin+1);
  dissv=zeros(Lmax-Lmin+1,Mmax-Mmin+1);

  sustr=rnt_loadvar(ctl,it,'sustr')./rho0;
  sustr=perm(u2rho(perm(sustr)));
  sustr=squeeze(sustr(Lmin:Lmax,Mmin:Mmax));
  svstr=rnt_loadvar(ctl,it,'svstr')./rho0;
  svstr=perm(v2rho(perm(svstr)));
  svstr=squeeze(svstr(Lmin:Lmax,Mmin:Mmax));

  dissu=sustr./dzr(:,:,2);;
  dissv=svstr./dzr(:,:,2);;
%
  [u,v,dissu,dissv]=windowing(u,v,dissu,dissv,window);

  [M,L]=size(u);
  if it==lit(1);
     tmpamp=zeros(M,L); 
  end
  fcoefdissu=fft2(dissu); fcoefdissv=fft2(dissv);
  fcoefu=fft2(u); fcoefv=fft2(v);
  tmpamp=tmpamp+real(conj(fcoefu).*fcoefdissu+conj(fcoefv).*fcoefdissv);

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
