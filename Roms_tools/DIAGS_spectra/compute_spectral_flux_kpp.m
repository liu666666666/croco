function [amp,count,ktmp,dk]=compute_spectral_flux_kpp(ctl,model,dx,kchoice,window,lims,lit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [amp,count,ktmp]=compute_spectral_flux_kpp(ctl,dx,kchoice,window,lims,lit);
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
Kv_constant=1;

grd=rnt_gridload(model);
pm=grd.pm;pm=squeeze(pm(Lmin:Lmax,Mmin:Mmax));
pn=grd.pn;pn=squeeze(pn(Lmin:Lmax,Mmin:Mmax));

kmin=kchoice-1;kmax=kchoice+1;
hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;h=grd.h;h0=h(1,1);method=grd.Method;
zr=zlevs(h0,0,thetas,thetab,hc,N,'r');
zw=zlevs(h0,0,thetas,thetab,hc,N,'w');
pz=zeros(size(zr));
pz(2:N-1)=0.5.*(zr(3:N)-zr(1:N-2));pz(1)=pz(2);pz(N)=pz(N-1);
pz=1./pz; pz=pz(kmin:kmax);zraux=zr(kmin:kmax);zwaux=zw(kmin:kmax+1);

if nargin==5
%  it1=1; it2=length(ctl.time);
  lit=[1:length(ctl.time)];
end

for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])
  if Kv_constant==1,
   Kv=1.e-4*ones(Lmax-Lmin+1,Mmax-Mmin+1,2);
  else
   Kv=rnt_loadvar_partialz(ctl,it,'AKv',kchoice,kchoice+1);
   Kv=Kv(Lmin:Lmax,Mmin:Mmax,:);
  end

  u=rnt_loadvar_partialz(ctl,it,'u',kchoice-1,kchoice+1); % u is 3d
  u=perm(u2rho_3d(perm(u))); 
  u=u(Lmin:Lmax,Mmin:Mmax,:);
  dudz=partialz_r2w(u,zraux);
  dissu=Kv.*dudz;
  dissu=partialz_w2r(dissu,squeeze(zwaux));
  u=squeeze(u(:,:,2));dudz=squeeze(mean(dudz,3));         % v is 2d now

  v=rnt_loadvar_partialz(ctl,it,'v',kchoice-1,kchoice+1); % v is 3d
  v=perm(v2rho_3d(perm(v))); 
  v=v(Lmin:Lmax,Mmin:Mmax,:);
  dvdz=partialz_r2w(v,zraux);
  dissv=Kv.*dvdz;
  dissv=partialz_w2r(dissv,squeeze(zwaux));
  v=squeeze(v(:,:,2));dvdz=squeeze(mean(dvdz,3));         % v is 2d now

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
