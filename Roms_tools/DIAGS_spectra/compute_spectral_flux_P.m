function [amp,count,ktmp,dk]=compute_spectral_flux_P(ctl,model,dx,kchoice,window,lims,lit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  HORIZONTAL PRESSURE WORK OF SPECTRAL EKE BUDGET: 
%  fft(u)* .  fft(grad P) 
%  Capet 2008
%  modif: Capet, Marchesiello 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flatbottom=1;
salinity=0;
g=9.81;rho0=1000; 

Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
L=Lmax-Lmin+1; M=Mmax-Mmin+1;

grd=rnt_gridload(model);
% horizontal grid
pm=grd.pm;pm=squeeze(pm(Lmin:Lmax,Mmin:Mmax));;
pn=grd.pn;pn=squeeze(pn(Lmin:Lmax,Mmin:Mmax));;
% Vertical grid
if flatbottom,
 hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;method=grd.Method;
 h=grd.h; h=h(1,1); %h=cut_var(h,lims);
 zr=zlevs(h,0,thetas,thetab,hc,N,'r');
 zw=zlevs(h,0,thetas,thetab,hc,N,'w');
 zr0=zr; zw0=zw;
 dz(2:N)=zw(3:N+1)-zw(2:N);
 dz(1)=dz(2);
 dz=dz(kchoice:N)';
else
 hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;method=grd.Method;
 h=grd.h; h=squeeze(h(Lmin:Lmax,Mmin:Mmax));
 h0=h(1,1)*ones(size(h));
 lon=grd.lonr; lon=squeeze(lon(Lmin:Lmax,Mmin:Mmax));
 lat=grd.latr; lat=squeeze(lat(Lmin:Lmax,Mmin:Mmax));
 zr=zlevs(h0,0,thetas,thetab,hc,N,'r');
 zw=zlevs(h0,0,thetas,thetab,hc,N,'w');
 zr=permute(zr,[2,3,1]);
 zw=permute(zw,[2,3,1]);
 dz=zeros(size(zr));
 dz(:,:,1:N-1)=zr(:,:,2:N)-zw(:,:,1:N-1);
 dz(:,:,N)=zw(:,:,N+1)-zr(:,:,N);
 dz=dz(:,:,kchoice:N);
 zr=zr(:,:,kchoice:N);
 zw=zw(:,:,kchoice:N+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% windowing: there are 2 options because all the signal is known
%            remove the trend signal(end)-signal(1)
%            apply a filter (rectangular 1, triangular 2 or hanning 3)
%
for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])

  u=rnt_loadvar_partialz(ctl,it,'u',kchoice,kchoice); % u is 3d
  u=perm(u2rho(perm(u)));
  u=squeeze(u(Lmin:Lmax,Mmin:Mmax));

  v=rnt_loadvar_partialz(ctl,it,'v',kchoice,kchoice); % u is 3d
  v=perm(v2rho(perm(v)));  % there was a bug there
  v=squeeze(v(Lmin:Lmax,Mmin:Mmax));

  temp=rnt_loadvar_partialz(ctl,it,'temp',kchoice,N);
  temp=squeeze(temp(Lmin:Lmax,Mmin:Mmax,:));
  if salinity
    salt=rnt_loadvar_partialz(ctl,it,'salt',kchoice,N);
    salt=squeeze(salt(Lmin:Lmax,Mmin:Mmax,:));
  else
    rho=rnt_loadvar_partialz(ctl,it,'rho',kchoice,N);
    rho=squeeze(rho(Lmin:Lmax,Mmin:Mmax,:));
  end
  zeta=rnt_loadvar(ctl,it,'zeta');
  zeta=squeeze(zeta(Lmin:Lmax,Mmin:Mmax));

  if flatbottom,
    if salinity,
      rho=rnt_rho_potential_nomex(temp,salt);
    end
    rho=rho+rho0;
    [L,M]=size(zeta);
    dz3d=permute(repmat(dz,[1,L,M]),[2,3,1]);
    dz3d(:,:,end)=squeeze(dz3d(:,:,end))+zeta;
    rhow=rho;
    rhow(:,:,1:end-1)=0.5*(rho(:,:,1:end-1)+rho(:,:,2:end));
    pressure=g.*sum(dz3d.*rhow,3);
    dxP=-partialx(pressure,pm)./rho0;  % watch out the - sign is included here
    dyP=-partialy(pressure,pn)./rho0;
  else
    zr=zlevs(h,zeta,thetas,thetab,hc,N,'r');
    zw=zlevs(h,zeta,thetas,thetab,hc,N,'w');
    zr=permute(zr,[2,3,1]);
    zw=permute(zw,[2,3,1]);
    zr=zr(:,:,kchoice:N);
    zw=zw(:,:,kchoice:N+1);
    rho=rho_eos(temp,salt,zr);
    [dxP,dyP]=prsgrd(rho,zw,zr,pm,pn);
    dxP=perm(u2rho(perm(squeeze(dxP(:,:,1)))));
    dyP=perm(v2rho(perm(squeeze(dyP(:,:,1)))));
  end
%  figure(9); pcolor(lon,lat,dxP);shading flat; colorbar;

  [u,v,dxP,dyP]=windowing(u,v,dxP,dyP,window);

  [M,L]=size(u);

  fftdim=2;
  if fftdim==1,
    if it==lit(1);
     tmpamp=zeros(L,1);
    end
    for i=1:M;
      fcoefu=fft(u(i,:),[],2)'; 
      fcoefv=fft(v(i,:),[],2)';
      fcoefpx=fft(dxP(i,:),[],2)';
      fcoefpy=fft(dyP(i,:),[],2)';
      tmpamp=tmpamp+real(conj(fcoefu).*fcoefpx+conj(fcoefv).*fcoefpy);
    end
  else
    if it==lit(1);
     tmpamp=zeros(M,L);
    end
    fcoefpx=fft2(dxP); fcoefpy=fft2(dyP);
    fcoefu=fft2(u); fcoefv=fft2(v);
    tmpamp=tmpamp+real(conj(fcoefu).*fcoefpx+conj(fcoefv).*fcoefpy);
  end

end % temporal loop

if fftdim==1,
  amp=1/(length(lit)*(L*M))*tmpamp;
  k=2*pi/dx*(1:L/2)/L;
  imin=1; imax=length(k);
  amp=0.05*amp(imin:imax)';
  k=k(imin:imax);
  dk=ones(size(k)).*(k(3)-k(2));
  count=zeros(size(k));
  ktmp=k;
else

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

end
return

