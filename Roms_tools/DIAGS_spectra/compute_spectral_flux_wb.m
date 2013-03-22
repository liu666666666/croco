function [amp,count,ktmp,dk]=compute_spectral_flux_wb(ctl,model,dx,kchoice,window,lims,lit,ilevel);

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
g=9.81;rho0=1025;
%
submeso      = 0;
salinity     = 0;

grd=rnt_gridload(model);
pm=grd.pm;pm=squeeze(pm(Lmin:Lmax,Mmin:Mmax));
pn=grd.pn;pn=squeeze(pn(Lmin:Lmax,Mmin:Mmax));
%lon=grd.lonr;lon=squeeze(lon(Lmin:Lmax,Mmin:Mmax));
%lat=grd.latr;lat=squeeze(lat(Lmin:Lmax,Mmin:Mmax));
xr=grd.xr;xr=squeeze(xr(Lmin:Lmax,Mmin:Mmax));
yr=grd.yr;yr=squeeze(yr(Lmin:Lmax,Mmin:Mmax));
hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;h=grd.h;h0=h(1,1);method=grd.Method;
zr=zlevs(h,0,thetas,thetab,hc,N,'r');
zr=permute(zr,[2,3,1]);

if nargin==5
  lit=[1:length(ctl.time)];
end

cff_scale=0;

for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])
  w=rnt_loadvar_partialz(ctl,it,'w',kchoice,kchoice); % w is 2d
  if isempty(w)
    w=squeeze(mean(rnt_loadvar_partialz(ctl,it,'omega',kchoice,kchoice+1),3)); % w is 3d
  end
  w=squeeze(w(Lmin:Lmax,Mmin:Mmax));
  T=rnt_loadvar_partialz(ctl,it,'temp',kchoice,kchoice); % temp is 2d
  if salinity,
    S=rnt_loadvar_partialz(ctl,it,'salt',kchoice,kchoice); % salt is 2d
    %rho=rho_eos(T,S,zr(:,:,kchoice));
    rho=rnt_rho_potential_nomex(T,S);
    rho=squeeze(rho(Lmin:Lmax,Mmin:Mmax));
  else
    rho=rnt_loadvar_partialz(ctl,it,'rho',kchoice,kchoice); % rho is 2d
    rho=squeeze(rho(Lmin:Lmax,Mmin:Mmax));
  end
  T=-g/rho0*rho;

  if submeso==1,
    nhan =  10.^(ilevel-1); % nb of hanning filtering
    navg =  1; % avg period = navg*2+1 days
    tmp=low_filter1(ctl,'rho',it,kchoice,nhan,navg);
    rfilt=tmp(Lmin:Lmax,Mmin:Mmax);
    T=-g/rho0*(rho-rfilt);
    %T=-g/rho0*rfilt;
    tmp=low_filter1(ctl,'w',it,kchoice,nhan,navg);
    wfilt=tmp(Lmin:Lmax,Mmin:Mmax);
    w=w-wfilt;
    %w=wfilt;
  end

  %figure; pcolor(lon,lat,w.*T); shading flat;colorbar;
  %stop

  [w,w,T,T]=windowing(w,w,T,T,window);
  cff_scale=cff_scale+mean(mean(w.*T));

  [M,L]=size(w);

  fftdim=2;
  if fftdim==1,
    if it==lit(1);
     tmpamp=zeros(L,1);
    end
    for i=1:M;
      fcoefT=fft(T(i,:),[],2)'; fcoefw=fft(w(i,:),[],2)';
      tmpamp=tmpamp+real(conj(fcoefT).*fcoefw);
    end
  else
    if it==lit(1);
     tmpamp=zeros(M,L);
    end
    fcoefT=fft2(T); fcoefw=fft2(w);
    tmpamp=tmpamp+real(conj(fcoefT).*fcoefw);
  end

end % temporal loop

cff_scale=1/(length(lit))*cff_scale;

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

% Correct for true integrated value
amp=amp*cff_scale/sum(amp);

end
return
