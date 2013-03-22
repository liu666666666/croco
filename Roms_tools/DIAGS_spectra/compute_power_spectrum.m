function [amp,count,ktmp,dk]=compute_power_spectrum(ctl,dx,kchoice,window,npts,lims,lit,grd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [amp,count,ktmp]=compute_power_spectrum(ctl,dx,kchoice,window,lims,lit,grd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%           lit    : time vector in the series. Default: [1:length(ctl.time)]
%           grd    : grd structure so that geostrophic velocities can be computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
pi=3.1415926535;

if nargin<7
  lit=[1:length(ctl.time)];
end
if nargin<8 
   var=1;
else
   var=2;
end

for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])
  if var==1
    u=rnt_loadvar_partialz(ctl,it,'u',kchoice,kchoice); % u is 3d
    v=rnt_loadvar_partialz(ctl,it,'v',kchoice,kchoice); % v is 3d
  else
    [u,v]=geostrophic_uv(ctl,grd,it,kchoice);
  end
  u=perm(u2rho(perm(u)));
  u=squeeze(u(Lmin:Lmax,Mmin:Mmax));
  v=perm(v2rho(perm(v)));
  v=squeeze(v(Lmin:Lmax,Mmin:Mmax));
  u=u-mean(mean(u));
  v=v-mean(mean(v));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% windowing: there are 2 options because all the signal is known
%            remove the trend signal(end)-signal(1)
%            apply a filter (rectangular 1, triangular 2 or hanning 3)
%            apply a dct filter transformation (0)
  [M,L]=size(u);
  if window==1
    wdw1=0.5-0.5*cos(2*pi*(1:L)/(L-1));
    wdw2=0.5-0.5*cos(2*pi*(1:M)/(M-1));
    wdw=wdw2'*wdw1;u=u.*wdw; v=v.*wdw;
  elseif window==2
    %wdw1=gausswin(L,3)'; wdw2=gausswin(M,3)';
    wdw1=tukeywin(L,0.25)'; wdw2=tukeywin(M,0.25)';
    %wdw1=blackman(L,'periodic')'; wdw2=blackman(M,'periodic')';
    %wdw1=hann(L,'periodic')'; wdw2=hann(M,'periodic')';
    %wdw1=barthannwin(L)'; wdw2=barthannwin(M)';
    %wdw1=blackmanharris(L)'; wdw2=blackmanharris(M)';
    wdw=wdw2'*wdw1;u=u.*wdw; v=v.*wdw;
  elseif window==3
    u=mirror_dct(u);v=mirror_dct(v);
  elseif window==4
    u=mirror_dct_1d(u);v=mirror_dct_1d(v);
  end

  [L,M]=size(u);
  if npts>0 & npts>L & npts>M,
   M=npts; L=npts;
  end

  if it==lit(1);
     tmpamp=zeros(M,L);
  end

  fcoefu=fft2(u,npts,npts); fcoefv=fft2(v,npts,npts);
  tmpamp=0.5*(tmpamp+real(conj(fcoefu).*fcoefu+conj(fcoefv).*fcoefv));

end % temporal loop

%-----------------------------------------------------------------
% if check parseval's relationship
% lhswt=mean(mean(w.*T));
% rhswt=1/(Lu*Mu)^2*sum(sum(tmpampwT));

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
%  leng=round(0.5*sqrt(L^2+M^2));
  aux=leng/Kh(M,L);ktmp=[1:leng]/aux;ktmp=perm(ktmp);
  amp=zeros(leng,1);count=zeros(leng,1);
  dk=ktmp(10)-ktmp(9);
  ik = min(round(aux*Kh)+1,leng);
  for ind=1:leng
    II=find(ik==ind);
    amp(ind)=sum(tmpamp(II));
    count(ind)=length(II);
  end
end
return

