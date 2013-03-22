function [amp,count,ktmp,dk]=compute_power_spectrum_tracer(ctl,dx,kchoice,windo,lims,lit,name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [amp,count,ktmp]=compute_power_spectrum(ctl,dx,kchoice,windo,lims,it1,it2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  this function computes velocity power spectrum for a netcdf file series. 
%%%% inputs: 
%           ctl    : structure with netcdf files loaded
%           dx     : meshgrid size (!!! supposed to be uniform !!!)
%           kchoice: sigma level to be processed (z level not implemented)
%           windo : windoing option. There are 2 options. 
%                       apply a filter (rectangular 1, triangular 2 or hanning 3)
%                       apply a dct    (0) or do nothing (-1)
%           lims   : limits of the domain to be processed (4 elements vector with
%                       lims(1) and lims(2) mins and max in xsi direction. 
%           lit    : list of it times to process. used to be
%%%%%                            t1,t2  : time limits in the series. Default: 1 and length(ctl.time)
%           name of the tracer (consistent with netcdf file). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
pi=3.1415926535;

if nargin==5
%  it1=1; it2=length(ctl.time);
  lit=[1:length(ctl.time)];
end

for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])
%  disp(['u=rnt_loadvar_partialz(ctl,it,''' name ''',kchoice,kchoice);']); % u is 3d
  eval(['u=rnt_loadvar_partialz(ctl,it,''' name ''',kchoice,kchoice);']); % u is 3d
  u=squeeze(u(Lmin:Lmax,Mmin:Mmax));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% windowing: there are 2 options because all the signal is known
%            remove the trend signal(end)-signal(1)
%            apply a filter (rectangular 1, triangular 2 or hanning 3)
%            apply a dct filter transformation (0)
  [M,L]=size(u);
  if window==1
    wdw1=0.5-0.5*cos(2*pi*(1:L)/(L-1));
    wdw2=0.5-0.5*cos(2*pi*(1:M)/(M-1));
    wdw=wdw2'*wdw1;u=u.*wdw;
  elseif window==3
    u=mirror_dct(u);
    [M,L]=size(u);
  elseif window==4
    u=mirror_dct_1d(u);
    [M,L]=size(u);
  end
  if it==lit(1)
     tmpamp=zeros(M,L);
  end
  fcoefu=fft2(u);
  tmpamp=tmpamp+real(conj(fcoefu).*fcoefu);

end % temporal loop

% if check parseval's relationship
% lhswt=mean(mean(w.*T));
% rhswt=1/(Lu*Mu)^2*sum(sum(tmpampwT));

% do the reorganization and processing of results. 
tmpamp=1/(length(lit)*(L*M)^2)*tmpamp; % parseval equality
tmpamp=reorganize_fft2d(tmpamp);       % reorganize to be centered

% new version of reorganization
new=1;
if new==1
  method=2;
  [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,method);
else
  Lh=L/2;Mh=M/2;i2L=1/L^2;i2M=1/M^2;
  kk=repmat(perm([1:L])-Lh,[1 M]); ll=repmat([1:M]-Mh,[L 1]);
  Kh=i2L*kk.^2+i2M*ll.^2; Kh=6.28/(dx)*sqrt(Kh);Kh=perm(Kh);
  leng=L;% leng=round(sqrt(L^2+M^2))+1;
  %leng=round(0.5*sqrt(L^2+M^2));
  aux=leng/Kh(M,L);ktmp=[1:leng]/aux;ktmp=perm(ktmp);
  amp=zeros(leng,1);count=zeros(leng,1);
  ik = min(round(aux*Kh)+1,leng);
  for ind=1:leng
    II=find(ik==ind);
    amp(ind)=sum(tmpamp(II));
    count(ind)=length(II);
  end
%  amp=modif_amp(amp,count,ktmp,M,L,dx);
end
return
