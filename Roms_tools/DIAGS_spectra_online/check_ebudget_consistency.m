%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  this program checks the consistence of diagnostics terms 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%addpath('../DIAGS_spectra')

window=1;
root='./ROMS_FILES/';
pi=3.1415926535;dx=36000;

imin=1;imax=106; jmin=1;jmax=62;  % case where full domain is used 
imin=5;imax=102; jmin=5;jmax=58; % case where part domain is used

removem=1;                        % removes average or not. 

iminu=imin;imaxu=imax-1;jminu=jmin;jmaxu=jmax;
iminv=imin;imaxv=imax;jminv=jmin;jmaxv=jmax-1;

diaguv=[root 'roms_diaM_avg.0000.nc'];ncdiag=netcdf(diaguv,'r');
hisuv=[root 'roms_his.0000.nc'];nchis=netcdf(hisuv,'r');
it=80;kchoice=28; it2=81;

res=-ncdiag{'u_rate'}(it,kchoice,:,:)+ncdiag{'u_xadv'}(it,kchoice,:,:)+ncdiag{'u_yadv'}(it,kchoice,:,:);
res=res+ncdiag{'u_vadv'}(it,kchoice,:,:)+ncdiag{'u_cor'}(it,kchoice,:,:)+ncdiag{'u_Prsgrd'}(it,kchoice,:,:);
res=res+ncdiag{'u_vmix'}(it,kchoice,:,:)+ncdiag{'u_hmix'}(it,kchoice,:,:);

u=nchis{'u'}(it2,kchoice,jminu:jmaxu,iminu:imaxu);
v=nchis{'v'}(it2,kchoice,jminv:jmaxv,iminv:imaxv);

figure;imagesc(res);colorbar;caxis([-1e-9 1e-9]);title(['residual in energy online energy budget']);colorbar

advu=-ncdiag{'u_xadv'}(it,kchoice,jminu:jmaxu,iminu:imaxu)-ncdiag{'u_yadv'}(it,kchoice,jminu:jmaxu,iminu:imaxu)...
  -ncdiag{'u_vadv'}(it,kchoice,jminu:jmaxu,iminu:imaxu);
advv=-ncdiag{'v_xadv'}(it,kchoice,jminv:jmaxv,iminv:imaxv)-ncdiag{'v_yadv'}(it,kchoice,jminv:jmaxv,iminv:imaxv)...
  -ncdiag{'v_vadv'}(it,kchoice,jminv:jmaxv,iminv:imaxv);
mixu=-ncdiag{'u_hmix'}(it,kchoice,jminu:jmaxu,iminu:imaxu)-ncdiag{'u_vmix'}(it,kchoice,jminu:jmaxu,iminu:imaxu);
mixv=-ncdiag{'v_hmix'}(it,kchoice,jminv:jmaxv,iminv:imaxv)-ncdiag{'v_vmix'}(it,kchoice,jminv:jmaxv,iminv:imaxv);
pcfu=-ncdiag{'u_cor'}(it,kchoice,jminu:jmaxu,iminu:imaxu)-ncdiag{'u_Prsgrd'}(it,kchoice,jminu:jmaxu,iminu:imaxu);
pcfv=-ncdiag{'v_cor'}(it,kchoice,jminv:jmaxv,iminv:imaxv)-ncdiag{'v_Prsgrd'}(it,kchoice,jminv:jmaxv,iminv:imaxv);
ratu=ncdiag{'u_rate'}(it,kchoice,jminu:jmaxu,iminu:imaxu);
ratv=ncdiag{'v_rate'}(it,kchoice,jminv:jmaxv,iminv:imaxv);

if removem==1
  advu=advu-mean(mean(advu));
  advv=advv-mean(mean(advv));
  mixu=mixu-mean(mean(mixu));
  mixv=mixv-mean(mean(mixv));
  pcfu=pcfu-mean(mean(pcfu));
  pcfv=pcfv-mean(mean(pcfv));
  ratu=ratu-mean(mean(ratu));
  ratv=ratv-mean(mean(ratv));
end
if window==1 
  [Mu,Lu]=size(u);
  M=Mu; L=Lu;
  wdw1=0.5-0.5*cos(2*pi*(0:L-1)/(L-1));
  wdw2=0.5-0.5*cos(2*pi*(0:M-1)/(M-1));
  wdw=wdw2'*wdw1;
  u=u.*wdw;       
  advu=advu.*wdw; 
  mixu=mixu.*wdw; 
  pcfu=pcfu.*wdw; 
  ratu=ratu.*wdw; 

  [Mv,Lv]=size(v);
  M=Mv; L=Lv;
  wdw1=0.5-0.5*cos(2*pi*(0:L-1)/(L-1));
  wdw2=0.5-0.5*cos(2*pi*(0:M-1)/(M-1));
  wdw=wdw2'*wdw1;
  v=v.*wdw;
  advv=advv.*wdw;
  mixv=mixv.*wdw;
  pcfv=pcfv.*wdw;
  ratv=ratv.*wdw;

elseif window==2 % windowing for tranposed es that has periodicity in eta
  [M,L]=size(u);
  wdw1=ones(1,L);
  wdw2=0.5-0.5*cos(2*pi*(0:M-1)/(M-1));
  wdw=wdw2'*wdw1;u=u.*wdw; v=v.*wdw;
  advu=advu.*wdw; advv=advv.*wdw;
end

option=1;
if option==1
  M=Mu;L=Lv;
elseif option==2
  M=Mu;L=Lu;
elseif option==3
  M=Mv;L=Lv;
end

tmpadv=zeros(M,L); 
tmpmix=zeros(M,L); 
tmppcf=zeros(M,L); 
tmprat=zeros(M,L); 

if option==1
  fcoefadvu=fft2(u2rho(advu)); fcoefadvv=fft2(v2rho(advv));
  fcoefmixu=fft2(u2rho(mixu)); fcoefmixv=fft2(v2rho(mixv));
  fcoefpcfu=fft2(u2rho(pcfu)); fcoefpcfv=fft2(v2rho(pcfv));
  fcoefratu=fft2(u2rho(ratu)); fcoefratv=fft2(v2rho(ratv));
  fcoefu=fft2(u2rho(u)); fcoefv=fft2(v2rho(v));

  tmpadv=tmpadv+real(conj(fcoefu).*fcoefadvu+conj(fcoefv).*fcoefadvv);
  tmpmix=tmpmix+real(conj(fcoefu).*fcoefmixu+conj(fcoefv).*fcoefmixv);
  tmppcf=tmppcf+real(conj(fcoefu).*fcoefpcfu+conj(fcoefv).*fcoefpcfv);
  tmprat=tmprat+real(conj(fcoefu).*fcoefratu+conj(fcoefv).*fcoefratv);

elseif option==2

  fcoefadvu=fft2((advu)); 
  fcoefmixu=fft2((mixu)); 
  fcoefpcfu=fft2((pcfu)); 
  fcoefratu=fft2((ratu)); 
  fcoefu=fft2((u)); fcoefv=fft2(v);

  tmpadv=tmpadv+real(conj(fcoefu).*fcoefadvu);
  tmpmix=tmpmix+real(conj(fcoefu).*fcoefmixu);
  tmppcf=tmppcf+real(conj(fcoefu).*fcoefpcfu);
  tmprat=tmprat+real(conj(fcoefu).*fcoefratu);

elseif option==3

  fcoefadvv=fft2(advv);
  fcoefmixv=fft2(mixv);
  fcoefpcfv=fft2(pcfv);
  fcoefratv=fft2(ratv);
  fcoefv=fft2(v);

  tmpadv=tmpadv+real(conj(fcoefv).*fcoefadvv);
  tmpmix=tmpmix+real(conj(fcoefv).*fcoefmixv);
  tmppcf=tmppcf+real(conj(fcoefv).*fcoefpcfv);
  tmprat=tmprat+real(conj(fcoefv).*fcoefratv);

end

tmpadv=1/(L*M)^2*tmpadv; % parseval equality
tmpmix=1/(L*M)^2*tmpmix; % parseval equality
tmppcf=1/(L*M)^2*tmppcf; % parseval equality
tmprat=1/(L*M)^2*tmprat; % parseval equality

tmpadv=reorganize_fft2d(tmpadv);       % reorganize to be centered
tmpmix=reorganize_fft2d(tmpmix);       % reorganize to be centered
tmppcf=reorganize_fft2d(tmppcf);       % reorganize to be centered
tmprat=reorganize_fft2d(tmprat);       % reorganize to be centered

method=2;
[ampadv,count,ktmp,dk]=integ_fft2d(tmpadv,dx,M,L,method);
[ampmix,count,ktmp,dk]=integ_fft2d(tmpmix,dx,M,L,method);
[amppcf,count,ktmp,dk]=integ_fft2d(tmppcf,dx,M,L,method);
[amprat,count,ktmp,dk]=integ_fft2d(tmprat,dx,M,L,method);

windo=0;ilevel=1;

figure;
plot(ktmp,ampadv,'b');
hold on;
plot(ktmp,ampadv+ampmix+amppcf+amprat,'r');
set(gca,'Yscale','linear','Xscale','log');



