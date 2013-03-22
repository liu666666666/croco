function [k,E]=get_e_spectra(u2d,v2d,ut2d,vt2d,npts,dx,method);

[L,M]=size(u2d);

if method==1,
%
% Get 1D energy spectrum 
%
 demeaning=1;
 windowing=1;
 wdw=tukeywin(L,1);
 E=zeros(1,npts)';
 n=0;
 for j=1:M
  n=n+1;
  u=squeeze(u2d(:,j));
  v=squeeze(v2d(:,j));
  ut=squeeze(ut2d(:,j));
  vt=squeeze(vt2d(:,j));
  if demeaning,
   u=u-mean(u); 
   v=v-mean(v);
   ut=ut-mean(ut); 
   vt=vt-mean(vt);
  end
  if windowing,
   wdw=tukeywin(L,1);
   u=u.*wdw; 
   v=v.*wdw;
   ut=ut.*wdw; 
   vt=vt.*wdw;
  end
  FU=fft(u,npts);
  FV=fft(v,npts);
  FUT=fft(ut,npts);
  FVT=fft(vt,npts);
  E=real(FU.*conj(FUT) +  FV.*conj(FVT))./npts;
  k=2*pi*(0:npts/2)/(dx*npts);
  imin=max(find(k<2/(dx*(L-1))));
  imax=length(k);
  E=E(imin:imax);
  k=k(imin:imax);
 end
 E=E/n;
end
%
% 
%
return
