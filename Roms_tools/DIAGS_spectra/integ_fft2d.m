function [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,method);

%M=length(xi); L=length(eta);

Lh=L/2;Mh=M/2;
cff=2*pi/dx;
i2L=1/L^2;i2M=1/M^2;
kk=repmat(perm([1:L])-Lh,[1 M]); ll=repmat([1:M]-Mh,[L 1]);
Kh=i2L*kk.^2+i2M*ll.^2; Kh=6.28/(dx)*sqrt(Kh);Kh=perm(Kh);

%if M==L;method=1;end

if method==1 % old case
  %
  % there are several issues with this old part. 
  % most importantly the first non zero wavenumber corresponds to 
  % a smaller value than any computed by fft (Kh(M,L)/L)
  %
  leng=L;
  aux=leng/Kh(M,L);ktmp=[1:leng]/aux;ktmp=perm(ktmp);
  amp=zeros(leng,1);count=zeros(leng,1);
  
  ik = min(round(aux*Kh)+1,leng);
  for ind=1:leng
    II=find(ik==ind);
    amp(ind)=sum(tmpamp(II));
    count(ind)=length(II);
  end
  dk=ktmp(10)-ktmp(9);
  amp=modif_amp(amp,count,ktmp,dk);
 
elseif method==2 % new case
  %
  % Ring integration method with non constant dk
  % accounts for different lengths in xi and eta
  %
  if M>L
    dcase=1;Lmin=L;Lmax=M;% %% xi dim is bigger
  else
    dcase=0;Lmin=M;Lmax=L;
  end
%  zob=ones(size(tmpamp)); % for debug purposes. 
  %%% first sample the long wavelengths permitted in longest 
  %%% direction with wavenumber 0 in shortest direction
  indaux=0;
  for ind=1:floor(Lmax/Lmin)-1
    amp(ind)=tmpamp(Mh+dcase*ind,Lh+(1-dcase)*ind);
    count(ind)=1;
    ktmp(ind)=cff/Lmax*ind;
    dk(ind)=cff/Lmax;
%    zob(Mh+dcase*ind,Lh+(1-dcase)*ind)=zob(Mh+dcase*ind,Lh+(1-dcase)*ind)-1;
    indaux=indaux+1;
  end
  %%% second do the regular summation over rings with 
  %%% step that is function of k step for shortest dimension 
  %%% ie we use the largest step.
%  Khtmp=Kh;
%  Khtmp(II)=nan;
%  pcolor(Khtmp);shading flat;
  for ind=floor(Lmax/Lmin):Lmin/2
    ktmp(ind)=cff/Lmin*(ind-indaux);
    II=find(abs(Kh-ktmp(ind))<cff/(2*Lmin));
    [xI,yI]=find(abs(Kh-ktmp(ind))<cff/(2*Lmin));
    amp(ind)=sum(tmpamp(II));
    count(ind)=length(II);
    dk(ind)=cff/Lmin;
%    zob(II)=zob(II)-1;
  end
%  figure;imagesc(zob);shading flat;title(['should be nan']);
  %%% third correct for spectral defects due to summation over 
  %%% finite number of points that imperfectly sample the wavenumber 
  %%% space. 
  amp=modif_amp(amp,count,ktmp,dk);
end
return
