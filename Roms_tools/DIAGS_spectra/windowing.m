function [u,v,termu,termv]=windowing(u,v,termu,termv,window)

%======================================================================
% windowing: there are 2 options because all the signal is known
%            remove the trend signal(end)-signal(1)
%            apply a filter (hanning or others in matlab toolbox)
%            apply a dct filter transformation (0)
%======================================================================

mmean=1;
if mmean,
 u=u-mean(mean(u));
 v=v-mean(mean(v));
 termu=termu-mean(mean(termu));
 termv=termv-mean(mean(termv));
end

[M,L]=size(u);

%figure; pcolor(termu); colorbar; shading flat; caxis([-5.e-6 5.e-6])

if window==0            % no windowing!
      %%
elseif window==1        % Hann (or Hanning) window 
  wdw1=0.5-0.5*cos(2*pi*(1:L)/(L-1));
  wdw2=0.5-0.5*cos(2*pi*(1:M)/(M-1));
  wdw=wdw2'*wdw1;
  u=u.*wdw;
  v=v.*wdw;
  termu=termu.*wdw;
  termv=termv.*wdw;
elseif window==2        % other matlab windowing methods:
  filt_cff=0.5;           % 1=Hann; 0=rectangle; 0.5=default
  wdw1=tukeywin(L,filt_cff)'; wdw2=tukeywin(M,filt_cff)';
  wdw=wdw2'*wdw1;
  u=u.*wdw;
  v=v.*wdw;
  termu=termu.*wdw;
  termv=termv.*wdw;
elseif window==3
  u=mirror_dct(u);
  v=mirror_dct(v);
  termu=mirror_dct(termu);
  termv=mirror_dct(termv);
  %[M,L]=size(u);
  %wdw1=tukeywin(L,1)'; wdw2=tukeywin(M,1)';
  %wdw=wdw2'*wdw1;
  %u=u.*wdw;
  %v=v.*wdw;
  %termu=termu.*wdw;
  %termv=termv.*wdw;
elseif window==4
  u=mirror_dct_1d(u);
  v=mirror_dct_1d(v);
  termu=mirror_dct_1d(termu);
  termv=mirror_dct_1d(termv);
end

%figure; pcolor(termu); colorbar; shading flat; caxis([-5.e-6 5.e-6])
%figure; pcolor(wdw); colorbar; shading flat;
%stop

return
end

