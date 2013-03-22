function [pt]=plot_spectrum(ktmp,spec,N0,llim,ln);

%function fig=plot_spectrum(ktmp,spec,figin,N0,llim);
%  THIS FUNCTION PLOTS SPECTRUM EG COMING OUT OF SPECSMOZ

pro={'b';'r';'g';'c';'m';'y';'k'};
%pro={'r';'b';'m';'k';'k';'k';'k'};

if nargin==5
  prop=pro(ln);
end
if nargin<4 | isnan(llim)
  llim=length(ktmp); 
end
if nargin<3 | isnan(N0)
mid=floor(llim/8);
else
mid=max(min(floor(llim/N0),llim),1);
end

pt=plot(ktmp(1:llim),spec(1:llim),pro{ln}); hold on;
set(pt,'Linewidth',2)
set(gca,'Yscale','log','Xscale','log');

if ln==2
  cff0=1;
  spec1=ceil(9/10*mid);  %ceil(8/10*mid)
  spec2=floor(12/10*mid);
  cff0=cff0.*max(spec(spec1:spec2));
%
  %cff=cff0./(ktmp(mid).^(-5/3));
  %hh=plot(ktmp(1:llim),cff.*ktmp(1:llim).^(-5/3),'k:');
  %set(hh,'Linewidth',2);
  %text(ktmp(end),cff*ktmp(end).^(-5/3),'k^{-5/3}')
%
  cff=cff0./(ktmp(mid).^(-2));
  hi=plot(ktmp(1:llim),cff.*ktmp(1:llim).^(-2),'k--');
  set(hi,'Linewidth',1);
  text(ktmp(end),cff*ktmp(end).^(-2),'k^{-2}')
%
  cff=cff0./(ktmp(mid).^(-3));
  hj=plot(ktmp(1:llim),cff.*ktmp(1:llim).^(-3),'k-.'); 
  set(hj,'Linewidth',1);
  text(ktmp(end),cff*ktmp(end).^(-3),'k^{-3}')
end
