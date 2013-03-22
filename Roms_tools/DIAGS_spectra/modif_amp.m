function amp=modif_amp(amp,count,ktmp,dk);


% function amp=modif_amp(amp,count,ktmp,dk);
% adjusts the amp values to correct for the 
% non exactness of discrete versus continuous 
% ring integration
if nargin<=3
dk=ktmp(5)-ktmp(4);
end
cff=2*pi*dk(1:end).*ktmp(1:end)/(pi*ktmp(end).^2).*sum(count)./count(1:end);
amp=amp.*cff;
%tmp=pi/dx;
%amp(ktmp>tmp)=nan;
return


