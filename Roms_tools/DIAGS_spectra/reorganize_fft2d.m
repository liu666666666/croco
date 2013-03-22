function out=reorganize_fft2d(inp)
% this function takes a fft2d as computed from 
% matlab (constant terms and the lower and left borders) 
% and reorganizes its quadrants so as to have the coefficient 
% corresponding to 0 wavenumber right in the middle. 

[M,L]=size(inp);
Mh=M/2;Lh=L/2;
if (Mh~=round(Mh) | Mh~=round(Mh)) 
  disp(['size of input array not divisible by 2 ????'])
  stop
end
out=zeros(M,L);
out(Mh:M,Lh:L)=inp(1:Mh+1,1:Lh+1);    % lower left quadrant into upper right
out(Mh:M,1:Lh-1)=inp(1:Mh+1,Lh+2:L);  % lower right quadrant into upper left
out(1:Mh-1,Lh:L)=inp(Mh+2:M,1:Lh+1);  % upper left quadrant into lower right
out(1:Mh-1,1:Lh-1)=inp(Mh+2:M,Lh+2:L);% upper right quadrant into lower left
return
