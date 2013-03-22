function Z=filter_pp_1D(x,z,R);
%
% get the anomalies between a field and a smoothed field
%

pio2=pi/2;
cff1=pio2/R;
R2=R^2;
L=length(z);
dx=mean(x(2:L)-x(1:L-1));
%
% build the stencil
%
x1=(dx:dx:R);
x1=[fliplr(-x1) 0 x1];
[ni]=length(x1);
%
% get the filter coefficients
%
R0=x1.^2;
cosX=0.*R0;
D=find(R0<R2);
cosX(D)=cos(cff1*sqrt(R0(D)));
%
% extend the domain
%
bi=(ni-1)/2;
zbig=zeros(L+ni-1);
L=length(zbig);
zbig(1+bi:L-bi)=z;
%
% get the mask
%
warning off
mask=zbig./zbig;
warning on
mask(isnan(mask))=0;
%
% apply the filter
%
Z=0.*z;
cff=Z;
[i]=find(cosX~=0);
npts=length(i);
for n=1:npts
  Z=Z+cosX(i(n))*zbig(i(n):L-ni+i(n));
  cff=cff+cosX(i(n))*zbig(i(n):L-ni+i(n));
end
Z=Z./cff;
Z=Z-z;


