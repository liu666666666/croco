function Z=filter_pp(x,y,z,R,flag);
%
% filter a field z at location (x,y) [meters]
% or get anomalies (if flag=0).
% R is the decorrelation scale [meters]

if nargin<5, flag=1; end;

pio2=pi/2;
cff1=pio2/R;
R2=R^2;
[M,L]=size(z);
dx=mean(mean((x(:,2:L)-x(:,1:L-1))));
dy=mean(mean((y(2:M,:)-y(1:M-1,:))));
%
% build the stencil
%
x1=(dx:dx:R);
y1=(dy:dy:R);
x1=[fliplr(-x1) 0 x1];
y1=[fliplr(-y1) 0 y1];
[X1,Y1]=meshgrid(x1,y1);
[nj,ni]=size(X1);
%
% get the filter coefficients
%
R0=X1.^2+Y1.^2;
cosX=0.*R0;
D=find(R0<R2);
cosX(D)=cos(cff1*sqrt(R0(D)));
%
% extend the domain
%
bj=(nj-1)/2;
bi=(ni-1)/2;
zbig=zeros(M+nj-1,L+ni-1);
[M,L]=size(zbig);
zbig(1+bj:M-bj,1+bi:L-bi)=z;
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
[j,i]=find(cosX~=0);
npts=length(i);
for n=1:npts
  Z=Z+cosX(j(n),i(n))*zbig(j(n):M-nj+j(n),i(n):L-ni+i(n));
  cff=cff+cosX(j(n),i(n))*mask(j(n):M-nj+j(n),i(n):L-ni+i(n));
end
Z=Z./cff;

if flag==0,
  Z=Z-z;
end;



