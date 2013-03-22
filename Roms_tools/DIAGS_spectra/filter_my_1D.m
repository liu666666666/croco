function [Z]=filter_my_1D(x,z,R);

warning off;

if R==0, Z=z; return; end;

pio2=pi/2;
cff1=pio2/R;

L=length(x);

for i=1:L; 

%  disp([i,L])

  X=zeros(L,1);
  x0=x-x(i);
  D=find(x0.^2<R^2 & z~=0);
  X(D)=cff1*sqrt(x0(D).^2);
  cff=1/sum(cos(X(D)));
  Z(i)=cff*sum(z(D).*cos(X(D)));

end;

warning on;

return


