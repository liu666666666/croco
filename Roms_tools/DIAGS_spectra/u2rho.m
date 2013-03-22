function [var_rho]=u2rho(var_u,transp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [var_rho]=u2rho(var_u,transp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% transfert a 2d or 3d field at u points to a field at rho points
% transp=1 works on fields loaded with manu's convention
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==2
  if transp==1
     var_u=perm(var_u);
  end
end
if length(size(var_u))==2
  [Mp,L]=size(var_u);
  Lp=L+1;
  Lm=L-1;
  M=Mp-1;
  var_rho=zeros(Mp,Lp);

  var_rho(:,2:L)=0.5*(var_u(:,1:Lm)+var_u(:,2:L));

  %var_rho(:,2:L)=var_u(:,2:L);
  %var_rho(2:Mp,1:L)=0.5*(var_u(1:M,:)+var_u(2:Mp,:));
  %var_rho(1,:)=var_rho(2,:);

  var_rho(:,1)=var_rho(:,2);
  var_rho(:,Lp)=var_rho(:,L);

elseif length(size(var_u))==3
  [N,Mp,L]=size(var_u);
  Lp=L+1;
  Lm=L-1;
  M=Mp-1;
  var_rho=zeros(N,Mp,Lp);

  var_rho(:,:,2:L)=0.5*(var_u(:,:,1:Lm)+var_u(:,:,2:L));

  %var_rho(:,:,2:L)=var_u(:,:,2:L);
  %var_rho(:,2:Mp,1:L)=0.5*(var_u(:,1:M,:)+var_u(:,2:Mp,:));
  %var_rho(:,1,:)=var_rho(:,2,:);

  var_rho(:,:,1)=var_rho(:,:,2);
  var_rho(:,:,Lp)=var_rho(:,:,L);

end
if nargin==2
  if transp==1
     var_rho=perm(var_rho);
  end
end
return

