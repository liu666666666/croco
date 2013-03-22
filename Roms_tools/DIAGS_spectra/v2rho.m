function [var_rho]=v2rho(var_v,transp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [var_rho]=v2rho(var_v,transp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% transfert a field at v points to a field at rho points
% transp=1 works on fields loaded with manu's convention
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==2
  if transp==1
     var_v=perm(var_v);
  end
end

if length(size(var_v))==2
  [M,Lp]=size(var_v);
  Mp=M+1;
  Mm=M-1;
  L=Lp-1;
  var_rho=zeros(Mp,Lp);

  var_rho(2:M,:)=0.5*(var_v(1:Mm,:)+var_v(2:M,:));

  %var_rho(2:M,:)=var_v(2:M,:);
  %var_rho(1:M,2:Lp)=0.5*(var_v(:,1:L)+var_v(:,2:Lp));
  %var_rho(:,1)=var_rho(:,2);
   
  var_rho(1,:)=var_rho(2,:);
  var_rho(Mp,:)=var_rho(M,:);

elseif length(size(var_v))==3
  [N,M,Lp]=size(var_v);
  Mp=M+1;
  Mm=M-1;
  L=Lp-1;
  var_rho=zeros(N,Mp,Lp);

  var_rho(:,2:M,:)=0.5*(var_v(:,1:Mm,:)+var_v(:,2:M,:));

  %var_rho(:,2:M,:)=var_v(:,2:M,:);
  %var_rho(:,1:M,2:Lp)=0.5*(var_v(:,:,1:L)+var_v(:,:,2:Lp));
  %var_rho(:,:,1)=var_rho(:,:,2);

  var_rho(:,1,:)=var_rho(:,2,:);
  var_rho(:,Mp,:)=var_rho(:,M,:);

end
if nargin==2
  if transp==1
     var_rho=perm(var_rho);
  end
end

return
