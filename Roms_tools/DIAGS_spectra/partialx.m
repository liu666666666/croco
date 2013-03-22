function [dxvar]=partialx(var,pm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%   compute the relative x derivative of var, pm is the metrics   %
%   array of size [Mp,Lp,1]. This is done for a 2d or 3d variable %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Mp,Lp]=size(pm);
siz=size(var);
if length(siz)==2
[Mv,Lv]=size(var);
dxvar=zeros(Mv,Lv);
dxvar(2:Mv-1,:)=0.5.*pm(2:Mv-1,:).*(var(3:Mv,:)-var(1:Mv-2,:));
dxvar(1,:)=dxvar(2,:);
dxvar(Mv,:)=dxvar(Mv-1,:);
elseif length(siz)==3  % 3d variable
[Mv,Lv,Kv]=size(var);
dxvar=zeros(Mv,Lv,Kv);
pm=repmat(pm,[1 1 Kv]);
dxvar(2:Mv-1,:,:)=0.5.*pm(2:Mv-1,:,:).*(var(3:Mv,:,:)-var(1:Mv-2,:,:));
dxvar(1,:,:)=dxvar(2,:,:);
dxvar(Mv,:,:)=dxvar(Mv-1,:,:);
elseif length(siz)==4 % time variable
[Mv,Lv,Kv,Tv]=size(var);
dxvar=zeros(Mv,Lv,Kv,Tv);
pm=repmat(pm,[1 1 Kv Tv]);
dxvar(2:Mv-1,:,:,:)=0.5.*pm(2:Mv-1,:,:,:).*(var(3:Mv,:,:,:)-var(1:Mv-2,:,:,:));
dxvar(1,:,:,:)=dxvar(2,:,:,:);
dxvar(Mv,:,:,:)=dxvar(Mv-1,:,:,:);
end
return
