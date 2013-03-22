function [dzvar]=partialz_w2r(var,zw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%   compute the relative z derivative of var, zw is the  w level height    %
%   array of size [Mp,Lp,1]. This is done for a 2d or 3d variable          %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Mv,Lv,Kv]=size(var);
dzvar=zeros(Mv,Lv,Kv-1);

if length(size(zw))<3
  zw=perm(repmat(zw,[1,Lv,Mv]));
end
  dzvar(:,:,1:Kv-1)=(var(:,:,2:Kv)-var(:,:,1:Kv-1))./(zw(:,:,2:Kv)-zw(:,:,1:Kv-1));
return
