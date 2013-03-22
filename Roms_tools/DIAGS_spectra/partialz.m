function [dzvar]=partialz(var,pz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%   compute the relative z derivative of var, pz is the metrics   %
%   array of size Kv (pzr or pzw depending on whether var is      % 
%   defined at w or r levels.                                     %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Mv,Lv,Kv]=size(var);
dzvar=zeros(Mv,Lv,Kv);
pz=perm(repmat(pz,[1,Lv,Mv]));
dzvar(:,:,2:Kv-1)=0.5.*pz(:,:,2:Kv-1).*(var(:,:,3:Kv)-var(:,:,1:Kv-2));

dzvar(:,:,1)=dzvar(:,:,2);
dzvar(:,:,Kv)=dzvar(:,:,Kv-1);
return
