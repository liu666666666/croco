function [dzvar]=partialz_r2w(var,zr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   compute the relative z derivative of var (available at rho levels      %
%   zr is the rho level height. The result will be provided at w levels    %
%   with no empty spot ie the output has one less point                    %
%   THIS ONE IS ACTUALLY IDENTICAL TO partialz_w2r : it just needs to be   % 
%   called with zr instead of zw. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Mv,Lv,Kv]=size(var);
dzvar=zeros(Mv,Lv,Kv-1);

zr=perm(repmat(zr,[1,Lv,Mv]));
%for ii=1:Mv
%  for jj=1:Lv
%    for kk=2:Kv
%      dzvar(ii,jj,kk)=(zr(kk)-zw(kk))*var(ii,jj,kk-1)+(zw(kk)-zr(kk-1))*var(ii,jj,kk);
%      dzvar(ii,jj,kk)=dzvar(ii,jj,kk)./(zr(kk)-zr(kk-1))
%    end
%  end
%end
%dzvar(:,:,2:Kv)=(zr(:,:,2:Kv)-zw(:,:,2:Kv)).*var(:,:,1:Kv-1)+(zw(:,:,2:Kv)-zr(:,:,1:Kv-1)).*var(:,:,2:Kv);
%dzvar(:,:,2:Kv)=dzvar(:,:,2:Kv)./(zr(:,:,2:Kv)-zr(:,:,1:Kv-1));
dzvar(:,:,1:Kv-1)=(var(:,:,2:Kv)-var(:,:,1:Kv-1))./(zr(:,:,2:Kv)-zr(:,:,1:Kv-1));
return
