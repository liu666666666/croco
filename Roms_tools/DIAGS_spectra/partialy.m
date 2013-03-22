function [dyvar]=partialy(var,pn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%   compute the y derivative of var, pm is the metrics            %
%   array of size [Mp,Lp,1]. This is done for a 2d or 3d variable %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
siz=size(var);
if length(siz)==2
[Mv,Lv]=size(var);
dyvar=zeros(Mv,Lv);
dyvar(:,2:Lv-1)=0.5.*pn(:,2:Lv-1).*(var(:,3:Lv)-var(:,1:Lv-2));
dyvar(:,1)=dyvar(:,2);
dyvar(:,Lv)=dyvar(:,Lv-1);
elseif length(siz)==3  % 3d variable
[Mv,Lv,Kv]=size(var);
dyvar=zeros(Mv,Lv,Kv);
pm=repmat(pn,[1 1 Kv]);
dyvar(:,2:Lv-1,:)=0.5.*pm(:,2:Lv-1,:).*(var(:,3:Lv,:)-var(:,1:Lv-2,:));

dyvar(:,1,:)=dyvar(:,2,:);
dyvar(:,Lv,:)=dyvar(:,Lv-1,:);
elseif length(siz)==4  % 4d variable
[Mv,Lv,Kv,Tv]=size(var);
dyvar=zeros(Mv,Lv,Kv,Tv);
pm=repmat(pn,[1 1 Kv Tv]);
dyvar(:,2:Lv-1,:,:)=0.5.*pm(:,2:Lv-1,:,:).*(var(:,3:Lv,:,:)-var(:,1:Lv-2,:,:));

dyvar(:,1,:,:)=dyvar(:,2,:,:);
dyvar(:,Lv,:,:)=dyvar(:,Lv-1,:,:);
end

return
