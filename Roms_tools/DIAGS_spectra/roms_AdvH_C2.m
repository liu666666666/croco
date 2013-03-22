function [advu,advv]=roms_hadv_C2(u,v,Hz,pn,pm)

%=======================================================================
%  ROMS 3rd Upstream-biaised horizontal advection of momentum
%=======================================================================
%  Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
%  of tensor of momentum flux due to horizontal advection; after that
%  compute horizontal advection terms in units of [m/s2]
%
%  Note: The difficulty is in dealing with 0 array indices not allowed by 
%  matlab. This requires attention in the computation of off-diagonal 
%  components which mixes staggered U and V type of variables.
%========================================================================

u=permute(u,[2,1]);
v=permute(v,[2,1]);
Hz=permute(Hz,[2,1]);
pm=permute(pm,[2,1]);
pn=permute(pn,[2,1]);

[Mu Lu]=size(u);
[Mv Lv]=size(v);
[Mr Lr]=size(pm);
L=Lr-1; M=Mr-1; Lm=L-1; Mm=M-1; 

Huon(1:Mu,1:Lu)=(Hz(1:Mr,2:Lr)+Hz(1:Mr,1:Lr-1))./ ...
                    (pn(1:Mr,2:Lr)+pn(1:Mr,1:Lr-1)).* ...
                     u(1:Mu,1:Lu);
Hvom(1:Mv,1:Lv)=(Hz(2:Mr,1:Lr)+Hz(1:Mr-1,1:Lr))./ ...
                    (pm(2:Mr,1:Lr)+pm(1:Mr-1,1:Lr)).* ...
                     v(1:Mv,1:Lv);
mnoHu(1:Mu,1:Lu)=0.5*(pm(1:Mr,2:Lr)+pm(1:Mr,1:Lr-1)).* ...
                     (pn(1:Mr,2:Lr)+pn(1:Mr,1:Lr-1))./ ...
                     (Hz(1:Mr,2:Lr)+Hz(1:Mr,1:Lr-1));
mnoHv(1:Mv,1:Lv)=0.5*(pm(2:Mr,1:Lr)+pm(1:Mr-1,1:Lr)).* ...
                     (pn(2:Mr,1:Lr)+pn(1:Mr-1,1:Lr))./ ...
                     (Hz(2:Mr,1:Lr)+Hz(1:Mr-1,1:Lr));

cffX=zeros(Mr,Lr);
cffE=zeros(Mr,Lr);
UFx=zeros(Mr,Lr);
UFe=zeros(Mr,Lr);
VFx=zeros(Mr,Lr);
VFe=zeros(Mr,Lr);
MU_Xadv=zeros(Mu,Lu);
MU_Yadv=zeros(Mu,Lu);
MV_Xadv=zeros(Mv,Lv);
MV_Yadv=zeros(Mv,Lv);

%===============================
% Diagonal component UFx
%
cffX(:,1:Lm)=u(:,1:Lm)+u(:,2:L);
UFx(2:M,1:Lm)=0.25*cffX(2:M,1:Lm).* ...
                   (Huon(2:M,1:Lm)+Huon(2:M,2:L));

%================================
% Diagonal component VFe
%
cffE(1:Mm,:)=v(1:Mm,:)+v(2:M,:);
VFe(1:Mm,2:L)=0.25*cffE(1:Mm,2:L).* ...
                   (Hvom(1:Mm,2:L)+Hvom(2:M,2:L));

%=================================
% off-diagonal component UFe
%
cffX(2:M+1,1:L)=u(2:M+1,:)+u(1:M,:);
cffE(2:M+1,2:Lm)=Hvom(:,3:L)+Hvom(:,2:Lm);
UFe(2:M+1,2:Lm)=0.25*cffX(2:M+1,2:Lm).* ...
                     cffE(2:M+1,2:Lm); 

%==================================
% off-diagonal VFx
%
cffE(1:M,2:L+1)=v(:,2:L+1)+v(:,1:L);
cffX(2:Mm,2:L+1)=Huon(3:M,:)+Huon(2:Mm,:);
VFx(2:Mm,2:L+1)=0.25*cffE(2:Mm,2:L+1).* ...
                     cffX(2:Mm,2:L+1); 

%===================================
% Compute advection terms.
% Then Divide advection terms by the cell volume Hz/(pm*pn).
% There after the unit of diag terms are :
% (unit of velocity) / s  =  [m/s2]
%
MU_Xadv(2:M,2:Lm) = -UFx(2:M,2:Lm)+UFx(2:M,1:Lm-1);
MU_Yadv(2:M,2:Lm) = -UFe(3:M+1,2:Lm)+UFe(2:M,2:Lm);
MV_Xadv(2:Mm,2:L) = -VFx(2:Mm,3:L+1)+VFx(2:Mm,2:L);
MV_Yadv(2:Mm,2:L) = -VFe(2:Mm,2:L)+VFe(1:Mm-1,2:L);

advu = MU_Xadv.*mnoHu+MU_Yadv.*mnoHu;
advv = MV_Xadv.*mnoHv+MV_Yadv.*mnoHv;

adv=0;
if adv,
  divH=zeros(Mr,Lr);
  divH(2:Mm,2:Lm)=Huon(2:Mm,3:L)-Huon(2:Mm,2:Lm)+Hvom(3:M,2:Lm)-Hvom(2:Mm,2:Lm);
  divH_u=rho2u_2d(divH);
  divH_v=rho2v_2d(divH);
  advu(2:M,2:Lm)=advu(2:M,2:Lm)+u(2:M,2:Lm).*divH_u(2:M,2:Lm).*mnoHu(2:M,2:Lm);
  advv(2:Mm,2:L)=advv(2:Mm,2:L)+v(2:Mm,2:L).*divH_v(2:Mm,2:L).*mnoHv(2:Mm,2:L);
end

advu=permute(advu,[2,1]);
advv=permute(advv,[2,1]);

return
