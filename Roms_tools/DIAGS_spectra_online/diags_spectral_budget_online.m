%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  this program checks the consistence of diagnostics terms using
% /home3/fastnet/pklein/ES_JUNE2006/s1_2km_100l_fpos_525d_t5_diags

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

addpath('../DIAGS_spectra_online')
addpath('../DIAGS_spectra')

mmean=0;
window=0; cff_tukey=1;
salinity=0;
wind=0;

root=['/data/models/JET/5K_ISO/'];
%root=[''];
model=['jet'];
dx=5.e3;

kchoicelist=[9:9];
lims=[1 100 1 400];
it1=60;it2=60;
rec_int_days=6;

g=9.81;rho0=1000;
cff_scale=1.e7;
xmin=0.9e-5; xmax=6.e-4; ymin=-1; ymax=1;  % limits for plots
print_and_keep=0;

his=-10;
dia=-9;
nstrt=0;nend=0;ninc=60;

loadmode = input('loadmode? (0/no 1/yes) ')

online_T=1;
online_P=1;
online_Cor=1;
online_Dv=1;
%=============================================================

if loadmode

  eval(['load SPECTRAL_KE_PK_2d_t' num2str(it1) '_to_t' num2str(it2) '.mat'])

else

kchoice=kchoicelist(1);

if mod(lims(2)-lims(1)+1,2)==1
   lims(2)=lims(2)-1;
end
if mod(lims(4)-lims(3)+1,2)==1
   lims(4)=lims(4)-1;
end

% Lateral grid
grd=rnt_gridload(model);

pm=grd.pm;pm=cut_var(pm,lims);
pn=grd.pn;pn=cut_var(pn,lims);
%lon=grd.lonr;lon=cut_var(lon,lims);
%lat=grd.latr;lat=cut_var(lat,lims);
f=grd.f;f=cut_var(f,lims);
% Vertical grid
hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;method=grd.Method;
h0=grd.h; h=h0(1,1); %h=cut_var(h,lims);
zr=zlevs(h,0,thetas,thetab,hc,N,'r');
zw=zlevs(h,0,thetas,thetab,hc,N,'w');
kk=kchoicelist;ll=length(kk);
dz(2:N)=zw(3:N+1)-zw(2:N);
dz(1)=dz(2);
dz_klist=dz(kk(1):kk(end))';
dz=dz(kchoice:N)';
%
kmin=kchoice-1;kmax=min(N,kchoice+1);
Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
zr=zlevs(h0,0,thetas,thetab,hc,N,'r');
zw=zlevs(h0,0,thetas,thetab,hc,N,'w');
zr=permute(zr,[2,3,1]);
zw=permute(zw,[2,3,1]);
dzw=zeros(size(zw));
dzw(:,:,2:N)=zr(:,:,2:N)-zr(:,:,1:N-1);
dzw(:,:,1)=dzw(:,:,2);
dzw=dzw(Lmin:Lmax,Mmin:Mmax,kmin:kmax);
dzr=zeros(size(zr));
dzr(:,:,2:N)=zw(:,:,3:N+1)-zw(:,:,2:N);
dzr(:,:,1)=dzr(:,:,2);
dzr=dzr(Lmin:Lmax,Mmin:Mmax,kmin:kmax);

ctldiags = ctlload(root(2:end),dia,model,nstrt,nend,ninc);
ctlhis   = ctlload(root(2:end),his,model,nstrt,nend,ninc+1);

ik=0;
for kchoice=kchoicelist
ik=ik+1;

%
% TEMPORAL LOOP ================================
%

for it=it1:it2

  if his==-10,
   ithis=it+1;
  else
   ithis=it;
  end;

  disp(['Treating rec ' num2str(it)])

%%%% CENTERED DIFFERENCE EVALUATION OF THE NL TERMS %%%%
%%%% computed as if in the rhs
  u=rnt_loadvar_partialz(ctlhis,ithis,'u',kchoice,kchoice); % v is 3d
  v=rnt_loadvar_partialz(ctlhis,ithis,'v',kchoice,kchoice); % v is 3d

  u=perm(u2rho(perm(u)));u=cut_var(u,lims);
  v=perm(v2rho(perm(v)));v=cut_var(v,lims);

  advu1=-u.*partialx(u,pm)-v.*partialy(u,pn);
  advv1=-u.*partialx(v,pm)-v.*partialy(v,pn);

%%%% UPSTREAM BIASED NL TERMS %%%%
%%%% computed as if in the rhs (- sign already included). 

  advu2=rnt_loadvar_partialz(ctldiags,it,'u_xadv',kchoice,kchoice);
  advu2=advu2+rnt_loadvar_partialz(ctldiags,it,'u_yadv',kchoice,kchoice);
  advu2=perm(u2rho(perm(advu2)));advu2=cut_var(advu2,lims);

  advv2=rnt_loadvar_partialz(ctldiags,it,'v_xadv',kchoice,kchoice);
  advv2=advv2+rnt_loadvar_partialz(ctldiags,it,'v_yadv',kchoice,kchoice);
  advv2=perm(v2rho(perm(advv2)));advv2=cut_var(advv2,lims);

%%%% implicit dissipation 
%%%% computed as in the rhs
  dissu=advu2-advu1; 
  dissv=advv2-advv1; 

%%%% explicit dissipation
%%%% computed as if in the rhs (- sign already included).
  dissu2=rnt_loadvar_partialz(ctldiags,it,'u_hmix',kchoice,kchoice);
  dissu2=perm(u2rho(perm(dissu2)));dissu2=cut_var(dissu2,lims);

  dissv2=rnt_loadvar_partialz(ctldiags,it,'v_hmix',kchoice,kchoice);
  dissv2=perm(v2rho(perm(dissv2)));dissv2=cut_var(dissv2,lims);

  dissu=dissu2;
  dissv=dissv2;

  if wind,
%%%% WIND WORK  %%%%
%%%% computed as if in the rhs
   sustr=rnt_loadvar(ctlhis,ithis,'sustr')./(rho0*sum(dz_klist));
   svstr=rnt_loadvar(ctlhis,ithis,'svstr')./(rho0*sum(dz_klist));
   us=rnt_loadvar_partialz(ctlhis,ithis,'u',N,N); 
   vs=rnt_loadvar_partialz(ctlhis,ithis,'v',N,N);
   sustr=perm(u2rho(perm(sustr)));sustr=cut_var(sustr,lims);
   svstr=perm(v2rho(perm(svstr)));svstr=cut_var(svstr,lims);
   us=perm(u2rho(perm(us)));us=cut_var(us,lims);
   vs=perm(u2rho(perm(vs)));vs=cut_var(vs,lims);
  end

%%%% temporal derivative
  if online_T==1,
   dudt=rnt_loadvar_partialz(ctldiags,it,'u_rate',kchoice,kchoice); % 
   dudt=perm(u2rho(perm(dudt)));dudt=cut_var(dudt,lims);

   dvdt=rnt_loadvar_partialz(ctldiags,it,'v_rate',kchoice,kchoice); %
   dvdt=perm(v2rho(perm(dvdt)));dvdt=cut_var(dvdt,lims);
  else
   uit1m=rnt_loadvar_partialz(ctlhis,it1-1,'u',kchoice,kchoice); % u is 3d
   uit1 =rnt_loadvar_partialz(ctlhis,it1  ,'u',kchoice,kchoice);
   uit2 =rnt_loadvar_partialz(ctlhis,it2  ,'u',kchoice,kchoice);
   uit2p=rnt_loadvar_partialz(ctlhis,it2+1,'u',kchoice,kchoice);
   vit1m=rnt_loadvar_partialz(ctlhis,it1-1,'v',kchoice,kchoice); % v is 3d
   vit1 =rnt_loadvar_partialz(ctlhis,it1  ,'v',kchoice,kchoice); 
   vit2 =rnt_loadvar_partialz(ctlhis,it2  ,'v',kchoice,kchoice); 
   vit2p=rnt_loadvar_partialz(ctlhis,it2+1,'v',kchoice,kchoice);
   dudt=0.5*(uit2+uit2p-uit1m-uit1)./(rec_int_days*(it2-it1+1)*86400);
   dvdt=0.5*(vit2+vit2p-vit1m-vit1)./(rec_int_days*(it2-it1+1)*86400);
   dudt=perm(u2rho(perm(dudt)));dudt=cut_var(dudt,lims);
   dvdt=perm(v2rho(perm(dvdt)));dvdt=cut_var(dvdt,lims);
  end

%%%% wb term
  cffb=-g/rho0; % conversion from density to buoyancy
  w=rnt_loadvar_partialz(ctlhis,ithis,'w',kchoice,kchoice); % w 
  w=cut_var(w,lims);
  w=w-mean(mean(w));
  if salinity,
    T=rnt_loadvar_partialz(ctlhis,ithis,'temp',kchoice,kchoice); % temp is 2d
    S= zeros(size(T)); %rnt_loadvar_partialz(ctlhis,ithis,'salt',kchoice,kchoice); % temp is 2d
    rho=rnt_rho_potential_nomex(T,S);clear T S
  else
   rho=rnt_loadvar_partialz(ctlhis,ithis,'temp',kchoice,kchoice); % rho is 2d
  end
  rho=cut_var(rho,lims);
  b=cffb.*(rho-mean(mean(rho))); % we remove the mean of b because it is not clear how ...

%%%% pressure term
  if online_P==1,
   dxp=rnt_loadvar_partialz(ctldiags,it,'u_Prsgrd',kchoice,kchoice); %
   dxp=perm(u2rho(perm(dxp)));dxp=cut_var(dxp,lims);

   dyp=rnt_loadvar_partialz(ctldiags,it,'v_Prsgrd',kchoice,kchoice); %
   dyp=perm(v2rho(perm(dyp)));dyp=cut_var(dyp,lims);

  elseif online_P==2,
   if salinity,
    temp=rnt_loadvar_partialz(ctlhis,ithis,'temp',kchoice,N);
    salt=rnt_loadvar_partialz(ctlhis,ithis,'salt',kchoice,N);
    temp=cut_var(temp,lims);
    salt=cut_var(salt,lims);
    rho=rnt_rho_potential_nomex(temp,salt);
   else
    rho=rnt_loadvar_partialz(ctlhis,ithis,'rho',kchoice,N);
   end
   rho=cut_var(rho,lims);
   zeta=rnt_loadvar(ctlhis,ithis,'zeta');
   zeta=cut_var(zeta,lims);
   rho=rho+rho0;
   [L,M]=size(zeta);
   dz3d=permute(repmat(dz,[1,L,M]),[2,3,1]);
   dz3d(:,:,end)=squeeze(dz3d(:,:,end))+zeta;
   rhow=rho;
   rhow(:,:,1:end-1)=0.5*(rho(:,:,1:end-1)+rho(:,:,2:end));
   pressure=9.8*sum(dz3d.*rhow,3);
   dxp=-partialx(pressure,pm)./rho0;  % watch out the - sign is included here
   dyp=-partialy(pressure,pn)./rho0;
  else
   zeta=rnt_loadvar(ctlhis,ithis,'zeta');
   zeta=cut_var(zeta,lims);
   h0=cut_var(h0,lims);
   zr=zlevs(h0,zeta,thetas,thetab,hc,N,'r');
   zw=zlevs(h0,zeta,thetas,thetab,hc,N,'w');
   zr=permute(zr,[2,3,1]);
   zw=permute(zw,[2,3,1]);
   zr=zr(:,:,kchoice:N);
   zw=zw(:,:,kchoice:N+1);
   if salinity,
     temp=rnt_loadvar_partialz(ctlhis,ithis,'temp',kchoice,N);
     salt=rnt_loadvar_partialz(ctlhis,ithis,'salt',kchoice,N);
     temp=cut_var(temp,lims);
     salt=cut_var(salt,lims);
     rho=rho_eos(temp,salt,zr); 
   else
     rho=rnt_loadvar_partialz(ctlhis,ithis,'rho',kchoice,N);
   end
   rho=cut_var(rho,lims);
   [dxp,dyp]=prsgrd(rho,zw,zr,pm,pn);
   dxp=perm(u2rho(perm(squeeze(dxp(:,:,1)))));
   dyp=perm(v2rho(perm(squeeze(dyp(:,:,1)))));
  end

%
% Coriolis
%
  if online_Cor==1,
    ucor=rnt_loadvar_partialz(ctldiags,it,'u_cor',kchoice,kchoice); %
    ucor=perm(u2rho(perm(ucor)));ucor=cut_var(ucor,lims);
    vcor=rnt_loadvar_partialz(ctldiags,it,'v_cor',kchoice,kchoice); %
    vcor=perm(v2rho(perm(vcor)));vcor=cut_var(vcor,lims);
  else
    dndx=zeros(size(pm)); dmde=zeros(size(pm));
    dndx(2:end-1,:)=0.5./pn(3:end,:)-0.5./pn(1:end-2,:);
    dmde(:,2:end-1)=0.5./pm(:,3:end)-0.5./pm(:,1:end-2);
    dndx(1,:)=dndx(2,:); dndx(end,:)=dndx(end-1,:);
    dmde(:,1)=dmde(:,2); dmde(:,end)=dmde(:,end-1);
    cff=f+(v.*dndx-u.*dmde).*pm.*pn;
    ucor=cff.*v;
    vcor=-cff.*u;
    %ucor(2:end-1,:)=0.5*ucor(2:end-1,:)+0.25*(ucor(1:end-2,:)+ucor(3:end,:));
    %vcor(:,2:end-1)=0.5*vcor(:,2:end-1)+0.25*(vcor(:,1:end-2)+vcor(:,3:end));
  end

%%%% vertical dissipation
  if online_Dv==1,
   ukpp=rnt_loadvar_partialz(ctldiags,it,'u_vmix',kchoice,kchoice); %
   ukpp=perm(u2rho(perm(ukpp)));ukpp=cut_var(ukpp,lims);

   vkpp=rnt_loadvar_partialz(ctldiags,it,'v_vmix',kchoice,kchoice); %
   vkpp=perm(v2rho(perm(vkpp)));vkpp=cut_var(vkpp,lims);
  elseif online_Dv==2,
   Kv=rnt_loadvar_partialz(ctlhis,ithis,'AKv',kchoice-1,kchoice+1);
   Kv=cut_var(Kv,lims);
   u=rnt_loadvar_partialz(ctlhis,ithis,'u',kchoice-1,kchoice+1); % u is 3d
   u=perm(u2rho_3d(perm(u))); u=cut_var(u,lims);
   v=rnt_loadvar_partialz(ctlhis,ithis,'v',kchoice-1,kchoice+1); % v is 3d
   v=perm(v2rho_3d(perm(v))); v=cut_var(v,lims);
   k=2;
   if kchoice==N,
    sustr=rnt_loadvar(ctl,it,'sustr')./rho0;
    sustr=perm(u2rho(perm(sustr))); sustr=cut_var(sustr,lims);
    svstr=rnt_loadvar(ctl,it,'svstr')./rho0;
    svstr=perm(v2rho(perm(svstr))); svstr=cut_var(svstr,lims);
    ukpp=[sustr - ...
           Kv(:,:,k  ).*(u(:,:,k  )-u(:,:,k-1))./dzw(:,:,k  )] ...
                                               ./dzr(:,:,k);
    vkpp=[svstr - ...
           Kv(:,:,k  ).*(v(:,:,k  )-v(:,:,k-1))./dzw(:,:,k  )] ...
                                               ./dzr(:,:,k);
   else
    ukpp=[Kv(:,:,k+1).*(u(:,:,k+1)-u(:,:,k  ))./dzw(:,:,k+1) - ...
          Kv(:,:,k  ).*(u(:,:,k  )-u(:,:,k-1))./dzw(:,:,k  )] ...
                                               ./dzr(:,:,k);
    vkpp=[Kv(:,:,k+1).*(v(:,:,k+1)-v(:,:,k  ))./dzw(:,:,k+1) - ...
          Kv(:,:,k  ).*(v(:,:,k  )-v(:,:,k-1))./dzw(:,:,k  )] ...
                                               ./dzr(:,:,k);
   end
   u=squeeze(u(:,:,k));
   v=squeeze(v(:,:,k));
  else
   Kv=rnt_loadvar_partialz(ctlhis,ithis,'AKv',kchoice,kchoice+1);
   Kv=cut_var(Kv,lims);
   uu=rnt_loadvar_partialz(ctlhis,ithis,'u',kchoice-1,kchoice+1); % u is 3d
   uu=perm(u2rho_3d(perm(uu)));
   uu=cut_var(uu,lims);
   zraux=squeeze(zr(1,1,kmin:kmax)); zwaux=squeeze(zw(1,1,kmin:kmax+1));
   dudz=partialz_r2w(uu,zraux);
   ukpp=Kv.*dudz;
   ukpp=partialz_w2r(ukpp,squeeze(zwaux));
   vv=rnt_loadvar_partialz(ctlhis,ithis,'v',kchoice-1,kchoice+1); % v is 3d
   vv=perm(v2rho_3d(perm(vv)));
   vv=cut_var(vv,lims);
   dvdz=partialz_r2w(vv,zraux);
   vkpp=Kv.*dvdz;
   vkpp=partialz_w2r(vkpp,squeeze(zwaux));
  end

%%%% vertical advection
  zadvu=rnt_loadvar_partialz(ctldiags,it,'u_vadv',kchoice,kchoice);
  zadvu=perm(u2rho(perm(zadvu)));zadvu=cut_var(zadvu,lims);

  zadvv=rnt_loadvar_partialz(ctldiags,it,'v_vadv',kchoice,kchoice);
  zadvv=perm(v2rho(perm(zadvv)));zadvv=cut_var(zadvv,lims);

%%%% residual 
  resu=-rnt_loadvar_partialz(ctldiags,it,'u_rate',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'u_xadv',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'u_yadv',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'u_cor',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'u_vmix',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'u_hmix',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'u_vadv',kchoice,kchoice)+ ...
        rnt_loadvar_partialz(ctldiags,it,'u_Prsgrd',kchoice,kchoice);
  resu=-resu; resu=perm(u2rho(perm(resu)));resu=cut_var(resu,lims);
  resv=-rnt_loadvar_partialz(ctldiags,it,'v_rate',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'v_xadv',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'v_yadv',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'v_cor',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'v_vmix',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'v_hmix',kchoice,kchoice)+...
        rnt_loadvar_partialz(ctldiags,it,'v_vadv',kchoice,kchoice)+ ...
        rnt_loadvar_partialz(ctldiags,it,'v_Prsgrd',kchoice,kchoice);
  resv=-resv; resv=perm(v2rho(perm(resv)));resv=cut_var(resv,lims);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  if it==it1 && kchoice==kchoicelist(1)
     e_rate=0;
     e_hdis=0;
     e_hadv=0;
     e_hdpr=0;
     e_zadv=0;
     e_zdis=0;
     e_cor=0;
     e_res=0;

     [L,M]=size(u); % works for non squared grids
     %cff_tukey=0.25;
     wdw1=tukeywin(L,cff_tukey)'; wdw2=tukeywin(M,cff_tukey)';
     wdw=wdw1'*wdw2;
     wdw=wdw;
  end

  if wind & kchoice==kchoicelist(1)
    % wind work only at the surface
    if mmean==1
     us=us-mean(mean(us));vs=vs-mean(mean(vs));
     sustr=sustr-mean(mean(sustr)); svstr=svstr-mean(mean(svstr));
    end
    if window==1     
     us=us.*wdw;vs=vs.*wdw;
     sustr=sustr.*wdw;svstr=svstr.*wdw;
    end
  end

  % apply window and remove spatial mean

%  e_rate=e_rate+mean(mean(u.*dudt+v.*dvdt));
%  e_hadv=e_hadv+mean(mean(u.*advu2+v.*advv2));
%  e_hdpr=e_hdpr+mean(mean(u.*dxp  +v.*dyp));
%  e_zadv=e_zadv+mean(mean(u.*zadvu+v.*zadvv));
%  e_zdis=e_zdis+mean(mean(u.*ukpp+v.*vkpp));
%  e_res=e_rate-(e_hadv+e_hdpr+e_zadv+e_zdis);

  if mmean==1    % b and w already done 
    u=u-mean(mean(u));v=v-mean(mean(v));
    advu1=advu1-mean(mean(advu1));advv1=advv1-mean(mean(advv1));
    advu2=advu2-mean(mean(advu2));advv2=advv2-mean(mean(advv2));
    dissu=dissu-mean(mean(dissu));dissv=dissv-mean(mean(dissv));
    dudt=dudt-mean(mean(dudt));dvdt=dvdt-mean(mean(dvdt));
    dxp=dxp-mean(mean(dxp));dyp=dyp-mean(mean(dyp));
    ukpp=ukpp-mean(mean(ukpp));vkpp=vkpp-mean(mean(vkpp));
    zadvu=zadvu-mean(mean(zadvu));zadvv=zadvv-mean(mean(zadvv));
    ucor=ucor-mean(mean(ucor));vcor=vcor-mean(mean(vcor));
    resu=resu-mean(mean(resu));resv=resv-mean(mean(resv));
  end
  if window==1
    u=u.*wdw;v=v.*wdw;
    advu1=advu1.*wdw;advv1=advv1.*wdw;
    advu2=advu2.*wdw;advv2=advv2.*wdw;
    dissu=dissu.*wdw;dissv=dissv.*wdw;
    dudt=dudt.*wdw;dvdt=dvdt.*wdw;
    dxp=dxp.*wdw;dyp=dyp.*wdw;
    ukpp=ukpp.*wdw;vkpp=vkpp.*wdw;
    zadvu=zadvu.*wdw;zadvv=zadvv.*wdw;
    ucor=ucor.*wdw;vcor=vcor.*wdw;
    resu=resu.*wdw;resv=resv.*wdw;
    b=b.*wdw;w=w.*wdw;
  end

  e_rate=e_rate+mean(mean(u.*dudt+v.*dvdt));
  e_hdis=e_hdis+mean(mean(u.*dissu+v.*dissv));
  e_hadv=e_hadv+mean(mean(u.*advu2+v.*advv2));
  e_hdpr=e_hdpr+mean(mean(u.*dxp  +v.*dyp));
  e_zadv=e_zadv+mean(mean(u.*zadvu+v.*zadvv));
  e_zdis=e_zdis+mean(mean(u.*ukpp+v.*vkpp));
  e_cor=e_cor+mean(mean(u.*ucor+v.*vcor));
  e_res=e_rate-(e_hadv+e_hdpr+e_zadv+e_zdis+e_cor);

%
%  SPECTRAL PROJECTION
%
  if it==it1 & kchoice==kchoicelist(1)
   cff=1/((L*M)^2);
   cff1=cff*dz_klist(ik)./sum(dz_klist);

   tmpamp0=zeros(M,L);
   tmpamp1=zeros(M,L); tmpamp2=zeros(M,L);
   tmpamp3=zeros(M,L); tmpamp4=zeros(M,L);
   tmpamp5=zeros(M,L); tmpamp6=zeros(M,L);
   tmpamp7=zeros(M,L); tmpamp8=zeros(M,L);
   tmpamp9=zeros(M,L); tmpamp10=zeros(M,L);
   tmpamp22=zeros(M,L);
  end

  fcoefu=fft2(perm(u),M,L); fcoefv=fft2(perm(v),M,L);
  tmpamp0=tmpamp0+cff1*real(conj(fcoefu).*fcoefu+conj(fcoefv).*fcoefv);         %  0  kinetic energy

  fcoeftmpu=fft2(perm(advu1),M,L); fcoeftmpv=fft2(perm(advv1),M,L);
  tmpamp1=tmpamp1+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  1  C4 horiz advection

  fcoeftmpu=fft2(perm(advu2),M,L); fcoeftmpv=fft2(perm(advv2),M,L);
  tmpamp2=tmpamp2+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  2  UP3 horiz advection

  fcoeftmpu=fft2(perm(dissu),M,L); fcoeftmpv=fft2(perm(dissv),M,L);
  tmpamp3=tmpamp3+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  3  horiz dissipation

  fcoeftmpu=fft2(perm(dudt),M,L);  fcoeftmpv=fft2(perm(dvdt),M,L);
  tmpamp4=tmpamp4+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  4  time rate of change

  fcoeftmpu=fft2(perm(dxp),M,L);  fcoeftmpv=fft2(perm(dyp),M,L);
  tmpamp5=tmpamp5+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  5  pressure gradient

  fcoeftmpu=fft2(perm(ukpp),M,L);  fcoeftmpv=fft2(perm(vkpp),M,L);
  tmpamp6=tmpamp6+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  6  vert dissipation

  fcoeftmpu=fft2(perm(zadvu),M,L);  fcoeftmpv=fft2(perm(zadvv),M,L);
  tmpamp7=tmpamp7+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  7  vert advection

  fcoeftmpu=fft2(perm(resu),M,L);  fcoeftmpv=fft2(perm(resv),M,L);
  tmpamp8=tmpamp8+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv);   %  8  residual

  fcoeftmpu=fft2(perm(b),M,L);  fcoeftmpv=fft2(perm(w),M,L);
  tmpamp9=tmpamp9+cff1*real(conj(fcoeftmpu).*fcoeftmpv);                        %  9  baroclinic conversion wb 

  fcoeftmpu=fft2(perm(ucor),M,L); fcoeftmpv=fft2(perm(vcor),M,L);
  tmpamp10=tmpamp10+cff1*real(conj(fcoefu).*fcoeftmpu+conj(fcoefv).*fcoeftmpv); % 10  coriolis

  if wind & kchoice==kchoicelist(1)
    % wind work only at the surface
    fcoefus=fft2(us); fcoefvs=fft2(vs);
    fcoeftmpu=fft2(sustr); fcoeftmpv=fft2(svstr);
    tmpamp22=tmpamp22+cff*real(conj(fcoefus).*fcoeftmpu+conj(fcoefvs).*fcoeftmpv);
    tmpamp22=1/(it2-it1+1)*tmpamp22;
    tmpamp22=reorganize_fft2d(tmpamp22);
    method=2;
    [amp22,count,ktmp,dk]=integ_fft2d(tmpamp22,dx,L,M,method);
  end

end % temporal loop

disp(['e_rate = ',num2str(e_rate)]);
disp(['e_hadv = ',num2str(e_hadv)]);
disp(['e_pgr  = ',num2str(e_hdpr)]);
disp(['e_vadv = ',num2str(e_zadv)]);
disp(['e_vdif = ',num2str(e_zdis)]);
disp(['e_cor  = ',num2str(e_cor)]);
disp(['e_res  = ',num2str(e_res)]);

end % vertical loop
disp(['sum tmpamp0 = ',num2str(sum(sum(tmpamp0)))]);
tmpamp0=1/(it2-it1+1)*tmpamp0;
tmpamp1=1/(it2-it1+1)*tmpamp1; tmpamp2=1/(it2-it1+1)*tmpamp2; 
tmpamp3=1/(it2-it1+1)*tmpamp3; tmpamp4=1/(it2-it1+1)*tmpamp4;
tmpamp5=1/(it2-it1+1)*tmpamp5; tmpamp6=1/(it2-it1+1)*tmpamp6; 
tmpamp7=1/(it2-it1+1)*tmpamp7; tmpamp8=1/(it2-it1+1)*tmpamp8; 
tmpamp9=1/(it2-it1+1)*tmpamp9; tmpamp10=1/(it2-it1+1)*tmpamp10;
disp(['tmpamp0 = ',num2str(mean(mean(tmpamp0)))]);
disp(['tmpamp1 = ',num2str(mean(mean(tmpamp1)))]);
disp(['tmpamp2 = ',num2str(mean(mean(tmpamp2)))]);
disp(['tmpamp3 = ',num2str(mean(mean(tmpamp3)))]);
disp(['tmpamp4 = ',num2str(mean(mean(tmpamp4)))]);
disp(['tmpamp5 = ',num2str(mean(mean(tmpamp5)))]);
disp(['tmpamp6 = ',num2str(mean(mean(tmpamp6)))]);
disp(['tmpamp7 = ',num2str(mean(mean(tmpamp7)))]);
tmpamp0=reorganize_fft2d(tmpamp0);
tmpamp1=reorganize_fft2d(tmpamp1); tmpamp2=reorganize_fft2d(tmpamp2);
tmpamp3=reorganize_fft2d(tmpamp3); tmpamp4=reorganize_fft2d(tmpamp4);
tmpamp5=reorganize_fft2d(tmpamp5); tmpamp6=reorganize_fft2d(tmpamp6);
tmpamp7=reorganize_fft2d(tmpamp7); tmpamp8=reorganize_fft2d(tmpamp8);
tmpamp9=reorganize_fft2d(tmpamp9); tmpamp10=reorganize_fft2d(tmpamp10);
disp(['tmpamp0 = ',num2str(mean(mean(tmpamp0)))]);
disp(['tmpamp1 = ',num2str(mean(mean(tmpamp1)))]);
disp(['tmpamp2 = ',num2str(mean(mean(tmpamp2)))]);
disp(['tmpamp3 = ',num2str(mean(mean(tmpamp3)))]);
disp(['tmpamp4 = ',num2str(mean(mean(tmpamp4)))]);
disp(['tmpamp5 = ',num2str(mean(mean(tmpamp5)))]);
disp(['tmpamp6 = ',num2str(mean(mean(tmpamp6)))]);
disp(['tmpamp7 = ',num2str(mean(mean(tmpamp7)))]);
% get 1D spectra
method=2;
[amp0,count,ktmp,dk]=integ_fft2d(tmpamp0,dx,M,L,method);
[amp1,count,ktmp,dk]=integ_fft2d(tmpamp1,dx,M,L,method);
[amp2,count,ktmp,dk]=integ_fft2d(tmpamp2,dx,M,L,method);
[amp3,count,ktmp,dk]=integ_fft2d(tmpamp3,dx,M,L,method);
[amp4,count,ktmp,dk]=integ_fft2d(tmpamp4,dx,M,L,method);
[amp5,count,ktmp,dk]=integ_fft2d(tmpamp5,dx,M,L,method);
[amp6,count,ktmp,dk]=integ_fft2d(tmpamp6,dx,M,L,method);
[amp7,count,ktmp,dk]=integ_fft2d(tmpamp7,dx,M,L,method);
[amp8,count,ktmp,dk]=integ_fft2d(tmpamp8,dx,M,L,method);
[amp9,count,ktmp,dk]=integ_fft2d(tmpamp9,dx,M,L,method);
[amp10,count,ktmp,dk]=integ_fft2d(tmpamp10,dx,M,L,method);

eval(['save SPECTRAL_KE_PK_nowindow_2d_t' num2str(it1) '_to_t' num2str(it2) '.mat amp0 amp1 amp2 amp3 amp4 amp5 amp6 amp7 amp8 amp9 amp10 ktmp count kchoicelist']);

end; % loadmode

amp11=amp4-(amp2+amp3+amp5+amp6+amp7+amp10);

istr=1;
disp(['e_rate rebuilt = ',num2str(sum(amp4(istr:end)))]);
disp(['e_hdif rebuilt = ',num2str(sum(amp3(istr:end)))]);
disp(['e_hadv rebuilt = ',num2str(sum(amp2(istr:end)))]);
disp(['e_pgr  rebuilt = ',num2str(sum(amp5(istr:end)))]);
disp(['e_vadv rebuilt = ',num2str(sum(amp7(istr:end)))]);
disp(['e_vdif rebuilt = ',num2str(sum(amp6(istr:end)))]);
disp(['e_cor  rebuilt = ',num2str(sum(amp10(istr:end)))]);
disp(['e_res  rebuilt = ',num2str(sum(amp11(istr:end)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% PLOT
%

lx=30; ly=15;
figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')


dk=ktmp(10)-ktmp(9);
istr=1;
K=ktmp(istr:end);
A(1,:) =cff_scale*amp1(istr:end)./dk;
A(2,:) =cff_scale*amp2(istr:end)./dk;
A(3,:) =cff_scale*amp3(istr:end)./dk;
A(4,:) =cff_scale*amp4(istr:end)./dk;
A(5,:) =cff_scale*amp5(istr:end)./dk;
A(6,:) =cff_scale*amp6(istr:end)./dk;
A(7,:) =cff_scale*amp7(istr:end)./dk;
A(8,:) =cff_scale*amp8(istr:end)./dk;
A(9,:) =cff_scale*amp9(istr:end)./dk;
A(10,:)=cff_scale*amp10(istr:end)./dk;
A(11,:)=cff_scale*amp11(istr:end)./dk;
if wind,
 A(22,:)=cff_scale*amp22(istr:end)./dk;
end

istr=1;
K(1:istr-1)=[];
hl1=plot(K,A(3,istr:end),K,A(2,istr:end),K,A(4,istr:end),K,A(5,istr:end)); hold on;
hl2=plot(K,A(6,istr:end),K,A(7,istr:end),K,A(10,istr:end),K,A(8,istr:end),K,A(11,istr:end)); hold on
set(hl1,'Linewidth',2)
set(hl2,'Linewidth',2,'Linestyle','--')
legend('D_H','A_H','T','P','D_V','A_V','COR','Res','Res(spectra)')

number1=K(1);
[coef,expon] =strread(strrep(sprintf('%E',number1),'E','#'),'%f#%f');
number2=floor(coef)*10^(expon);
set(gca,'Yscale','linear','Xscale','log');
%set(gca,'Ylim',[-1 1],'Xlim',[number2 4e-4]);
%set(gca,'Ylim',[ymin ymax],'Xlim',[xmin xmax]);
set(gca,'Xlim',[xmin xmax]);
xlabel('k [rad/m]');
[coef,expon] =strread(strrep(sprintf('%E',cff_scale),'E','#'),'%f#%f');
ylabel(['KE Tendencies [10^{-',num2str(expon),'} m^3/s^3]']);
set(gca,'fontsize',16)

hl=line([3e-6 2e-3],[0 0]); hold off
set(hl,'Linewidth',1,'Linestyle','-','Color','g')
%
x=[5000:-1000:2000,1000:-500:500,400:-100:100,80:-10:10,5];
LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');

if print_and_keep,
 outname=['KE_budget_spectra_online'];
 warning off
 eval(['print -painter -depsc2 ',outname,'.eps;'])
 warning on
 eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
 eval(['! mv ',outname,'.eps ',dirout_EPS]);
 eval(['! mv ',outname,'.jpg ',dirout_JPG]);
 eval(['! rm -f ',outname,'.eps']);
 eval(['! rm -f ',outname,'.jpg']);
end

