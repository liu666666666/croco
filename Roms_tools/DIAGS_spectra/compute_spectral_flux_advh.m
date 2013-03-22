function [amp,count,ktmp,dk]=compute_spectral_flux_advh(ctl,model,dx,kchoice,window,lims,lit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [amp,count,ktmp]=compute_spectral_flux_adv(ctl,dx,kchoice,window,lims,lit);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  this function computes velocity power spectrum for a netcdf file series. 
%%%% inputs: 
%           ctl    : structure with netcdf files loaded
%           dx     : meshgrid size (!!! supposed to be uniform !!!)
%           kchoice: sigma level to be processed (z level not implemented)
%           window : windowing option. There are 2 options. 
%                       apply a filter (rectangular 1, triangular 2 or hanning 3)
%                       apply a 2d dct    (0) or 1d dct (xsi dir, -1)
%           lims   : limits of the domain to be processed (4 elements vector with
%                       lims(1) and lims(2) mins and max in xsi direction. 
%           t1,t2  : time limits in the series. Default: 1 and length(ctl.time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
pi=3.1415926535;

grd=rnt_gridload(model);
% horizontal grid
pm0=grd.pm;pm=squeeze(pm0(Lmin:Lmax,Mmin:Mmax));
pn0=grd.pn;pn=squeeze(pn0(Lmin:Lmax,Mmin:Mmax));
h0=grd.h; h=squeeze(h0(Lmin:Lmax,Mmin:Mmax));
f0=grd.f; f=squeeze(f0(Lmin:Lmax,Mmin:Mmax));
% Vertical grid
hc=grd.hc;thetas=grd.thetas;thetab=grd.thetab;N=grd.N;h=grd.h;method=grd.Method;
zw=zlevs(h0,0,thetas,thetab,hc,N,'w');
Hz0(:,:)=zw(kchoice+1,:,:)-zw(kchoice,:,:);
Hz=squeeze(Hz0(Lmin:Lmax,Mmin:Mmax));

if nargin==5
%  it1=1; it2=length(ctl.time);
  lit=[1:length(ctl.time)];
end

for iit=1:length(lit)
  it=lit(iit);
  disp(['Treating rec ' num2str(it)])

  % Vertical grid
  %zeta0=rnt_loadvar(ctl,it,'zeta');
  %zw=zlevs(h0,zeta0,thetas,thetab,hc,N'w');
  %Hz0(:,:)=zw(kchoice+1,:,:)-zw(kchoice,:,:);
  %Hz=squeeze(Hz0(Lmin:Lmax,Mmin:Mmax));

  u0=rnt_loadvar_partialz(ctl,it,'u',kchoice,kchoice); % u is 3d
  u=perm(u2rho(perm(u0)));
  u=squeeze(u(Lmin:Lmax,Mmin:Mmax));

  v0=rnt_loadvar_partialz(ctl,it,'v',kchoice,kchoice); % v is 3d
  v=perm(v2rho(perm(v0)));
  v=squeeze(v(Lmin:Lmax,Mmin:Mmax));
%
% advection with 4th order centered scheme
%
  advC2_xav=0;
  if advC2_xav==1,
    advu=-u.*partialx(u,pm)-v.*partialy(u,pn);
    advv=-u.*partialx(v,pm)-v.*partialy(v,pn);
  else
    %[advu,advv]=roms_AdvH_C2(u0,v0,Hz0,pn0,pm0);
   advC2=0; % 0 --> 4th order scheme
   [Lu,Mu]=size(u0); L=Lu;M=Mu-1;Lm=L-1;Mm=M-1;Lm2=L-2;Mm2=M-2;
   advu=zeros(size(u0)); advv=zeros(size(v0));
   if advC2,
    advu(2:Lm2,2:Mm2)=-0.25*u0(2:Lm2,2:Mm2) ...
                    .*(u0(3:Lm,2:Mm2)-u0(1:Lm2-1,2:Mm2)) ...
                    .*(pm0(2:Lm2,2:Mm2)+pm0(3:Lm,2:Mm2)) ...
                    -0.0625*(v0(2:Lm2  ,2:Mm2  )+v0(1:Lm2-1,2:Mm2  )+ ...
                             v0(1:Lm2-1,1:Mm2-1)+v0(2:Lm2  ,1:Mm2-1)) ...
                    .*(u0(2:Lm2,3:Mm)-u0(2:Lm2,1:Mm2-1)) ...
                    .*(pn0(2:Lm2,2:Mm2)+pn0(3:Lm,2:Mm2));
    advv(2:Lm2,2:Mm2)=-0.25*v0(2:Lm2,2:Mm2) ...
                    .*(v0(2:Lm2,3:Mm)-v0(2:Lm2,1:Mm2-1)) ...
                    .*(pn0(2:Lm2,2:Mm2)+pn0(2:Lm2,3:Mm)) ...
                    -0.0625*(u0(2:Lm2,2:Mm2  )+u0(3:Lm ,2:Mm2  )+ ...
                             u0(3:Lm ,1:Mm2-1)+u0(2:Lm2,1:Mm2-1)) ...
                    .*(v0(3:Lm,2:Mm2)-v0(1:Lm2-1,2:Mm2)) ...
                    .*(pm0(2:Lm2,2:Mm2)+pm0(2:Lm2,3:Mm));
   else,
    advu(2:Lm2,2:Mm2)=-0.25*u0(2:Lm2,2:Mm2) ...
                    .*( (u0(3:Lm,2:Mm2)-u0(1:Lm2-1,2:Mm2)) ...
                        -1/6*(u0(4:L,2:Mm2)-3.*u0(3:Lm,2:Mm2)+ ...       % 4th order term
                              3.*u0(2:Lm2,2:Mm2)-u0(1:Lm2-1,2:Mm2)) ...  %
                      ) ...
                    .*(pm0(2:Lm2,2:Mm2)+pm0(3:Lm,2:Mm2)) ...
                    -0.0625*(v0(2:Lm2  ,2:Mm2  )+v0(1:Lm2-1,2:Mm2  )+ ...
                             v0(1:Lm2-1,1:Mm2-1)+v0(2:Lm2  ,1:Mm2-1)) ...
                    .*( (u0(2:Lm2,3:Mm)-u0(2:Lm2,1:Mm2-1)) ...
                        -1/6*(u0(2:Lm2,4:M)-3.*u0(2:Lm2,3:Mm)+ ...       % 4th order term
                              3.*u0(2:Lm2,2:Mm2)-u0(2:Lm2,1:Mm2-1)) ...  %
                      ) ...
                     .*(pn0(2:Lm2,2:Mm2)+pn0(3:Lm,2:Mm2));

    advv(2:Lm2,2:Mm2)=-0.25*v0(2:Lm2,2:Mm2) ...
                    .*( (v0(2:Lm2,3:Mm)-v0(2:Lm2,1:Mm2-1)) ...
                        -1/6*(v0(2:Lm2,4:M)-3.*v0(2:Lm2,3:Mm)+ ...       % 4th order term
                              3.*v0(2:Lm2,2:Mm2)-v0(2:Lm2,1:Mm2-1)) ...  %
                      ) ...
                     .*(pn0(2:Lm2,2:Mm2)+pn0(2:Lm2,3:Mm)) ...
                    -0.0625*(u0(2:Lm2,2:Mm2  )+u0(3:Lm ,2:Mm2  )+ ...
                             u0(3:Lm ,1:Mm2-1)+u0(2:Lm2,1:Mm2-1)) ...
                    .*( (v0(3:Lm,2:Mm2)-v0(1:Lm2-1,2:Mm2)) ...
                        -1/6*(v0(4:L,2:Mm2)-3.*v0(3:Lm,2:Mm2)+ ...       % 4th order term
                              3.*v0(2:Lm2,2:Mm2)-v0(1:Lm2-1,2:Mm2)) ...  %
                      ) ...
                     .*(pm0(2:Lm2,2:Mm2)+pm0(2:Lm2,3:Mm));
   end
   advu=perm(u2rho(perm(advu)));
   advv=perm(v2rho(perm(advv)));
   advu=squeeze(advu(Lmin:Lmax,Mmin:Mmax));
   advv=squeeze(advv(Lmin:Lmax,Mmin:Mmax));

   curvgrid=1;
   if curvgrid,
    dndx=zeros(size(pm)); dmde=zeros(size(pm));
    dndx(2:end-1,:)=0.5./pn(3:end,:)-0.5./pn(1:end-2,:);
    dmde(:,2:end-1)=0.5./pm(:,3:end)-0.5./pm(:,1:end-2);
    dndx(1,:)=dndx(2,:); dndx(end,:)=dndx(end-1,:);
    dmde(:,1)=dmde(:,2); dmde(:,end)=dmde(:,end-1);
    cff=(v.*dndx-u.*dmde).*pm.*pn;
    ucor=cff.*v;
    ucor(2:end-1,:)=0.5*ucor(2:end-1,:)+0.25*(ucor(1:end-2,:)+ucor(3:end,:));
    vcor=-cff.*u;
    vcor(:,2:end-1)=0.5*vcor(:,2:end-1)+0.25*(vcor(:,1:end-2)+vcor(:,3:end));

    advu=advu+ucor;
    advv=advv+vcor;
   end

  end

  [u,v,advu,advv]=windowing(u,v,advu,advv,window);

  [M,L]=size(u);
  if it==lit(1);
     tmpamp=zeros(M,L); 
  end
  fcoefadvu=fft2(advu); fcoefadvv=fft2(advv);
  fcoefu=fft2(u); fcoefv=fft2(v);
  tmpamp=tmpamp+real(conj(fcoefu).*fcoefadvu+conj(fcoefv).*fcoefadvv);

end % temporal loop

%----------------------------------------------------------------
% if check parseval's relationship
% lhswt=mean(mean(w.*T));
% rhswt=1/(Lu*Mu)^2*sum(sum(tmpampwT));
%
% do the reorganization and processing of results.
tmpamp=1/(length(lit)*(L*M)^2)*tmpamp; % parseval equality
tmpamp=reorganize_fft2d(tmpamp);       % reorganize to be centered

%
% Changes to get 1D spectra
%
new=1;
if new==1
  method=2;
  [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,method);
else
  Lh=L/2;Mh=M/2;i2L=1/L^2;i2M=1/M^2;
  kk=repmat(perm([1:L])-Lh,[1 M]); ll=repmat([1:M]-Mh,[L 1]);
  Kh=i2L*kk.^2+i2M*ll.^2; Kh=6.28/(dx)*sqrt(Kh);Kh=perm(Kh);
  leng=L;
  aux=leng/Kh(M,L);ktmp=[1:leng]/aux;ktmp=perm(ktmp);
  amp=zeros(leng,1);count=zeros(leng,1);

  ik = min(round(aux*Kh)+1,leng);
  for ind=1:leng
    II=find(ik==ind);
    amp(ind)=sum(tmpamp(II));
    count(ind)=length(II);
  end
end

return
