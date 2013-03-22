%==============================================================
% SPECTRAL KE BUDGET
% convention manu with first dim being xsi and second eta
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================
%close all
clear all

diags_params

loadmode = input('loadmode? (0/no 1/yes) ');

print_and_keep=0;
ilevel=1;
%================================================================

%for kchoice=klist;
  disp(' ')
  disp(['Vertical level K = ',num2str(kchoice)])
  disp(' ')

  dx=dx_lev(ilevel);
  if Nlev==1,
    lims=lims_lev;
    model=model0;
  else
    lims=squeeze(lims_lev(:,ilevel+1));
    model=[model0 num2str(ilevel)];
  end
  if mod(lims(2)-lims(1)+1,2)==1
     lims(2)=lims(2)+1;
  end
  if mod(lims(4)-lims(3)+1,2)==1
     lims(4)=lims(4)+1;
  end
  ctl=ctlload(root(2:end),his,model,0,lastfile_indx,Nrec);
  ctlhis=ctlload(root(2:end),his,model,0,lastfile_indx,Nrec);

  lx=30; ly=15;
  figure('units','centimeters','position', ...
         [0 0 lx ly],'paperpositionmode','auto')


% BAROCLINIC CONVERSION
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_I.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    [amp,count,ktmp,dk]=compute_spectral_flux_wb(ctl,model,dx,kchoice,window,lims,lit);
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_I.mat'])
  end
%
  amp_conv=amp; %save for 3D pressure work computation
  amptot=amp./dk;
  K=ktmp(istr:end);
  A=cff_scale*amp(istr:end)./dk(istr:end);
  if filtamp
    A=hanning_1d(A,filtcff);
  end
  K=interp(K,4);
  A=interp(A,4);
  pt1=plot(K,A); hold on;
  set(gca,'fontsize',16)
%
% HORIZONTAL ADVECTION
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_Ah.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    [amp,count,ktmp,dk]=compute_spectral_flux_advh(ctl,model,dx,kchoice,window,lims,lit);
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_Ah.mat'])
  end
%
  amptot=amptot+amp./dk;
  ampAH=amp;
%
% VERTICAL ADVECTION
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_Av.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    [amp,count,ktmp,dk]=compute_spectral_flux_advw(ctl,model,dx,kchoice,window,lims,lit);
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_Av.mat'])
  end
%
  amptot=amptot+amp./dk;
  amp=amp+ampAH;
  A=cff_scale*amp(istr:end)./dk(istr:end);
  if filtamp
    A=hanning_1d(A,filtcff);
  end
  A=interp(A,4);
  pt2=plot(K,A); hold on;
%
% HORIZONTAL DISSIPATION (x10)
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_Dh.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    [amp,count,ktmp,dk]=compute_spectral_flux_diffh(ctl,model,dx,kchoice,window,lims,lit);
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_Dh.mat'])
  end
%
  amptot=amptot+amp./dk;
  A=10*cff_scale*amp(istr:end)./dk(istr:end);
  if filtamp
    A=hanning_1d(A,filtcff);
  end
  A=interp(A,4);
  pt3=plot(K,A); hold on;
  set(gca,'fontsize',16)
%
% VERTICAL DISSIPATION
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_Dv.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    [amp,count,ktmp,dk]=compute_spectral_flux_diffv(ctl,model,dx,kchoice,window,lims,lit);
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_Dv.mat'])
  end
%
  amptot=amptot+amp./dk;
  ampDV=amp;
  A=cff_scale*amp(istr:end)./dk(istr:end);
  if filtamp
    A=hanning_1d(A,filtcff);
  end
  A=interp(A,4);
  pt4=plot(K,A); hold on;
  set(gca,'fontsize',16)
%
% 3D PRESSURE WORK
%
  if loadmode==1,
    eval(['!cp ',root,specname,'_KE_P.mat spectrum_KE.mat'])
    load spectrum_KE.mat
  else
    [amp,count,ktmp,dk]=compute_spectral_flux_P(ctl,model,dx,kchoice,window,lims,lit);
    amp=amp-amp_conv;
    save spectrum_KE amp count ktmp dk
    eval(['!mv spectrum_KE.mat ',root,specname,'_KE_P.mat'])
  end

  amptot=amptot+amp./dk;
  ampP=amp;
  A=cff_scale*amp(istr:end)./dk(istr:end);
  if filtamp
    A=hanning_1d(A,filtcff);
  end
  A=interp(A,4);
  pt5=plot(K,A); hold on;
  set(gca,'fontsize',16)

%
% TIME TENDENCY
%
  Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
  it1=lit(1); it2=lit(end)+1;
  for it=[it1 it2]
    disp(['Treating rec ' num2str(it)])
    u=rnt_loadvar_partialz(ctlhis,it,'u',kchoice,kchoice); % u is 3d
    v=rnt_loadvar_partialz(ctlhis,it,'v',kchoice,kchoice); % v is 3d
    u=perm(u2rho(perm(u)));
    u=squeeze(u(Lmin:Lmax,Mmin:Mmax));
    v=perm(v2rho(perm(v)));
    v=squeeze(v(Lmin:Lmax,Mmin:Mmax));
    u=u-mean(mean(u));
    v=v-mean(mean(v));
    [u,v,u,v]=windowing(u,v,u,v,window);
    [M,L]=size(u);
    fcoefu=fft2(u); fcoefv=fft2(v);
    tmpamp=0.5*(real(conj(fcoefu).*fcoefu+conj(fcoefv).*fcoefv));
    if it==it1,
      amp_ti=tmpamp;
    else
      amp_tf=tmpamp;
    end
  end % temporal loop
  tmpamp=(amp_tf-amp_ti)/((it2-it1)*86400);
  tmpamp=1/(L*M)^2*tmpamp;     
  tmpamp=reorganize_fft2d(tmpamp);  
  [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,2);

  amptot=amptot-amp./dk;
  pt7=plot(ktmp(istr:end),amp(istr:end)./dk(istr:end),'b');hold on;
  set(gca,'fontsize',16)

%%%% WIND WORK  %%%%
  windwork=0;
  if windwork
  for it=lit
  disp(['Treating rec ' num2str(it)])

  N=30;
  sustr=rnt_loadvar(ctlhis,it,'sustr')./(rho0*5);
  svstr=rnt_loadvar(ctlhis,it,'svstr')./(rho0*5);
  us=rnt_loadvar_partialz(ctlhis,it,'u',N,N);
  vs=rnt_loadvar_partialz(ctlhis,it,'v',N,N);
  sustr=perm(u2rho_2d(perm(sustr)));sustr=cut_var(sustr,lims);
  svstr=perm(v2rho_2d(perm(svstr)));svstr=cut_var(svstr,lims);
  us=perm(u2rho_2d(perm(us)));us=cut_var(us,lims);
  vs=perm(u2rho_2d(perm(vs)));vs=cut_var(vs,lims);
  us=us-mean(mean(us));vs=vs-mean(mean(vs));
  sustr=sustr-mean(mean(sustr)); svstr=svstr-mean(mean(svstr));
  [us,vs,sustr,svstr]=windowing(us,vs,sustr,svstr,window);
  fcoefus=fft2(us); fcoefvs=fft2(vs);
  fcoeftmpu=fft2(sustr); fcoeftmpv=fft2(svstr);
  tmpamp=tmpamp+real(conj(fcoefus).*fcoeftmpu+conj(fcoefvs).*fcoeftmpv);
  end
  tmpamp=1/((it2-it1+1)*(L*M)^2)*tmpamp;
  tmpamp=reorganize_fft2d(tmpamp);       % reorganize to be centered
  [amp,count,ktmp,dk]=integ_fft2d(tmpamp,dx,M,L,2);
  amp=cff_scale*amp;
  pt10=plot(ktmp(istr:end),amp(istr:end)./dk(istr:end),'b');
  hold on;
  set(gca,'fontsize',16)
  end

%
%----------------------------------------------------------------------------------
%
% FINALIZE PLOT
%
  line([3e-6 2e-3],[0 0],'Color','r'); hold off;
  set(pt1,'Linewidth',2,'Linestyle','-' ,'Color','k')
  set(pt2,'Linewidth',2,'Linestyle','--','Color','k')
  set(pt3,'Linewidth',2,'Linestyle','--' ,'Color',.7*[1 1 1])
  set(pt4,'Linewidth',2,'Linestyle','-' ,'Color',.7*[1 1 1])
  set(pt5,'Linewidth',2,'Linestyle',':' ,'Color',.7*[1 1 1])
  legend('I','A','D_H (x10)','D_V','P','Location','SouthEast')

  number1=K(1);
  [coef,expon] =strread(strrep(sprintf('%E',number1),'E','#'),'%f#%f');
  number2=floor(coef)*10^(expon);
  set(gca,'Yscale','linear','Xscale','log');
  set(gca,'Ylim',[-1 1],'Xlim',[number2 4e-4]);
  xlabel('k [rad/m]'); 
  [coef,expon] =strread(strrep(sprintf('%E',cff_scale),'E','#'),'%f#%f');
  ylabel(['KE Tendencies [10^{-',num2str(expon),'} m^3/s^3]'])
  set(gca,'fontsize',16)
%
  x=[2000,1500,1000,700,500,400:-100:100,80:-10:10,5];
  LinkTopAxisData((6.28*1.e-3)./x,x,'Wavelength [km]');
%
  if print_and_keep,
   outname=['KE_budget_spectra'];
   warning off
   eval(['print -painter -depsc2 ',outname,'.eps;'])
   warning on
   eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])
   eval(['! mv ',outname,'.eps ',dirout_EPS]);
   eval(['! mv ',outname,'.jpg ',dirout_JPG]);
   eval(['! rm -f ',outname,'.eps']);
   eval(['! rm -f ',outname,'.jpg']);
  end

%end % kchoice


