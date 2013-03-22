%==============================================================
% COMPUTE T FLUX DIAGNOSTICS OF EMBEDDED GRID SOLUTIONS
% Capet 2008
% modif: Capet, Marchesiello 2010
%==============================================================

clear all
close all

diags_params

print_and_keep = 0;
kchoice        = 1;
filter_smart   = 1;
navg           = 1;      % time filter: avg period =navg*2+1 days
R              = 100.e3; % space filter: decorrelation radius [m]

loaddata       = 0;
%=============================================================

for ilevel=1:Nlev;

 nhan=3.^(ilevel-1);
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
 ctlavg=ctlload(root(2:end),his,model,0,lastfile_indx,Nrec);

 imin=lims(1);imax=lims(2);jmin=lims(3);jmax=lims(4);
 ikm=2*min(lims(2)-lims(1),lims(4)-lims(3));
 L=imax-imin+1;M=jmax-jmin+1;

 grd=rnt_gridload(model);
 theta_s=grd.thetas;theta_b=grd.thetab;hc=grd.hc;Method=grd.Method;N=grd.N;
 kmin=kchoice; kmax=N;

 pm=grd.pm;    pm=pm(imin:imax,jmin:jmax);
 pn=grd.pn;    pn=pn(imin:imax,jmin:jmax);
 %lon=grd.lonr; lon=lon(imin:imax,jmin:jmax);
 %lat=grd.latr; lat=lat(imin:imax,jmin:jmax);
 xr=grd.xr;    xr=xr(imin:imax,jmin:jmax);
 yr=grd.yr;    yr=yr(imin:imax,jmin:jmax);
 h=grd.h;      h=h(imin:imax,jmin:jmax);
 zr=perm(zlevs(h(1,1),0,theta_s,theta_b,hc,N,'r'));
 pz=zeros(size(zr));
 pz(2:N-1)=0.5.*(zr(3:N)-zr(1:N-2));pz(1)=pz(2);pz(N)=pz(N-1);
 pz=1./pz; pz=pz(kmin:kmax); zr=zr(kmin:kmax);

 if loaddata,

  if ilevel==1,
   load dwtdz_l1
  elseif ilevel==2,
   load dwtdz_l2
  elseif ilevel==3,
   load dwtdz_l3
  end

 else  % compute and save data

  N=kmax-kmin+1;
  wt=zeros(L,M,N);
  hbl=zeros(L,M);

  it1=lit(1); it2=lit(end);

  for it=it1+navg:it2-navg

    disp(['processing rec ' num2str(it) ]);

    % Total T, Filtered and subfiltered T
    tmp=rnt_loadvar_partialz(ctlavg,it,'temp',kmin,kmax);
    ttot=tmp(imin:imax,jmin:jmax,:);
    tfilt=zeros(size(ttot));
    for k=kmin:kmax-1
      if filter_smart,
        tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'temp',k,k)),3);
        tmp=tmp(imin:imax,jmin:jmax);
        tfilt(:,:,k-kmin+1)=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
      else
        tmp=low_filter1(ctlavg,'temp',it,k,nhan,navg);
        tfilt(:,:,k-kmin+1)=tmp(imin:imax,jmin:jmax);
      end
    end
    tsub=ttot-tfilt;

    % Total, filtered and subfiltered  w
    tmp=rnt_loadvar_partialz(ctlavg,it,'w',kmin,kmax);
    tmp(tmp>1e30)=NaN;
    wtot=tmp(imin:imax,jmin:jmax,:);
    wfilt=zeros(L,M,N);
    for k=kmin:kmax-1
      if filter_smart,
        tmp=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,'w',k,k)),3);
        tmp=tmp(imin:imax,jmin:jmax);
        wfilt(:,:,k-kmin+1)=perm(filter_pp(perm(xr),perm(yr),perm(tmp),R));
      else
        tmp=low_filter1(ctlavg,'w',it,k,nhan,navg);
        wfilt(:,:,k-kmin+1)=tmp(imin:imax,jmin:jmax);
      end
    end
    wsub=wtot-wfilt;

    % T flux associated with submesoscale
    tmp=wsub.*tsub;
    wt=wt+tmp;

    % HBL
    %tmp=rnt_loadvar(ctl,it,'hbl');
    hbl=hbl+tmp(imin:imax,jmin:jmax);

  end

  tmp=length([1+navg:length(lit)-navg]);
  wt=wt./tmp;
  hbl=hbl./tmp;
  dwtdz=-partialz(wt,pz');

  dwtdz=mean(dwtdz,1);
  hbl=mean(hbl,1);

  x=squeeze(yr(1,:));
  z=zr;
  var=perm(squeeze(dwtdz));

  if ilevel==1,
   save dwtdz_l1 x z var hbl
  elseif ilevel==2,
   save dwtdz_l2 x z var hbl
  elseif ilevel==3,
   save dwtdz_l3 x z var hbl
  end
  
 end % load or save data

 yr1d=squeeze(yr(1,:));
 for k=1:length(z)
    var(k,:)=filter_my_1D(yr1d',squeeze(var(k,:))',50.e3);
 end


 %  =============  PLOT ====================
 lx=30; ly=15;
 figure('units','centimeters','position', ...
          [0 0 lx ly],'paperpositionmode','auto')

 cff_scale=30*86400; cmin=-0.1; cint=0.01; cmax=0.1;

 var=cff_scale*var;
 contourf(x,z,var,[cmin:cint:cmax],'k:'); hold on;
 caxis([cmin cmax]); colorbar
 hl=line(x,-hbl); hold off
 set(hl,'Linewidth',2,'Linestyle','-','Color','k')
 set(gca,'fontsize',16)

 set(gca,'Yscale','linear','Xscale','linear');
% set(gca,'Xlim',[-3.5 8],'Ylim',[-85 -2.5]);
 ylabel('Depth')
 xlabel('X')
 set(gca,'fontsize',16)

 if print_and_keep,
   outname=['dwtdz_TIW_l',num2str(ilevel),'_s',num2str(season)];
   warning off
   eval(['print -painter -depsc2 ',outname,'.eps;'])
   warning on
   eval(['! convert -quality 100 ',outname,'.eps ',outname,'.jpg'])

   eval(['! mv ',outname,'.eps ',dirout_EPS]);
   eval(['! mv ',outname,'.jpg ',dirout_JPG]);
   eval(['! rm -f ',outname,'.eps']);
   eval(['! rm -f ',outname,'.jpg']);
 end

end




