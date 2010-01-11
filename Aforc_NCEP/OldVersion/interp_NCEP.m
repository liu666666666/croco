function interp_NCEP(NCEP_dir,Y,M,Roa,interp_method,...
                     lon1,lat1,mask1,tin,...
		     nc_frc,nc_blk,lon,lat,angle,tout)

%
% Read the local NCEP files and perform the interpolations
%
% Pierrick 2005
% Menkes 2007
% Update, Feb-2008, J. Lefevre to use Nomads3 server
% Update, March-2009, G. Cambon & F. Marin Compute net solar short wave
%---------------------------------------------------------------------------------
%
% 1: Air temperature: Convert from Kelvin to Celsius
%
nc=netcdf([NCEP_dir,'tmp2m_Y',num2str(Y),'M',num2str(M),'.nc']);
tair=squeeze(nc{'tmp2m'}(tin,:,:));
close(nc);
%tair=get_missing_val(lon1,lat1,mask1.*tair);
tair=get_missing_val(lon1,lat1,mask1.*tair,nan,Roa,nan);
tair=tair-273.15;
tair=interp2(lon1,lat1,tair,lon,lat,interp_method);
%
% 2: Relative humidity: Convert from % to fraction
%
% Get Specific Humidity [Kg/Kg]
%
nc=netcdf([NCEP_dir,'spfh2m_Y',num2str(Y),'M',num2str(M),'.nc']);
shum=squeeze(nc{'spfh2m'}(tin,:,:));
close(nc);
%shum=get_missing_val(lon1,lat1,mask1.*shum);
shum=get_missing_val(lon1,lat1,mask1.*shum,nan,Roa,nan);
shum=interp2(lon1,lat1,shum,lon,lat,interp_method);
%
% computes specific humidity at saturation (Tetens  formula)
% (see air_sea tools, fonction qsat)
%
rhum=shum./qsat(tair);
%
% 3: Precipitation rate: Convert from [kg/m^2/s] to cm/day
%
nc=netcdf([NCEP_dir,'pratesfc_Y',num2str(Y),'M',num2str(M),'.nc']);
prate=squeeze(nc{'pratesfc'}(tin,:,:));
close(nc);
prate=get_missing_val(lon1,lat1,mask1.*prate,nan,Roa,nan);
prate=prate*0.1*(24*60*60.0);
prate=interp2(lon1,lat1,prate,lon,lat,interp_method);
prate(abs(prate)<1.e-4)=0;
%
% 4: Net shortwave flux: [W/m^2]
%      ROMS convention: downward = positive
%
% Downward solar shortwave
%
nc=netcdf([NCEP_dir,'dswrfsfc_Y',num2str(Y),'M',num2str(M),'.nc']);
  dswrs=squeeze(nc{'dswrfsfc'}(tin,:,:));
  close(nc);
%  
% Upward solar shortwave
% 
  nc=netcdf([NCEP_dir,'uswrfsfc_Y',num2str(Y),'M',num2str(M),'.nc']);
  uswrs=squeeze(nc{'uswrfsfc'}(tin,:,:));
  close(nc);
  
%%%%%%  dswrs=get_missing_val(lon1,lat1,mask1.*dswrs);
  dswrs=get_missing_val(lon1,lat1,mask1.*dswrs,nan,Roa,nan);
  uswrs=get_missing_val(lon1,lat1,mask1.*uswrs,nan,Roa,nan);
%  
%Net solar shortwave radiation  
%
   nswrs=dswrs - uswrs;
%----------------------------------------------------  
% GC le 31 03 
%  radsw is NET soalr shortwave radiation
%  no more downward only solar radiation
% GC  bug fixe by F. Marin IRD/LEGOS
%-----------------------------------------------------

% $$$   radsw=interp2(lon1,lat1,dswrs,lon,lat,interp_method);
% $$$   radsw(radsw<1.e-10)=0;
  
  radsw=interp2(lon1,lat1,nswrs,lon,lat,interp_method);
  radsw(radsw<1.e-10)=0;
  
%
% 5: Net outgoing Longwave flux:  [W/m^2]
%      ROMS convention: positive upward (opposite to nswrs)
%
% get downward longwave flux [W/m^2]
  nc=netcdf([NCEP_dir,'dlwrfsfc_Y',num2str(Y),'M',num2str(M),'.nc']);
  dlwrf=squeeze(nc{'dlwrfsfc'}(tin,:,:));
  close(nc);
%  dlwrf=get_missing_val(lon1,lat1,mask1.*dlwrf);
  dlwrf=get_missing_val(lon1,lat1,mask1.*dlwrf,nan,Roa,nan);
%
% get skin temperature [K].
%
  nc=netcdf([NCEP_dir,'tmpsfc_Y',num2str(Y),'M',num2str(M),'.nc']);
  skt=squeeze(nc{'tmpsfc'}(tin,:,:));
  close(nc);
%  skt=get_missing_val(lon1,lat1,mask1.*skt);
  skt=get_missing_val(lon1,lat1,mask1.*skt,nan,Roa,nan);
  skt=skt-273.15; 
%
% computes net longwave heat flux following Dickey et al (1994) 
% (see air_sea tools, fonction lwhf)
%
  nlwf = -lwhf(skt,dlwrf); 
  radlw=interp2(lon1,lat1,nlwf,lon,lat,interp_method);

%
% get the  downward longwave heat flux  
%
  radlw_in=interp2(lon1,lat1,dlwrf,lon,lat,interp_method);

%
% 6: Wind & Wind stress [m/s]
%
nc=netcdf([NCEP_dir,'ugrd10m_Y',num2str(Y),'M',num2str(M),'.nc']);
uwnd=squeeze(nc{'ugrd10m'}(tin,:,:));
close(nc)
uwnd=get_missing_val(lon1,lat1,mask1.*uwnd,nan,Roa,nan);
%uwnd=get_missing_val(lon1,lat1,mask1.*uwnd);
uwnd=interp2(lon1,lat1,uwnd,lon,lat,interp_method);
%
nc=netcdf([NCEP_dir,'vgrd10m_Y',num2str(Y),'M',num2str(M),'.nc']);
vwnd=squeeze(nc{'vgrd10m'}(tin,:,:));
close(nc)
vwnd=get_missing_val(lon1,lat1,mask1.*vwnd,nan,Roa,nan);
%vwnd=get_missing_val(lon1,lat1,mask1.*vwnd);
vwnd=interp2(lon1,lat1,vwnd,lon,lat,interp_method);
%
% Compute the stress
%
wspd=sqrt(uwnd.^2+vwnd.^2);
[Cd,uu]=cdnlp(wspd,10.);
rhoa=air_dens(tair,rhum*100);
tx=Cd.*rhoa.*uwnd.*wspd;
ty=Cd.*rhoa.*vwnd.*wspd;
%
% Rotations on the ROMS grid
%
cosa=cos(angle);
sina=sin(angle);
%
sustr=rho2u_2d(tx.*cosa+ty.*sina);
svstr=rho2v_2d(ty.*cosa-tx.*sina);
%
% uwnd et vwnd sont aux points 'rho'
%
u10=rho2u_2d(uwnd.*cosa+vwnd.*sina);
v10=rho2v_2d(vwnd.*cosa-uwnd.*sina);
%
% Fill the ROMS files
%
if ~isempty(nc_frc)
  nc_frc{'sustr'}(tout,:,:)=sustr;
  nc_frc{'svstr'}(tout,:,:)=svstr;
end
if ~isempty(nc_blk)
  nc_blk{'tair'}(tout,:,:)=tair;
  nc_blk{'rhum'}(tout,:,:)=rhum;
  nc_blk{'prate'}(tout,:,:)=prate;
  nc_blk{'wspd'}(tout,:,:)=wspd;
  nc_blk{'radlw'}(tout,:,:)=radlw;
  nc_blk{'radlw_in'}(tout,:,:)=radlw_in;
  nc_blk{'radsw'}(tout,:,:)=radsw;
  nc_blk{'uwnd'}(tout,:,:)=u10;
  nc_blk{'vwnd'}(tout,:,:)=v10;
%  nc_blk{'sustr'}(tout,:,:)=sustr;
%  nc_blk{'svstr'}(tout,:,:)=svstr;
end

