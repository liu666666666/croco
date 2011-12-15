function download_SODA2(Ymin,Ymax,Mmin,Mmax,lonmin,lonmax,latmin,latmax,...
                       OGCM_dir,OGCM_prefix,url,Yorig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extract a subgrid from SODA to get a ROMS forcing
% Store that into monthly files.
% Take care of the Greenwitch Meridian.
% 
%
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    6-Sep-2006 by Pierrick Penven
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
%close all
if nargin < 1
  Ymin=2000;
  Ymax=2000;
  Yorig=2000;
  Mmin=1;
  Mmax=1;
  lonmin=16;
  lonmax=19;
  latmin=-34;
  latmax=-32;
  OGCM_dir='DATA/SODA/';
  OGCM_prefix='SODA2_';
  url='http://soda.tamu.edu/opendap/SODA_2.1.6/MONTHLY/SODA_2.1.6_';
end
%
disp([' '])
disp(['Get data from Y',num2str(Ymin),'M',num2str(Mmin),...
      ' to Y',num2str(Ymax),'M',num2str(Mmax)])
disp(['Minimum Longitude: ',num2str(lonmin)])
disp(['Maximum Longitude: ',num2str(lonmax)])
disp(['Minimum Latitude: ',num2str(latmin)])
disp(['Maximum Latitude: ',num2str(latmax)])
disp([' '])
%
% Create the directory
%
disp(['Making output data directory ',OGCM_dir])
eval(['!mkdir ',OGCM_dir])
%
% Start 
%
disp(['Process the dataset: ',url])
%
% Find a subset of the SODA grid
%
if Mmin<10;
  path=[url,num2str(Ymin),'0',num2str(Mmin),'.cdf'];
else
  path=[url,num2str(Ymin),num2str(Mmin),'.cdf'];
end
[i1min,i1max,i2min,i2max,i3min,i3max,jrange,krange,lon,lat,depth]=...
 get_SODA2_subgrid(path,lonmin,lonmax,latmin,latmax);
%
% Loop on the years
%
for Y=Ymin:Ymax
  disp(['Processing year: ',num2str(Y)])
%
% Loop on the months
%
  if Y==Ymin
    mo_min=Mmin;
  else
    mo_min=1;
  end
  if Y==Ymax
    mo_max=Mmax;
  else
    mo_max=12;
  end
  for M=mo_min:mo_max
    disp(['  Processing month: ',num2str(M)])
%
% Get the time for this year and month
%
    trange='';
    SODA_time=datenum(Y,M,15)-datenum(Yorig,1,1);
%
    if M<10;
      path=[url,num2str(Y),'0',num2str(M),'.cdf'];
    else
      path=[url,num2str(Y),num2str(M),'.cdf'];
    end
%
% Extract SODA data
%
    extract_SODA(OGCM_dir,OGCM_prefix,path,Y,M,...
                 lon,lat,depth,SODA_time,...
                 trange,krange,jrange,...
                 i1min,i1max,i2min,i2max,i3min,i3max,...
                 Yorig)
  end
end
return