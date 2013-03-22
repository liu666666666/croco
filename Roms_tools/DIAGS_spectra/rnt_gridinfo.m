% (R)oms (N)umerical (T)oolbox
% 
% FUNCTION grdinfo = rnt_gridinfo(gridid)
%
% Loads the grid configuration for gridid
% To add new grid please edit this file.
% just copy an existing one and modify for
% your needs. It is simple.
%
% If you editing this file after using
% the Grid-pak scripts use the content 
% of variable "nameit" for gridid.
%
% Example: CalCOFI application
%
%    grdinfo = rnt_gridinfo('calc')
%
% RNT - E. Di Lorenzo (edl@ucsd.edu)

function gridindo=rnt_gridinfo(gridid)

% initialize to defaults
       gridindo.id      = gridid;
       gridindo.name    = '';
       gridindo.grdfile = '';	 	 
	gridindo.N       = 20;
       gridindo.thetas  = 5;  
       gridindo.thetab  = 0.4;  	 	 
       gridindo.tcline  = 0;
	 gridindo.cstfile = '/d6/edl/ROMS-pak/Grid-pak/WorldCstLine.mat';

    root='/work/CONFIGS/';
    switch gridid

%   case 'jet'
%        gridindo.id      = gridid;
%        gridindo.name    = 'baroclinic jet grid 10km';
%        gridindo.grdfile = 'jet_grd.nc';
%        gridindo.N       = 10;
%        gridindo.thetas  = 1;
%        gridindo.thetab  = 0;
%        gridindo.hc      = 400;
%        gridindo.Method  = 1;
%   
        
  case 'jet'
       gridindo.id      = gridid;
       gridindo.name    = 'baroclinic jet grid 5km';
       gridindo.grdfile = 'jet_grd.nc';
       gridindo.N       = 30;
       gridindo.thetas  = 5
       gridindo.thetab  = 0;
       gridindo.hc      = 400;
       gridindo.Method  = 1;
  
  case 'jet10'
       gridindo.id      = gridid;
       gridindo.name    = 'baroclinic jet grid 10km';
       gridindo.grdfile = 'jet10_grd.nc';
       gridindo.N       = 10;
       gridindo.thetas  = 1;
       gridindo.thetab  = 0;
       gridindo.hc      = 400;
       gridindo.Method  = 1;

  case 'jet5'
       gridindo.id      = gridid;
       gridindo.name    = 'baroclinic jet grid 5km';
       gridindo.grdfile = 'jet5_grd.nc';
       gridindo.N       = 10;
       gridindo.thetas  = 1;
       gridindo.thetab  = 0;
       gridindo.hc      = 400;
       gridindo.Method  = 1;

  case 'line90'
       gridindo.id      = gridid;
       gridindo.name    = 'CalCOFI obs. grid';
       gridindo.grdfile = which('Line90Grid.nc');
       gridindo.N       = 20;
       gridindo.thetas  = 7;
       gridindo.thetab  = 0;
       gridindo.tcline  = 200;
       gridindo.cstfile = which('rnt_calccoast.mat');

  case 'calc'
       gridindo.id      = gridid;
       gridindo.name    = 'CalCOFI grid - Manu';
       gridindo.grdfile = which('grid-calcofi.nc');
       gridindo.N       = 20;
       gridindo.thetas  = 5;
       gridindo.thetab  = 0.4;
       gridindo.tcline  = 200;
       gridindo.cstfile = which('rnt_calccoast.mat');

  case 'ee12'
       gridindo.id      = gridid;
       gridindo.name    = 'ee12';
       gridindo.grdfile = '/celtic/EDDY_FLUX/EE12/ee12_grd.nc';
       gridindo.N       = 40;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'ee10'
       gridindo.id      = gridid;
       gridindo.name    = 'ee10';
       gridindo.grdfile = '/celtic/EDDY_FLUX/EE10/ee10_grd.nc';
       gridindo.N       = 40;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'ee6'
       gridindo.id      = gridid;
       gridindo.name    = 'ee6';
       gridindo.grdfile = [root '/EDDY_FLUX/ee6/ee6_grd.nc']; 
       gridindo.N       = 40;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'ee3'
       gridindo.id      = gridid;
       gridindo.name    = 'ee3';
       gridindo.grdfile = '/celtic/EDDY_FLUX/EE3/ee3_grd.nc';
       gridindo.N       = 40;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'ee1'
       gridindo.id      = gridid;
       gridindo.name    = 'ee1';
       gridindo.grdfile = [root '/EDDY_FLUX/ee1/ee1_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'ee1p'
       gridindo.id      = gridid;
       gridindo.name    = 'ee1p';
       gridindo.grdfile = [root '/EDDY_FLUX/EE1/ee1_grd.nc'];
       gridindo.N       = 60;
       gridindo.hc      = 10;
       gridindo.thetas  = 7.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';


  case 'ee0'
       gridindo.id      = gridid;
       gridindo.name    = 'ee0';
       gridindo.grdfile = [root 'EDDY_FLUX/ee0/ee0_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'ee00'
       gridindo.id      = gridid;
       gridindo.name    = 'ee00';
       gridindo.grdfile = [root 'EDDY_FLUX/ee00/ee00_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'sa3'
       gridindo.id      = gridid;
       gridindo.name    = 'SA3';
       gridindo.grdfile = '/celtic/CONFIGS/Brazil/SA/sa3_grd.nc';
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';


  case 'brc1'
       gridindo.id      = gridid;
       gridindo.name    = 'BRC1';
       gridindo.grdfile = [root 'CONFIGS/Brazil/BRC1/brc1_grd.nc'];
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'brc2'
       gridindo.id      = gridid;
       gridindo.name    = 'BRC2';
       gridindo.grdfile = [root 'CONFIGS/Brazil/BRC2/brc2_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 100;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 2;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'brc3'
       gridindo.id      = gridid;
       gridindo.name    = 'BRC3';
       gridindo.grdfile = [root 'Brazil/BRC3/brc3_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'brc3red'
       gridindo.id      = gridid;
       gridindo.name    = 'BRC3';
       gridindo.grdfile = [root 'Brazil/BRC3/brc3red_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';


  case 'brc4'
       gridindo.id      = gridid;
       gridindo.name    = 'BRC4';
       gridindo.grdfile = [root 'CONFIGS/Brazil/BRC4/brc4_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';


  case 'mc2'
       gridindo.id      = gridid;
       gridindo.name    = 'MC2';
       gridindo.grdfile = [root 'Brazil/MC2/mc2_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'mc05'
       gridindo.id      = gridid;
       gridindo.name    = 'MC05';
       gridindo.grdfile = ['/home/capet/CONFIGS/Brazil/MC05/mc05_grd.nc'];
       gridindo.N       = 24;
       gridindo.hc      = 15;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';


  case 'mc2red'
       gridindo.id      = gridid;
       gridindo.name    = 'MC2';
       gridindo.grdfile = [root 'Brazil/MC2/mc2red_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'ob1'
       gridindo.id      = gridid;
       gridindo.name    = 'OB1';
       gridindo.grdfile = [root '/CONFIGS/Brazil/OB1/ob1_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'ob1h4000'
       gridindo.id      = gridid;
       gridindo.name    = 'OB1';
       gridindo.grdfile = [root 'CONFIGS/Brazil/OB1/ob1_h4000_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'mc0'
       gridindo.id      = gridid;
       gridindo.name    = 'MC0';
       gridindo.grdfile = [root 'CONFIGS/Brazil/MC0/mc0_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'mc01'
       gridindo.id      = gridid;
       gridindo.name    = 'MC01';
       gridindo.grdfile = [root 'CONFIGS/Brazil/MC01/mc01_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'cfrio5'
       gridindo.id      = gridid;
       gridindo.name    = 'CFRIO5';
       gridindo.grdfile = [root 'CONFIGS/Brazil/CFRIO/cfrio5_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';


  case 'cfrio3'
       gridindo.id      = gridid;
       gridindo.name    = 'CFRIO3';
       gridindo.grdfile = [root 'CONFIGS/Brazil/CFRIO/cfrio3_grd.nc'];
       gridindo.N       = 40;
       gridindo.hc      = 30;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorld.mat';

  case 'usw2'
       gridindo.id      = gridid;
       gridindo.name    = 'USW2';
       gridindo.grdfile = '/work/CONFIGS/USW2/usw2_grd.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'usw50'
       gridindo.id      = gridid;
       gridindo.name    = 'USW50';
       gridindo.grdfile = '/celtic/CONFIGS/USW50/usw50_grd.nc';
       gridindo.grdfileH= '/celtic/CONFIGS/USW50/usw50_grdH.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';


  case 'usw51'
       gridindo.id      = gridid;
       gridindo.name    = 'USW5 with new topo';
       gridindo.root    = [root 'usw51/'];
       gridindo.grdfile = [gridindo.root 'usw51_grd.nc'];
       gridindo.grdfileH= [gridindo.root 'usw51_grdH.nc'];
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'usw51H'
       gridindo.id      = gridid;
       gridindo.name    = 'USW5 with new topo';
       gridindo.grdfile = '/celtic/CONFIGS/USW5.1/usw51_grdH.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'usw52'
       gridindo.id      = gridid;
       gridindo.name    = 'BIG USW51 ';
       gridindo.grdfile = '/home/capet/celtic/USW52/usw52_grd.nc';
       gridindo.grdfileH= '/home/capet/celtic/USW52/usw52_grdH.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'beng1'
       gridindo.id      = gridid;
       gridindo.name    = 'Benguela Pierrick -Child 5km';
       gridindo.grdfile = '/home/capet/celtic/BENG_PIERRICK/beng1_grd.nc';
       gridindo.grdfileH= '/home/capet/celtic/BENG_PIERRICK/beng1_grdH.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'beng1H'
       gridindo.id      = gridid;
       gridindo.name    = 'Benguela Pierrick -Child 5km';
       gridindo.grdfile = '/home/capet/celtic/BENG_PIERRICK/beng1_grdH.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'cana1'
       gridindo.id      = gridid;
       gridindo.name    = 'Benguela Pierrick -Child 5km';
       gridindo.grdfile = '/home/capet/celtic/can1/cana1_grd.nc.1';
       gridindo.grdfileH= '/home/capet/celtic/can1/cana1_grdH.nc.1';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'cana1H'
       gridindo.id      = gridid;
       gridindo.name    = 'Benguela Pierrick -Child 5km';
       gridindo.grdfile = '/home/capet/celtic/can1/_grdH.nc.1';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'es2'
       gridindo.id      = gridid;
       gridindo.name    = 'Earth Simulator';
       gridindo.grdfile = '/celtic/CONFIGS/ES/PK/s1_2km_100l_fpos_525d_t5_diags/t1/history.out';
       gridindo.N       = 100;
       gridindo.hc      = 100;
       gridindo.thetas  = 5.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'es1'
       gridindo.id      = gridid;
       gridindo.name    = 'Earth Simulator';
       gridindo.grdfile = '/work/CONFIGS/ES/s1_1km_200l_strat0/es1_his.0001.nc';
       gridindo.N       = 200;
       gridindo.hc      = 100;
       gridindo.thetas  = 5.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'es4'
       gridindo.id      = gridid;
       gridindo.name    = 'Earth Simulator';
       gridindo.grdfile = '/work/CONFIGS/ES/ES4/history.out_Phillips';
       gridindo.N       = 50;
       gridindo.hc      = 100;
       gridindo.thetas  = 5.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'es8'
       gridindo.id      = gridid;
       gridindo.name    = 'Earth Simulator';
       gridindo.grdfile = '/work/CONFIGS/ES/ES8/es8_SCharney_his.nc';
       gridindo.N       = 40;
       gridindo.hc      = 100;
       gridindo.thetas  = 5.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'romsl3'
       gridindo.id      = gridid;
       gridindo.name    = 'CPacific';
       gridindo.grdfile = '/data2/matlab/FILES_season1/roms_grd.nc.3';
       gridindo.N       = 30;
       gridindo.hc      = 5;
       gridindo.thetas  = 6.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/data2/coastline_l.mat';

  case 'romsl2'
       gridindo.id      = gridid;
       gridindo.name    = 'CPacific';
       gridindo.grdfile = '/data2/matlab/FILES_season1/roms_grd.nc.2';
       gridindo.N       = 30;
       gridindo.hc      = 5;
       gridindo.thetas  = 6.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/data2/coastline_l.mat';

  case 'romsl1'
       gridindo.id      = gridid;
       gridindo.name    = 'CPacific';
       gridindo.grdfile = '/data2/matlab/FILES_season1/roms_grd.nc.1';
       gridindo.N       = 30;
       gridindo.hc      = 5;
       gridindo.thetas  = 6.;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/data2/coastline_l.mat';

  case 'pec50'
       gridindo.id      = gridid;
       gridindo.name    = 'USW50';
       gridindo.grdfile = '/celtic/CONFIGS/USW50/usw50_grd.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'uswx15'
       gridindo.id      = gridid;
       gridindo.name    = 'USWX15';
       gridindo.grdfile = '/celtic/CONFIGS/USWX15/uswx15_grd.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'usw17'
       gridindo.id      = gridid;
       gridindo.name    = 'USW17 grid patrick 5km';
       gridindo.grdfile = '/celtic/CONFIGS/USW17/usw17_grd.nc';
       gridindo.N       = 20;
       gridindo.thetas  = 7.0;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'usw15'
       gridindo.id      = gridid;
       gridindo.name    = 'USW15';
       gridindo.grdfile = '/celtic/CONFIGS/USW15/usw15_grd.nc';
       gridindo.N       = 32;
       gridindo.hc      = 50;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'usw8.1'
       gridindo.id      = gridid;
       gridindo.name    = 'USW8 child 1';
       gridindo.grdfile = '/celtic/CONFIGS/USW8/usw8_grd.nc.1';
       gridindo.N       = 32;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

  case 'gc'
       gridindo.id      = gridid;
       gridindo.name    = 'Gulf of California Simona';
       gridindo.grdfile = '/celtic/CONFIGS/GC/gc_grd.nc';
       gridindo.N       = 32;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0.0;
       gridindo.Method  = 1;
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

    case 'pac50'
       gridindo.id      = gridid;
       gridindo.name    = 'Pacific';
       gridindo.grdfile = '/celtic/PAC/Pac50/pacific_jpl_grd_xa.nc';
       gridindo.N       = 30;
       gridindo.hc      = 50;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0.0
       gridindo.Method  = 1;
       gridindo.cstfile ='/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

    case 'sacw'
       gridindo.id      = gridid;
       gridindo.name    = 'sacw';
       gridindo.grdfile = '/net/paracas/paracas/francois_work12/SACW/Interann/POP_ERS_COADS/sacw_grd.nc';
       gridindo.N       = 30;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0.0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile ='/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

    case 'pac25'
       gridindo.id      = gridid;
       gridindo.name    = 'Pacific';
       gridindo.grdfile = '/celtic/Pac25/pac25_grd.nc';
       gridindo.N       = 40;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0.0;
       gridindo.Method  = 1;
       gridindo.cstfile ='/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'uswc_z20'
       gridindo.id      = gridid;
       gridindo.name    = 'USWC 20 levs ';
       gridindo.grdfile = '/celtic/CONFIGS/USWC20/uswc_z20_grd.nc';
       gridindo.N       = 20;
       gridindo.hc      = 50;
       gridindo.thetas  = 7;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'uswc_z20.1'
       gridindo.id      = gridid;
       gridindo.name    = 'USWC 20 levs child 1';
       gridindo.grdfile = 'celtic/CONFIGS/USWC20/uswc_z20_grd.nc.1';
       gridindo.N       = 20;
       gridindo.hc      = 50;
       gridindo.thetas  = 7;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'uswc_z20.2'
       gridindo.id      = gridid;
       gridindo.name    = 'USWC 20 levs child 2';
       gridindo.grdfile = '/celtic/CONFIGS/USWC20/uswc_z20_grd.nc.2';
       gridindo.N       = 20;
       gridindo.hc      = 50;
       gridindo.thetas  = 7;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'usw15_z32'
       gridindo.id      = gridid;
       gridindo.name    = 'USW15 32 levs ';
       gridindo.grdfile = '/celtic/CONFIGS/grid15_z32/usw15_z32_grd.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'usw15_z32.1'
       gridindo.id      = gridid;
       gridindo.name    = 'USW15 32 levs child 1';
       gridindo.grdfile = '/celtic/CONFIGS/grid15_z32/usw15_z32_grd.nc.1';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';


   case 'usw15_z32.2'
       gridindo.id      = gridid;
       gridindo.name    = 'USW15 32 levs child 2';
       gridindo.grdfile = '/celtic/CONFIGS/grid15_z32/usw15_z32_grd.nc.2';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'usw15_strcst'
       gridindo.id      = gridid;
       gridindo.name    = 'USW15 20 levs  20km resolution straight coast, flat bottom';
       gridindo.grdfile = '/???????/usw15_grd.nc';
       gridindo.N       = 20;
       gridindo.thetas  = 7.0;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'peru'
       gridindo.id      = gridid;
       gridindo.name    = 'GRID PERU 10km (P.Penven)';
       gridindo.grdfile = '/celtic/CONFIGS/Peru/roms10_grd.nc';
       gridindo.N       = 32;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'pec1'
       gridindo.id      = gridid;
       gridindo.name    = 'GRID North East Pacific - (Xa) - 7km';
       gridindo.grdfile = '/celtic/CONFIGS/PEC/pec1_grd.nc';
       gridindo.grdfileH= '/celtic/CONFIGS/PEC/pec1_grdH.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'pec2'
       gridindo.id      = gridid;
       gridindo.name    = 'GRID North East Pacific - (Xa) - 35km';
       gridindo.grdfile = '/celtic/CONFIGS/PEC/pec2_grd.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'hfsa3'
       gridindo.id      = gridid;
       gridindo.name    = 'GRID PERU CHILE FRANCOIS';
       gridindo.grdfile = '/net/paracas/paracas/capet/hfsa3_grd.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'hfsa4'
       gridindo.id      = gridid;
       gridindo.name    = 'GRID PERU CHILE FRANCOIS';
       gridindo.grdfile = '/home/capet/celtic/hfsa4/hfsa4_grd.nc';
       gridindo.grdfileH= '/home/capet/celtic/hfsa4/hfsa4_grdH.nc';
       gridindo.N       = 32;
       gridindo.hc      = 10;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.Slon    = 'W';
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

   case 'nep'
       gridindo.id      = gridid;
       gridindo.name    = 'GRID North East Pacific - (Haidstrom)';
       gridindo.grdfile = '/celtic/CONFIGS/Alaska/NEP3_grid.nc';
       gridindo.N       = 40;
       gridindo.thetas  = 6.5;
       gridindo.thetab  = 0;
       gridindo.Method  = 1;
       gridindo.cstfile = '/home/capet/Roms/Matlab/Grid-pak/CoastlineWorldPacific.mat';

    case 'scb'
       gridindo.id      = gridid;
       gridindo.name    = 'Charles''s SCB';
       gridindo.grdfile = '/celtic/CONFIGS/Charles/scb_iswake_grid.nc';
       gridindo.N       = 40;
       gridindo.thetas  = 5;
       gridindo.thetab  = 0;
       gridindo.tcline  = 10;
       gridindo.Method  = 1;
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

    case 'as15'
       gridindo.id      = gridid;
       gridindo.name    = 'Arabian sea 15km';
       gridindo.grdfile = '/work/CONFIGS/Arabian_Sea/as15_grd.nc';
       gridindo.N       = 40;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.hc      = 30;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';
    case 'wio15'
       gridindo.id      = gridid;
       gridindo.name    = 'western indian ocean 15km';
       gridindo.grdfile = '/work/CONFIGS/Arabian_Sea/WIO15/wio15_grd.nc';
       gridindo.N       = 40;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.hc      = 30;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';
    case 'as3'
       gridindo.id      = gridid;
       gridindo.name    = 'arabian sea 3km';
       gridindo.grdfile = '/work/CONFIGS/Arabian_Sea/AS3/as3_grd.nc';
       gridindo.N       = 40;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.hc      = 30;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';
    case 'as1'
       gridindo.id      = gridid;
       gridindo.name    = 'arabian sea 3km';
       gridindo.grdfile = '/work/CONFIGS/Arabian_Sea/AS3/as3_grd.nc';
       gridindo.N       = 80;
       gridindo.thetas  = 6;
       gridindo.thetab  = 0;
       gridindo.hc      = 30;
       gridindo.Method  = 1;
       gridindo.Slon    = 'E';
       gridindo.cstfile = '/home/capet/matlab/Grid-pak/CoastlineWorldPacific.mat';

    otherwise
       gridindo.id      = gridid;
       gridindo.name    = 'null';
       gridindo.grdfile = '/dev/null';
       gridindo.N       = 0;
       gridindo.thetas  = 0;  
       gridindo.thetab  = 0;  	 	 
       gridindo.tcline  = 0;
	 disp([' RNT_GRIDINFO - ',gridid,' not configured']);
    end	 

