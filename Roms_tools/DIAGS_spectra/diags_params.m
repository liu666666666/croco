%========================================================================
%
%             PARAMETER FILE FOR SPECTRAL DIAGNOSTICS
%
% P. Marchesiello 2010
%========================================================================

model0='jet';
root=['/data/models/JET/5K_ISO/'];
datadir='/data/models/JET/5K_ISO/';
his=1;
avg=1;

lims_lev  = [1 60 1 100];  % horiz grid index range
dx        = 5;            % horiz grid scale (km)
lit       = [50:60];        % time grid index range
klist     = [9:9];         % vert grid index range
kchoice   = 9;             % vert grid index of choice

Nlev          =  1;        % nb nested levels
Nrec          =  0;        % nb record per files (0: 1 file only)
lastfile_indx =  0;        %                     (0: 1 file only)

window    = 2;             % windowing method
filtamp   = 1;             % filter spectra for plotting
filtcff   = 3;             % number of hanning filter steps
supersamp = 1;             % super-sampling of k

g=9.8;rho0=1025;

dx_lev=dx*3.^(Nlev-(1:Nlev))*1.e3;        % dx for nested grids 

cff_scale = 1.e4;                          % scaling factor
xmin=0.9e-5; xmax=6.e-4; ymin=-1; ymax=1;  % limits for plots
lonw=-3; lone=+3; lats=-5; latn=+5;
istr=2;                                    % index of first k plotted

specname=[datadir,model0,'_spectrum'];

dirout_EPS = [root,'EPS'];
dirout_JPG = [root,'JPG'];
dirout_PDF = [root,'PDF'];

