% (R)oms (N)umerical (T)oolbox
%
% FUNCTION [ctl]=rnt_timectl(files,timevar);
%
% Constructs a structure array called 'ctl' used to access
% variables whos time content is stored in more than
% one file.
%
% Example: you have 30 history files of the model
% and you want to load the temperature and you want to
% construct a timeindex which is the composite of all
% the 30 files for the temperature.
%
% ctl.time(:)   time value as a concatenated array for all 30 file
% ctl.file{:}   file names used in the composite
% ctl.ind(:)    indicies of the time array
% ctl.segm(:)   indicies that link  ctl.time   to the actual file
%               names. This is usefull for when you use rnt_loadvar
%               accessing multiple files.
%Example:
%
%ctl =
%
%    time: [504x1 double]
%    file: {1x14 cell}
%     ind: [504x1 double]
%    segm: [1 36 72 108 144 180 216 252 288 324 360 396 432 468 504]
%   dates: [504x7 double] (this is the ouput of rnt_date(ctl.time,'s')
%
% INPUT
%   files = {'file1.nc' 'file2.nc' .... }
%   timevar = 'scrum_time' or whatever the name of the time variable
%
% RNT - E. Di Lorenzo (edl@ucsd.edu)

function [ctl]=rnt_timectl(files,varargin)
  
%==========================================================
% % generate time control struct array ctl
%==========================================================
  
  if ~isa(files, 'cell')
    files{1} = files;
  end
  if nargin ==1
     timevar='scrum_time';
  end
  if nargin > 1
     timevar = varargin{1};
  end
  if nargin > 2
     offset=varargin{2};
   end     
  
  
  tmp=0;
  ctl.time=[]; ctl.file=[]; ctl.ind=[]; ctl.segm=0;
  
  for i=1:length(files)
    ctl.file{i} = files{i};
    nc=netcdf(files{i}); tmp=nc{timevar}(:);
    ctl.time = [ctl.time ; tmp ];
    tmp(:)=i; ctl.ind = [ctl.ind ; [1:length(tmp)]'];
%    ctl.segm = [ctl.segm length(tmp)*i];
    ctl.segm = [ctl.segm ctl.segm(end)+length(tmp)];
    close(nc);
  end  
% if strcmp(timevar,'scrum_time')
% ctl.date=rnt_date(ctl.time,'s');
% elseif strcmp(timevar,'ocean_time')
 ctl.date=rnt_date(ctl.time,'s');
% else
% ctl.date=rnt_date(ctl.time,'r');
% end
  
  if nargin > 2
  ctl.date(:,3) = ctl.date(:,3)+offset(1);
  ctl.date(:,2) = ctl.date(:,2)+offset(2);
  ctl.date(:,1) = ctl.date(:,1)+offset(3);
  end
      
  ctl.datenum=datenum(ctl.date(:,3),ctl.date(:,2),ctl.date(:,1));
  ctl.month=ctl.date(:,2);
  ctl.year=ctl.date(:,3);
  
return  
