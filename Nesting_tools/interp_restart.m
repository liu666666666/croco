function handles=interp_restart(h,handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get everything in order to compute the child restart file
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
%  Copyright (c) 2004-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(handles.parentgrid)
  handles=get_parentrstname(h,handles);
end
if isempty(handles.childgrid)
  handles=get_childgrdname(h,handles); 
end
if isempty(handles.parentrst)
  handles=get_parentrstname(h,handles);
end
lev=str2num(handles.parentrst(end));
if isempty(lev)
  childname=[handles.parentrst,'.1'];
else
  childname=[handles.parentrst(1:end-1),num2str(lev+1)];
end
Answer=questdlg(['Child restart name: ',childname,' OK ?'],'','Yes','Cancel','Yes');
switch Answer
 case {'Cancel'}
  return
 case 'Yes'
  handles.childrst=childname;
end
nested_restart(handles.childgrid,handles.parentrst,handles.childrst, ...
               handles.vertical_correc,handles.extrapmask)
return
