function wrf_colorbar(colaxis,loc,vname,fsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function my_colorbar(colaxis,loc,vname,fsize)
%
% put a colorbar at a fixed position
%
% input:
%
%  colpos   position of the colorbar [left, bottom, width, height]
%  colaxis  = [cmin cmax cint]: values assigned to the first and
%           last colors in the current colormap
%           (default: fit the min and max values of the variable)
%  vname    name of the variable (string)
%           (default: '')
%  fsize    font size
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2, 
  loc='v';
  vname='';
  fsize=10;
elseif nargin<3,
  vname='';
  fsize=10;
elseif nargin<4,
  fsize=10;
end;

dc=colaxis(3);

if loc=='v',
  colpos=[.9 .2 .03 .3];
  x=[0:1];
  y=[colaxis(1):dc:colaxis(2)];
  [X,Y]=ndgrid(x,y);
  subplot('position',colpos)
  contourf(X,Y,Y,[colaxis(1):dc:colaxis(2)])
  set(gca,'XTick',x,'XTickLabel',[' '])
%  set(gca,'YAxisLocation','right')
  text(x(length(x)),y(1)-.5, ...
       vname,'FontSize',fsize);
else,
  colpos=[.3 .07 .3 .03];
  y=[0:1];
  x=[colaxis(1):dc:colaxis(2)];
  [X,Y]=ndgrid(x,y);
  subplot('position',colpos)
  contourf(X,Y,X,[colaxis(1):dc:colaxis(2)])
  set(gca,'YTick',y,'YTickLabel',[' '])
  text(x(length(x))+.5,y(1), ...
       vname,'FontSize',fsize);
end;

set(gca,'Layer','top');

return
