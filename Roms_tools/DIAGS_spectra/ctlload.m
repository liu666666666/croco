function ctl=ctlload(root,his,model,nstrt,nend,ninc);
%%% 
% function ctl=ctlload(root,his,model,nstrt,nend,ninc);
% loads the outputs files for a model (ex: pac50 or uswc15)
% given a file type type (his or avg or something else). 
% root is the directory where these files can be found
% (and root has to exclude the first /)
% nstrt is the starting index for the files, nend the final 
% index and ninc the increment (ie the number of outputs in 
% a given file. 
% to load pac50_his.0000.nc, pac50_his.0010.nc ... til pac50_his.0400.nc  
%  use ctl=ctlload(root,1,pac50,0,400,10)
%  ctl=ctlload('net/inca/capet/CONFIGS/PEC/PEC1_IA/',0,'pec1_ia',100,120,5);
if his==1 || his==-10
 ftype='_his.';
elseif his==0 || his==-11
 ftype='_avg.';
elseif his==-9
 ftype='_diaM.';
else  % if soda
 ftype='_';
end
strg='';
ninc=max(ninc,1);
for ii=nstrt:ninc:nend
  if ii==0
    strg=[strg ' ' '''/' root model ftype 'nc'''];
  elseif ii<10 
    strg=[strg ' ' '''/' root model ftype '000' num2str(ii) '.nc'''];
  elseif ii<100
    strg=[strg ' ' '''/' root model ftype '00' num2str(ii) '.nc'''];
  elseif ii< 1000
    strg=[strg ' ' '''/' root model ftype '0' num2str(ii) '.nc'''];
  else
    strg=[strg ' ' '''/' root model ftype '' num2str(ii) '.nc'''];
  end
end
disp(strg)
if his==1 || his==0 || his==9 
eval(['ctl=rnt_timectl({' strg '}, ''ocean_time'');'])
elseif (his==-9 || his==-10 || his==-11)
eval(['ctl=rnt_timectl({' strg '}, ''scrum_time'');'])
else
eval(['ctl=rnt_timectl({' strg '}, ''TIME1'');'])
end
return
