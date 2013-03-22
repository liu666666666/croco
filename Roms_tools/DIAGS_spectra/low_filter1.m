function var=low_filter1(ctlavg,vname,it,kchoice,nhan,navg,zr,pm,pn);

%%%% loads and filters a field from a ctl series.
%%% this is done for a sigma level kchoice

%%% zr,pm,pn only if need to compute w from omega
switch vname
case {'hbls','zeta'}
  var=mean(rnt_loadvar(ctlavg,it-navg:it+navg,vname),3);
case 'rho'
  %vname1='temp';tmp1=squeeze(mean(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname1,kchoice-1,kchoice),4));
  %vname1='salt';tmp2=squeeze(mean(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname1,kchoice-1,kchoice),4));
  vname1='temp';tmp1=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname1,kchoice,kchoice)),3);
  vname1='salt';tmp2=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname1,kchoice,kchoice)),3);
  var=rnt_rho_potential_nomex(tmp1,tmp2);
case 'omega2w'
  vname1='u';tmp1=squeeze(mean(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname1,kchoice-1,kchoice),4));
  vname1='v';tmp2=squeeze(mean(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname1,kchoice-1,kchoice),4));
  vname1='omega';var=squeeze(mean(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname1,kchoice-1,kchoice+1),4));
  var=omega2w(var,tmp1,tmp2,pm,pn,zr);var=squeeze(var(:,:,2));
otherwise
  var=mean(squeeze(rnt_loadvar_partialz(ctlavg,it-navg:it+navg,vname,kchoice,kchoice)),3);
end
if nhan>0
  for ihan=1:nhan
    var=hanning(var);
  end
end

return
