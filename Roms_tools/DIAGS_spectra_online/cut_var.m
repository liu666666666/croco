
function var=cut_var(var,lims)

Lmin=lims(1); Lmax=lims(2); Mmin=lims(3); Mmax=lims(4);
if length(size(var))==3,
  var=squeeze(var(Lmin:Lmax,Mmin:Mmax,:));
elseif length(size(var))==2,
  var=squeeze(var(Lmin:Lmax,Mmin:Mmax));
else
  disp('wrong var size in CUT_VAR')
end

return
