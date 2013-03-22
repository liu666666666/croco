lims=complims(lon,lat,dx)

dx0=dx*0.5e-5;
Lmin=max(find(lon(:,1)<lonw+dx0 & lon(:,1)>lonw-dx0));
Lmax=min(find(lon(:,1)<lone+dx0 & lon(:,1)>lone-dx0));
Mmin=max(find(lat(1,:)<lats+dx0 & lat(1,:)>lats-dx0));
Mmax=min(find(lat(1,:)<latn+dx0 & lat(1,:)>latn-dx0));
lims=[Lmin Lmax Mmin Mmax]

return
end
