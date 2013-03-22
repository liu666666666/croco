function z = zlevs(h,zeta,theta_s,theta_b,hc,N,Method,type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%  pierrick 2002                                                            %
%                                                                           %
% function z = zlevs(h,zeta,theta_s,theta_b,hc,N,Method,type);              %
%                                                                           %
% this function compute the depth of rho or w points for ROMS               %
%                                                                           %
% On Input:                                                                 %
%                                                                           %
%    type    'r': rho point 'w': w point                                    %
%    Method   1 = old  -- 2 = new sasha                                     %
%                                                                           %
% On Output:                                                                %
%                                                                           %
%    z       Depths (m) of RHO- or W-points (3D matrix).                    %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(['COMPUTING ZLEV WITH thetas=' num2str(theta_s)])
%disp(['                        hc=' num2str(hc)])
%disp(['                    Method=' num2str(Method)])
[M,L]=size(h);
%
% Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
%
cff1=1./sinh(theta_s);
cff2=0.5/tanh(0.5*theta_s);
if type=='w'
  sc=((0:N)-N)/N;
  N=N+1;
else
  sc=((1:N)-N-0.5)/N;
end
Cs=(1.-theta_b)*cff1*sinh(theta_s*sc)...
    +theta_b*(cff2*tanh(theta_s*(sc+0.5))-0.5);
%
% Create S-coordinate system: based on model topography h(i,j),
% fast-time-averaged free-surface field and vertical coordinate
% transformation metrics compute evolving depths of of the three-
% dimensional model grid.
%    
z=zeros(N,M,L); z0=zeros(N,M,L);
if Method==1
  for k=1:N
    z0(k,:,:)=hc.*sc(k)+(h-hc).*Cs(k);
  end

elseif Method==2
  for k=1:N
    z0(k,:,:)=h.*(hc.*sc(k)+h.*Cs(k))./(h+hc);
  end
end

% twick so that we dont need to squeeze z0. Otherwise 
% it fails when the second dimension is also 1 (ie 
% when doing bry files for example). 
zetap=zeros(1,M,L); zetap(1,:,:)=zeta; zeta=zetap;
hp=zeros(1,M,L); hp(1,:,:)=h; h=hp;ih=1./h;
for k=1:N
    z(k,:,:)=z0(k,:,:)+zeta.*(1+z0(k,:,:).*ih);
end

%if type=='w'
%  hmin=min(min(h));
%  hmax=max(max(h));
%  for k=N:-1:1
%    cff1=sc(k)*hc+(hmin-hc)*Cs(k);
%    cff2=sc(k)*hc+(0.5*(hmin+hmax)-hc)*Cs(k);
%    cff3=sc(k)*hc+(hmax-hc)*Cs(k);
%    disp([num2str(k,6),' | ',num2str(sc(k),6),' | ',num2str(Cs(k)),' | ',...
%         num2str(cff1),' | ',num2str(cff2),' | ',num2str(cff3)])
%  end
%end

return

