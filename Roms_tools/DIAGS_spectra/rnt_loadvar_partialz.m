% (R)oms (N)umerical (T)oolbox
%
% FUNCTION [fieldout]=rnt_loadvar_partialz(ctl,ind,field,smin,smax)
%
% Extract the variable FIELD for the indices specified
% in the array IND from a composite netcdf file/variable which is defined 
% in the time controll arrays CTL.
% the variable must be 3d and will be extracted only for levels 
% comprised between smin and smax
%
%  example [fieldout]=rnt_loadvar(ctl,[38],'zeta');
%  load index 38 and variabe 'zeta'
%
% SEE rnt_timectl.m to generate a the CTL struct. array 
%      (it is easy ..no worry ok!)
%
% INPUT: example
%     ctl = rnt_timectl(files,timevar);
%     ind = [ 1:6 ] (get time indiceis 1:6 from composite field
%     field = 'temp'
% 
% OUTPUT:
%   fieldout(x,y,z,1:length(ind)) = theCompositeField(x,y,z,ind) 
%
% RNT - E. Di Lorenzo (edl@ucsd.edu)



function [tmp1]=rnt_loadvar_partialz(ctl,ind,field,smin,smax);


  tmp1=0;
  % find info about the variable and initialize arrays
    nc=netcdf(ctl.file{1});
    [s]=size(nc{field}(:,smin:smax,:,:));
    len=length(s);
%   if len==2
%     disp(['bricolo because pb between AKt and AKs']);
%     [s]=ncsize(nc{'AKt'}(:,smin:smax,:,:));
%     if len==2
%       [s]=ncsize(nc{'AKs'}(:,smin:smax,:,:));
%     end
%     len=length(s);
%   end
    lt=length(nc('time'));
    close(nc)
    if lt==1
      s=[length(ind) s];
    else
      s(1)=length(ind);
    end
    tmp1=zeros(s);

  if isempty(ind)
    disp(['rnt_loadvarsum - no time index match for ',field]);
    [s] = size(tmp1); order = [length(s) :-1:1];
    tmp1=permute(tmp1,order);
    tmp1(:)=NaN;
    return
  end
  
  j=0;
  % load array
  for istep=1:length(ctl.segm)-1
    in = find ( ind > ctl.segm(istep) & ind <= ctl.segm(istep+1));
    in_extr = ctl.ind(ind(in));
    
    if ~isempty(in_extr)
      jstart=j(end)+1;
      jend=j(end)+length(in_extr);
      nc=netcdf(ctl.file{istep});
      tmp2=nc{field}(in_extr,smin:smax,:,:) ;
      close(nc)
%size(tmp2)
%size(tmp1)
%size(squeeze(tmp2))
	
	if len==4 
          if lt==1
            if smin==smax
                tmp1(jstart:jend,:,:,:)=tmp2(:,:);     
            else
                 tmp1(jstart:jend,:,:,:)=tmp2(:,:,:);
            end
          else
             tmp1(jstart:jend,:,:,:)=tmp2(:,:,:,:); 
          end
	elseif len==3 
          if lt==1
            tmp1(jstart:jend,:,:,:)=tmp2(:,:,:);
          elseif lt~=1
            tmp1(jstart:jend,:,:)=squeeze(tmp2(:,:,:));
          end
        else
          tmp1(:,:)=squeeze(tmp2(:,:));
        end
      j=jend;
    end
    
  end
  
  tmpmean=tmp1/length(ind);
  [s] = size(tmp1); order = [length(s) :-1:1];
  tmp1=permute(tmp1,order);
  
  return
