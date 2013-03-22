function out=mirror_dct_1d(inp)
% this function takes a field and performs 1 
% symmetry in the xsi direction to get a periodic 
% field (assuming the field is already periodic 
% in eta direction). 
siz=size(inp);
if length(siz)==2
  [M,L]=size(inp);
  Dm=2*M;Dl=L;
  out=zeros(Dm,Dl);
  out(1:M,1:L)=inp;                            % fill lower left quadrant
  out(M+1:Dm,1:L)=flipdim(inp,1);              % fill lower right quadrant
elseif length(siz)==3
  [M,L,K]=size(inp);
  Dm=2*M;Dl=L;
  out=zeros(Dm,Dl,K);
  out(1:M,1:L,:)=inp;                            % fill lower left quadrant
  out(M+1:Dm,1:L,:)=flipdim(inp,1);              % fill lower right quadrant
else 
  disp([' SIZE OF INP ARG NOT OK '])
end
return
