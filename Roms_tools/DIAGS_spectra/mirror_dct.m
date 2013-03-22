function out=mirror_dct(inp)
% this function takes a field and performs 2 successive 
% symmetries in order to come up with a bi periodic field. 
siz=size(inp);
if length(siz)==2
[M,L]=size(inp);
Dm=2*M;Dl=2*L;
out=zeros(Dm,Dl);
out(1:M,1:L)=inp;                            % fill lower left quadrant
out(M+1:Dm,1:L)=flipdim(inp,1);              % fill lower right quadrant
out(1:M,L+1:Dl)=flipdim(inp,2);              % fill upper left quadrant
out(M+1:Dm,L+1:Dl)=flipdim(flipdim(inp,1),2);% fill upper right quadrant
elseif length(siz)==3
[M,L,K]=size(inp);
Dm=2*M;Dl=2*L;
out=zeros(Dm,Dl,K);
out(1:M,1:L,:)=inp;                            % fill lower left quadrant
out(M+1:Dm,1:L,:)=flipdim(inp,1);              % fill lower right quadrant
out(1:M,L+1:Dl,:)=flipdim(inp,2);              % fill upper left quadrant
out(M+1:Dm,L+1:Dl,:)=flipdim(flipdim(inp,1),2);% fill upper right quadrant
else 
disp([' SIZE OF INP ARG NOT OK '])
end
return
