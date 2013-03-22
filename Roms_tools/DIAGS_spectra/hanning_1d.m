function A=hanning_1d(A,N);

for k=1:N
  A(2:end-1)=0.5*A(2:end-1)+0.25*(A(1:end-2)+A(3:end));
end

return

