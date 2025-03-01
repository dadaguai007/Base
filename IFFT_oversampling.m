function [xt, time] = IFFT_oversampling(X,N,L)

%频域过采样
if nargin<3  
    L=1;  
end
NL=N*L; 
T=1/NL; 
time = [0:T:1-T];  
X = X(:).';
xt = L*ifft([X(1:N/2)  zeros(1,NL-N)  X(N/2+1:end)], NL);