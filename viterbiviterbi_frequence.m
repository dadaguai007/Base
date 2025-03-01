function [Eout, theta] = viterbiviterbi_frequence(E,M,N)

Eout=zeros(size(Ei));
nModes = size(Ei, 2);
zeroPad = zeros(N, nModes);
x = [zeroPad; E; zeroPad];  % Pad the signal with zeros
D=zeros(size(x));
theta=zeros(size(E));
for model = 1:nModes
for index = 2:length(x)
    d=x(index,model)*conj(x(index-1,model));
    D(index)=d.^M;

    if index >= 2*N
        % sum row
        sumD = mean(D);
        % get the best theata
        theta(k-2 * N+1,model) = angle(sumD)/4;
    end
end
Eout(:,model)=E(:,mode).* exp(-1i * theta(:,model));
end