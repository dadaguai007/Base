function xout=Duob(SpS, Nsamples, reverse, alpha)
% Dubo FTN格式
if nargin<1
    SpS=2;
end
if nargin<2
    Nsamples=16;
end
if nargin<3
    reverse='false';
end
if nargin<4
    alpha=0.01;
end

N=Nsamples+2*SpS;
% p = pulseShape('rrc', SpS, N, alpha);
p = rcosdesign(alpha,N,SpS,'sqrt');

% 向左滚动SpS 个符号
x = p + circshift(p, -SpS);

%翻转
if strcmp(reverse,'true')
    indrev = length(p):-1:1;
    x = x(indrev);
end

xout=x(SpS+1:end-SpS);


end