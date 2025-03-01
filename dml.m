function [Et, Pout]=dml(Pin,V,param)
%The Pin is the power ,The V is the signal of Electrical
if isfield(param, 'f')
    f = param.f;
end
if isfield(param, 'alpha')
    alpha = param.alpha;
end
if isfield(param, 'H')
    H = param.H;
end
if isfield(param, 'Fs')
    Fs = param.Fs;
end

if isfield(param, 'k')
    k = param.k;
end
if isfield(param, 'type')
    type = param.type;
end
% just give the PW(Pin) and the other ，such as f ，phi ，phasenoise add after
% the dml function


Pout = Pin*V;  
if strcmp(type,'on')
Pout = real(ifft(fft(Pout).*ifftshift(H(f))));
end
% adiabatic chirp
% k = 1e12;
%modulation depth.
% m=0.8;

% d = 1-m+m*V;  % modified according to the paper
% Phase = alpha/2*log(d) + alpha/2*k*Pin*cumsum(d)*1/fs;

% Pout = Pout.*exp(1j*alpha/2*log(abs(Pout).^2)); mistake
% Add transient chirp
dphi=alpha/2*log(Pout)+alpha/2*k*cumsum(Pout)*1/Fs;
% dphi=alpha/2*log(Pout)+alpha/2*k*Pout;
% Clip
% Pout = Pout.*exp(1j*dphi);



Pout(Pout < 0) = 0;

Et = sqrt(Pout);
Et = Et.*exp(1j*dphi);


end