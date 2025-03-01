function Hf = Himdd(f,L, D, alpha, type)
%% Fiber frequency response for an IM-DD system with transient chirp dominant
% This transfer function is for optical power not electic field
% i.e., Hfiber(f) = Pout(f)/Pin(f).
% Inputs:
% - f: frequency vector (Hz)
% - wavelength: wavelength (m)
% - alpha (optional, default alpha = 0): chirp parameter with
% sign convention such that for DML alpha > 0
% - type (optional, default type = 'small signal')
% in m), and alpha (optional, default zero) (chirp paramter).

if not(exist('alpha', 'var')) % if chirp parameter is not defined
    alpha = 0;
end
c = 299792458;
c_kms = c/ 1e3 ;   % speed of light (vacuum) in km/s
Fc = 193.1e12;
lamba = c_kms / Fc ;
% lamba=1.550e-9;
% Alpha =alpha;
beta = -(D * lamba.^2) / (2 * pi * c_kms);

% beta = beta2(wavelength);
theta = -1/2*beta*(2*pi*f).^2*L; % theta = -1/2*beta2*w.^2*L

if  strcmpi(type, 'large signal')
    %% Large signal 大信号模型
    mIM = 0.7; % modulation index is set to 70%
    Dphi = pi/2; % i.e., transient chirp dominant
    mFM = alpha/2*mIM;
    u = 2*mFM*sin(theta);
    Hf = cos(theta).*(besselj(0, u) - besselj(2, u)*exp(1j*Dphi)) - 2*exp(1j*Dphi)/(1j*mIM)*besselj(1, u);  % fiber large-signal frequency response
elseif strcmpi(type, 'small signal')
    %% Small signal 小信号模型
    Hf = cos(theta) - alpha*sin(theta);  % fiber small-signal frequency response
else
    error('fiber/Hf: undefined type of fiber frequency response')
end
end


function h = hdisp(L, t, lambda)
%% Dispersion impulse response inverseFourier(Hdisp(f))
% Inputs:
% - f = frequency vector (Hz)
% - lambda = wavelength (m)

b = 1/2*beta2(lambda)*L;
h = sqrt(-pi*1j/b)/(2*pi)*exp(1j*t.^2/(4*b));
% h = sqrt(1/(4*pi*b))*(cos(t.^2/(4*b) - pi/4) + 1j*sin(t.^2/(4*b) - pi/4));
end