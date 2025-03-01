function [E, z] = SSF_symmetric(E, hz, Lspan, alpha, gamma, D, Fc, Fs)
c = 299792458;           % 光速 (m/s)
lamba = c / Fc;
Alpha = 1e-3 * alpha / (10 * log10(exp(1)));
beta2 = -(D * lamba.^2) / (2 * pi * c);
Nfft = length(E);
%omega
w = 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);
z = 0;

E = fft(E);

while z <= Lspan
    %  linear
    E = E .* exp(-Alpha * (hz/2) + 1i * (beta2/2) * (w.^2) * (hz/2));

    %  non linear
    E = ifft(E);
    E = E .* exp(1i * gamma * (abs(E).^2) * hz);

    % linear
    E = fft(E);
    E = E .* exp(-Alpha * (hz/2) + 1i * (beta2/2) * (w.^2) * (hz/2));

    z = z + hz;
end

E = ifft(E);

end
