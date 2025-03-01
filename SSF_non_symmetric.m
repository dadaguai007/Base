function [E, z] = SSF_non_symmetric(E, hz, Lspan, alpha, gamma, D, Fc, Fs)
c = 299792458;           % 光速 (m/s)
lamba = c / Fc;
Alpha= 1e-3 * alpha / (10 * log10(exp(1)));
beta2 = -(D * lamba^2) / (2 * pi * c);

Nfft = length(E);
w= 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);

z = 0;
while z <= Lspan
    % 线性传输
    E = fft(E);
    E = E .* exp(-Alpha * hz + 1i * (beta2/2) * (w.^2) * hz);

    % 非线性传输
    E = ifft(E);
    E = E .* exp(1i * gamma * (abs(E).^2) * hz);

    z = z + hz;
end

end
