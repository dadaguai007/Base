function A = FGaussPulse(a0, t0, freq)
% Transformada Fourier Gaussiano
    A = a0 * t0 * sqrt(2 * pi) * exp(-0.5 * (2 * pi * freq * t0).^2);
end
