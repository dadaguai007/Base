function A = CFGaussPulse(a0, t0, freq, c)
% Frequency CGaussion
    A = a0 * sqrt((2 * pi * t0.^2) / (1 + 1i*c)) * exp(-0.5 * ((2 * pi * freq * t0).^2) / (1 + 1i*c));

end
