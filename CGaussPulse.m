function A = CGaussPulse(a0, t0, t, c)
% chirp gaussion
    A = a0 * exp(-0.5 * (1 + 1i*c) * (t/t0).^2);
end


