function A = GaussPulse(a0, t0, t)
% GaussPulse
    A = a0 * exp(-0.5 * (t./t0).^2);
end