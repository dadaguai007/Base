function A = SuperGaussPulse(a0, t0, t, m)
%Super Gaussion
    A = a0 * exp(-0.5 * (t/t0).^(2*m));
end
