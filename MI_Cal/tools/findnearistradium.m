function index = findnearistradium(amplitude, r2)
    [~, index] = min(abs(amplitude^2 - r2));
end