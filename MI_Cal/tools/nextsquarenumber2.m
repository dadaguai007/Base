function out = nextsquarenumber2(in)
    out = ceil(sqrt(in))^2;
    if mod(out, 2)
        slog('The modulation order is incorrect, and therefore, the correct constellation diagram cannot be generated.', 'ERR');
    end
end