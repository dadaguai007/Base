function out = binmap(in)
    in = num2str(in);
    in = in(:, 1:3:end);
    out = bin2dec(in);
end