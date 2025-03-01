function distance = hamming_dist(x, y)

    x = dec2bin(x,2);
    y = dec2bin(y,2);

    distance = sum(x ~= y, 2);

    
end