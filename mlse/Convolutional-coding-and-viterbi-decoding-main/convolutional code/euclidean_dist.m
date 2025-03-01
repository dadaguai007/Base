function distance = euclidean_dist(x, y)

    x = dec2bin(x,2);

    switch x
        case '00'
        distance = norm((1+1i)/sqrt(2)-y);

        case '01'
        distance = norm((1-1i)/sqrt(2)-y);

        case '10'
        distance = norm((-1+1i)/sqrt(2)-y);

        case '11'
        distance = norm((-1-1i)/sqrt(2)-y);

        otherwise
        warning('Wrong input');
    end

    
end