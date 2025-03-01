function distance = euclidean_dist_ampm(x, y)

    x = dec2bin(x,3);

    switch x
        case '000'
        distance = norm((1-1i)/sqrt(10)-y);

        case '001'
        distance = norm((-3+3i)/sqrt(10)-y);

        case '010'
        distance = norm((1+3i)/sqrt(10)-y);

        case '011'
        distance = norm((-3-1i)/sqrt(10)-y);
       
        case '100'
        distance = norm((3-3i)/sqrt(10)-y);
       
        case '101'
        distance = norm((-1+1i)/sqrt(10)-y);
        
        case '110'
        distance = norm((3+1i)/sqrt(10)-y);
         
        case '111'
        distance = norm((-1-3i)/sqrt(10)-y);

        otherwise
        warning('Wrong input');
    end

    
end