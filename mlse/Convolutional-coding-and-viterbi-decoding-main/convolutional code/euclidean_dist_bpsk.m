function distance = euclidean_dist_bpsk(x, y)

    x = dec2bin(x,2);
    y = y';

    switch x
        case '00'
        distance = (abs(y(1)-1))^2+(abs(y(2)-1))^2;

        case '01'
        distance = (abs(y(1)-1))^2+(abs(y(2)+1))^2;

        case '10'
        distance = (abs(y(1)+1))^2+(abs(y(2)-1))^2;

        case '11'
        distance = (abs(y(1)+1))^2+(abs(y(2)+1))^2;

        otherwise
        warning('Wrong input');
    end

    
end