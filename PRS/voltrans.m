function out = voltrans(in, order)
out = [];
len = length(in);
switch order
    case 2
        for i1 = 1:len
            for i2 = i1:len
                out = [out in(i1)*in(i2)];
            end
        end
    case 3
        for i1 = 1:len
            for i2 = i1:len
                for i3 = i2:len
                    out = [out in(i1)*in(i2)*in(i3)];
                end
            end
        end
    otherwise
        error('Invalid input order, 2 or 3 is accepted');
end
end

