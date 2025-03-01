function vblock = gen_vol_block(block, order)
L = size(block, 1);
vblock = [];
switch order
    case 1
        vblock = block;
    case 2
        for i = 1:L
            for j = 1:i
                vblock = [vblock; block(i, :).*block(j, :)];
            end
        end
    case 3
        for i = 1:L
            for j = 1:i
                for p = 1:j
                    vblock = [vblock; block(i, :).*block(j, :).*block(p, :)];
                end
            end
        end
    case 4
        for i = 1:L
            for j = 1:i
                for p = 1:j
                    for q = 1:p
                        vblock = [vblock; block(i, :).*block(j, :).*block(p, :).*block(q, :)];
                    end
                end
            end
        end
    otherwise
        error('仅支持最高四阶运算');
end
end

