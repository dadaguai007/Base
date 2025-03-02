function printprogress(k, iteration)
% total bar:50.
    progress = floor(k / iteration * 100);
    if ~mod(progress, 5)
        barlength = floor(progress * 50 / 100);
        equalbar = repmat('=', 1, barlength);
        fprintf(1,'%s%d%%%', equalbar, progress);
        fprintf(1, '\n');
    end
    if k == iteration
        fprintf(1, '=========================Over=========================\n');
    end
end
