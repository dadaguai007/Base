function out = decision(in, constellation)
    out = zeros(size(in));
    for idx = 1:size(in, 2)
        distance = abs(in(:, idx) - constellation(:).');
        [~, index] = min(distance');
        out(:, idx) = constellation(index);
    end
end