function p = maxwellBoltzmann(lamba, const)
    p = zeros(size(const));

    for ind = 1:length(const)
        p(ind) = exp(-lamba * abs(const(ind)).^2);
    end

    p = p ./ sum(p);
end
