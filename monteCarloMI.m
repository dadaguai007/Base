function MI = monteCarloMI(rx, tx, M, constType, px)
    % Monte Carlo based mutual information (MI) estimation.

    if nargin < 5
        px = [];
    end
    % create the px
    if isempty(px)
        px = (1 / M) * ones(M, 1);  % Assume uniform distribution
    end

    % Constellation parameters

    constSymb = GrayMapping(M, constType);
    Es = sum(abs(constSymb).^2 .* px);
    constSymb = constSymb / sqrt(Es);

    % Ensure that the signal sequences are disposed in columns
    if size(rx, 2) > size(rx, 1)
        rx = rx.';
    end
    if size(tx, 2) > size(tx, 1)
        tx = tx.';
    end
% polarization multi
    nModes = size(rx, 2);
    MI = zeros(1, nModes);

    for k = 1:nModes
        if ismember(constType, {'qam', 'psk'})
            % Correct possible phase ambiguity
            rot = mean(tx(:, k) ./ rx(:, k));
            rx(:, k) = rot * rx(:, k);
        end

        % Symbol normalization
        rx(:, k) = pnorm(rx(:, k));
        tx(:, k) = pnorm(tx(:, k));
    end

    % Estimate noise variance from the data
    noiseVar = var(rx - tx);

    for k = 1:nModes
        sigma2 = noiseVar(k);
        MI(k) = calcMI(rx(:, k), tx(:, k), sigma2, constSymb, px);
    end
end



% 计算接收信号和发送信号的比值可以得到相位值的原因是：
% 相位是描述信号在时域上的旋转程度的物理量，而比值可以反映出两个信号之间的相对相位关系。
% 具体来说，假设接收信号和发送信号分别为rx和tx，它们的比值为tx / rx。
% 当接收信号和发送信号相位一致时，比值为1；
% 当接收信号相位领先发送信号半个周期时，比值为e^(-jpi/2)=-j；
% 当接收信号相位滞后发送信号半个周期时，比值为e^(-jpi)=-1；
% 当接收信号相位领先发送信号一个周期时，比值为e^(-j2pi)=1。
% 
% 因此，计算接收信号和发送信号的比值可以得到一个复数，其幅角即为相位值。
% 通过计算一组接收信号和发送信号的比值的平均值，可以得到一个旋转因子，
% 用于纠正接收信号和发送信号之间的相位模糊。