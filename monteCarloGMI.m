function [GMI, NGMI] = monteCarloGMI(rx, tx, M, constType, px)
    % Monte Carlo based generalized mutual information (GMI) estimation.

    if nargin < 5
        px = [];
    end

    % Constellation parameters
    constSymb = GrayMapping(M, constType);
    % Normalize constellation
    Es = sum(abs(constSymb).^2 .* px);
    constSymb = constSymb / sqrt(Es);


    % Get bit mapping
    b = log2(M);
    bitMap = demodulateGray(constSymb, M, constType);
    bitMap = reshape(bitMap, [], b);

    % Ensure that the signal sequences are disposed in columns
    if size(rx, 2) > size(rx, 1)
        rx = rx.';
    end
    if size(tx, 2) > size(tx, 1)
        tx = tx.';
    end
    % polarization multi
    nModes = size(tx, 2);  % Number of signal modes
    GMI = zeros(1, nModes);
    NGMI = zeros(1, nModes);
    
    % create px
    if isempty(px)
        % If px is not defined, assume uniform distribution
        px = (1 / M) * ones(size(constSymb));
    end

    % Calculate source entropy
    H = -sum(px .* log2(px));

    % Symbol normalization
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

    for k = 1:nModes
        % Set the noise variance
        sigma2 = var(rx(:, k) - tx(:, k));

        % Demodulate transmitted symbol sequence
        btx = demodulateGray(sqrt(Es) * tx(:, k), M, constType);

        % Soft demodulation of the received symbols
        LLRs = calcLLR(rx(:, k), sigma2, constSymb, bitMap, px);

        % LLR clipping
        LLRs(LLRs == inf) = 500;
        LLRs(LLRs == -inf) = -500;

        % Compute bitwise MIs and their sum
        %确保MIperBitPosition数组的长度等于M的二进制表示的位数。这样可以确保对每个比特位都进行计算。
        b = log2(M);
        MIperBitPosition = zeros(1, b);

        for n = 1:b
            MIperBitPosition(n) = H / b - mean(log2(1 + exp((2 * btx(n:b:end) - 1) .* LLRs(n:b:end))));
        end
        GMI(k) = sum(MIperBitPosition);
        NGMI(k) = GMI(k) / H;
    end
end
% 其中，H是信源熵，b是比特数，btx是发送符号序列的每个比特位的子序列，LLRs是接收符号序列的每个比特位的LLR子序列。