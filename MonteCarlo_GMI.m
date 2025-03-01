function [GMI, NGMI] = MonteCarlo_GMI(rx, tx, M, constType, px)
% Monte Carlo based generalized mutual information (GMI) estimation.

if nargin < 5
    %如果没有提供先验符号概率 (px)，则将其设为空数组。
    px = [];
end

% Constellation parameters
%得到灰色映射的调制符号。
%     constSymb = GrayMapping(M, constType);
if strcmp(constType,'qam')
    z = (0:M-1)';
    constSymb=qammod(z,M);

    % Get bit mapping
    %计算每个调制符号的比特映射。
    b = log2(M);

    bitMap=qamdemod(constSymb, M,'OutputType','bit');
    bitMap = reshape(bitMap,b, []);
    bitMap=bitMap.';
elseif strcmp(constType,'pam')
    z = (0:M-1)';
    constSymb=pammod(z,M);
    % Get bit mapping
    %计算每个调制符号的比特映射。
    b = log2(M);

    bitMap=pamdemod(constSymb, M);
    % function [res,Res] = val2arr(val, bits)
    % 将整数转换为二进制，当然转换后的二进制结果需要进行一个前后颠倒
    % 第一位位于索引1处
    Res=zeros(length(bitMap),log2(M));
    for j=1:length(bitMap)
        res = zeros(1, b);
        for i = 1:b
            res(i) = bitand(bitshift(bitMap(j), -(i-1)), 1);
        end
        Res(j,:)=flip(res);
    end
    bitMap = Res;
end

% Ensure symbol sequences are in column form
%将接收和发送的符号序列调整为列向量的形式
if size(rx, 2) > size(rx, 1)
    rx = rx.';
end
if size(tx, 2) > size(tx, 1)
    tx = tx.';
end

nModes = size(tx, 2);  % Number of signal modes
GMI = zeros(1, nModes);
NGMI = zeros(1, nModes);
%如果未提供先验符号概率 (px)，则假定均匀分布。
if isempty(px)
    px = 1/M * ones(size(constSymb));
end

% Normalize constellation
%对调制符号进行归一化。
Es = sum(abs(constSymb).^2 .* px);
constSymb = constSymb / sqrt(Es);

% Calculate source entropy
%     计算源熵。
H = sum(-px .* log2(px));

% Symbol normalization
for k = 1:nModes
    if  strcmp(constType,'qam')
        % Correct (possible) phase ambiguity
        %解决了可能存在的相位模糊性
        rot = mean(tx(:, k) ./ rx(:, k));
        rx(:, k) = rot * rx(:, k);
    end
    % Symbol normalization
    rx(:, k) = pnorm(rx(:, k));
    tx(:, k) = pnorm(tx(:, k));
end

for k = 1:nModes
    % Set the noise variance
    % 噪声方程
    sigma2 = var(rx(:, k) - tx(:, k));

    % Demodulate transmitted symbol sequence
    %解调发送符号序列
    if strcmp(constType,'qam')
        btx=qamdemod(sqrt(Es) * tx(:, k), M,'OutputType','bit');
        btx=btx.';
    elseif strcmp(constType,'pam')
        btx=pamdemod(sqrt(Es) * tx(:, k), M);
        Res=zeros(length(btx),log2(M));
        for j=1:length(btx)
            res = zeros(1, b);
            for i = 1:b
                res(i) = bitand(bitshift(btx(j), -(i-1)), 1);
            end
            Res(j,:)=flip(res);
        end
        btx = Res;
    end
    % Soft demodulation of the received symbols
    %计算接收符号序列 rx 的 LLRs ，软解调
    LLRs = calcLLR(rx(:, k), sigma2, constSymb, bitMap, px);

    % LLR clipping
    LLRs(LLRs == inf) = 500;
    LLRs(LLRs == -inf) = -500;

    % Compute bitwise MIs and their sum
    %计算比特位置上的互信息
    MIperBitPosition = zeros(1, b);

    for n = 1:b
        MIperBitPosition(n) = H/b - mean(log2(1 + exp((2 * btx(n:b:end) - 1) .* LLRs(n:b:end))));
    end
    % 计算所有比特位置上的互信息并求和，得到广义互信息（GMI），然后通过除以源熵 H 得到归一化互信息（NGMI）。
    GMI(k) = sum(MIperBitPosition);
    NGMI(k) = GMI(k) / H;
end
end
