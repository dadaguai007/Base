function MI = calcMI(rx, tx, sigma2, constSymb, pX)
    % Mutual information (MI) calculation (circular AGWN channel).

    % rx received signal
    % tx Transmitted symbol signal
    % sigma2 varriance of noise
    % constSymb Constellation symbols. 直接形成的调制码
    % px  prob. mass function (p.m.f.) of the constellation symbols.星座符号的概率质量函数
    N = length(rx);
    H_XgY = 0;
    H_X = sum(-pX .* log2(pX)); %H_x=log2(M);

    for k = 1:N
        %找到最接近当前发送符号的星座符号索引
        [~, indSymb] = min(abs(tx(k) - constSymb));
        %条件熵
        %即给定发送符号条件下接收符号的概率对数形式，根据高斯噪声模型得到
        log2_pYgX = -(1 / sigma2) * abs(rx(k) - tx(k)).^2 * log2(exp(1));  % log2 p(Y|X)
        %根据接收符号和星座符号的概率质量函数，计算联合概率
        %
        pXY = exp(-(1 / sigma2) * abs(rx(k) - constSymb).^2) .* pX;  % p(Y,X) = p(Y|X) * p(X)
        % 根据联合概率计算边缘概率 
        pY = sum(pXY);  % p(Y) = sum(q(Y|X) * p(X)) in X
        % 迭代进行
        H_XgY = H_XgY - (log2_pYgX + log2(pX(indSymb)) - log2(pY));
    end
    H_XgY = H_XgY / N;

    MI = H_X - H_XgY;

end
%pXY代表接收信号取到rx[k]和发送信号取到constSymb的联合概率。
% 它是通过将接收信号与星座点的距离的权重与发送信号的概率质量函数相乘。
% 表示接收信号与星座点的距离的权重。权重越大，越接近星座点。
% pX是发送信号的概率质量函数，表示发送信号取到某个值的概率