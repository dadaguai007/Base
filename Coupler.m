function [sigRxo_aerfa_remaining,sigRxo_aerfa]=Coupler(aerfa,sigTxo_1,sigTxo_2)

if ~isrow(sigTxo_1)
    % 是否为行向量
    sigTxo_1=sigTxo_1.';
end

if nargin<3
    % 一路输入置零
    sigTxo_2=zeros(size(sigTxo_1));
end

if ~isrow(sigTxo_2)
    % 是否为行向量
    sigTxo_2=sigTxo_2.';
end

% 测试输入向量是否为行向量
% aerfa为分光比，一般选取较小路径输出的比例；1:9的分光比，aerfa应为0.1；

% Coulper矩阵
Coupler=[sqrt(1-aerfa),1j*sqrt(aerfa);1j*sqrt(aerfa),sqrt(1-aerfa)];

% 输入信号矩阵
sigTxo_Array=[sigTxo_1;sigTxo_2];

% 分光
sigRxo_Array=Coupler*sigTxo_Array;


sigRxo_aerfa_remaining=sigRxo_Array(1,:);
sigRxo_aerfa=sigRxo_Array(2,:);


end
