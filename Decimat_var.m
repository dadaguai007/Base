function [out, idx, symbols] = Decimat_var(x, Nss_in, Nss_out, offset)
%使用最大方差进行降采样(抽样)，找最佳采样点
if ~isvector(x)
    warning('Input signal should be a vector.');
end
if nargin<3 || isempty(Nss_out), Nss_out = 1; end

%要求输入采样点的数量必须是输出采样点数量的整数倍
%             r不是整数倍，需要报错
r = Nss_in/Nss_out;

%输入信号长度为Nss_in 的整数倍
N = Nss_in*fix(numel(x)/Nss_in);
%FIXME Two lines below: quick fix for out of memory errors
x = x(1:N);

%计算限制信号长度 LIM,取较小的那个值
SYMBOLS_LIMIT = 1e5;
%             LIM = Nss_in*fix(SYMBOLS_LIMIT/Nss_in);
LIM = SYMBOLS_LIMIT*Nss_in;
LIM = min(LIM,N);
%将信号按照每符号多少个采样点进行整形
%每一列表示一个符号（symbol）
symbols = reshape(x(1:LIM),Nss_in,[]).'; % reshape signal into columns (column=symbol)

%计算每个符号的方差，找到具有最大方差的采样点。这是用于选择降采样点的标准
[~,ptr] = max(var(symbols)); % find maximum variance point
ptr = mod(ptr-1+offset,r)+1;
%将r当做步长，也说的过去
idx = ptr:r:N;
out = x(idx);


end