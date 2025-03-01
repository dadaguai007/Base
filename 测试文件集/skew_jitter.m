%应用skew 和 jitter

%skew
% 长度是否等于1，如果等于1，表示只提供了一个偏斜值；把单一的skew复制为一个包含2倍输入信号维度数目的向量
%  长度不等于1，检查它是否等于2倍输入信号的维度数目。如果不等于2倍输入信号的维度数目，将在 skew 后面添加零值，以使其长度等于2倍输入信号的维度数目。
% Computing Skew
if length(skew) == 1
    obj.skew = ones(1,2*size(input,2))*obj.skew;
elseif length(skew) ~= 2*size(input,2)
    obj.skew(2*size(input,2)) = 0;
end

%jitter
% Computing Jitter
%随机值乘上jitter的方差，生成 jitter； 时钟偏差应用到每一个信号上。
%randn(floor(size(input,1)/obj.resamplingRate)-1,1)，随机信号生成。
%(in.Nss/obj.resamplingRate)*sqrt(obj.jitterVariance)，时钟抖动的标准差。Nss为每符号采样次数
%随机数序列再加上 (1+obj.clockError)*相应的长度，对应时钟偏差的部分，用于引入时钟偏差。
timing = cumsum([0 ; (in.Nss/obj.resamplingRate)*sqrt(obj.jitterVariance)* ...
    randn(floor(size(input,1)/obj.resamplingRate)-1,1)+(1+obj.clockError)* ...
    ones(floor(size(input,1)/obj.resamplingRate)-1,1)]);
%将时钟偏差赋值到每一列
timing = repmat(timing, 1, 2*size(input,2));

% Computing Timing
for ii = 1:2*size(input,2)
    %之前的偏差加上 Nss*obj.skew(ii)/obj.resamplingRate：来自偏斜参数的时间偏移。
    %20/obj.resamplingRate：一个常数时间偏移。
    timing(:,ii) = timing(:,ii) + in.Nss*obj.skew(ii)/obj.resamplingRate + 20/obj.resamplingRate;
end

%将时钟偏差应用到信号中：
% Skew, Jitter and Clock Error Insertion
output(:,ii) = interp1(1:size(real(signal),1), real(signal), timing(:,2*ii-1), 'spline') + ...
    + 1j*interp1(1:size(imag(signal),1), imag(signal), timing(:,2*ii), 'spline');
