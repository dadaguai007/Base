function [wout, error, w] = my_ffe_eda_lms(win, ntaps, step, ref, sps, iter, alpha)
%实现了一种基于最小均方（LMS）算法的前馈均衡器（FFE）与误差驱动自适应（EDA）相结合的信号处理算法。
w = zeros(ntaps, 1);
w(round(ntaps/2)) = 1;
L = floor((length(win)-ntaps)/sps+1);
wout = zeros(L, 1);
trainlen = length(ref);
error = zeros(iter*trainlen, 1);
tempIdx = 1:ntaps;
edatemp = zeros(length(alpha), 1);
% training
for i = 1:iter
    for j = 1:trainlen
        idx = tempIdx + (j-1)*sps;%计算当前输入信号子块的索引
        block = win(idx);
        temp = w.' * block - alpha*edatemp;
        error((i-1)*trainlen+j) = ref(j) - temp;
        edatemp = [error((i-1)*trainlen+j); edatemp(1:end-1)];
        w = w + step * error((i-1)*trainlen+j) * conj(block);
    end
end
% equalizing
for i = 1:L
    idx = tempIdx + (i-1)*sps;
    block = win(idx);
    wout(i) = w.' * block;
end
end