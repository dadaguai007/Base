function wout = my_dec(win, alpha, constellation)

wout = win;
Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
% 第一轮误差调整
for i = 2:length(win)-1
    if error_judge(i-1) == error_judge(i+1)
        temp = (error(i-1)+error(i+1))/2;
    else
        temp = 0;
    end
    wout(i) = wout(i) + alpha*temp;
end

% 第二轮误差调整
Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
% alpha = alpha / (1+alpha);
for i = 2:length(Q) - 1
    if error_judge(i-1) == error_judge(i+1)
        if error_judge(i-1) == error_judge(i)
            temp = error(i-1)+error(i+1); % +++/--- 判决有可能发生错误，更正。
        else
            temp = 0; % +-+/-+- 证明判决正确，不更正。
        end
    else
        % temp = 0;
        judge = sign(error(i-1)+error(i+1));
        if error_judge(i) ~= judge
            temp = 0;
        else
            temp = error(i-1)+error(i+1);
        end
    end
    wout(i) = wout(i) + alpha*temp/2;
end

end


% 这段Matlab代码定义了一个名为`my_dec`的函数，该函数实现了一种基于误差反馈的信号解码算法，可能用于通信系统中的信号解调或均衡。下面详细描述该函数的功能：
% 
% 1. **函数定义和输入输出**：
%定义了一个名为`my_dec`的函数，该函数接受三个输入参数：`win`（输入信号向量）、`alpha`（反馈增益参数）和`constellation`（星座图，用于信号判决），并返回一个输出向量`wout`。
% 
% 2. **初始化和初步判决**：
%     - `wout = win;`：将输入信号`win`复制到输出信号`wout`中。
%     - `Q = decision(wout, constellation);`：调用`decision`函数，根据星座图`constellation`对信号`wout`进行判决，得到判决结果`Q`。
%     - `error = wout - Q;`：计算信号`wout`与判决结果`Q`之间的误差`error`。
%     - `error_judge = sign(error);`：获取误差的符号，得到误差判决向量`error_judge`。
% 
% 3. **第一轮误差反馈调整**：
%     - `for i = 2:length(win)-1`：遍历信号向量`win`的中间元素（不包括第一个和最后一个元素）。
%     - `if error_judge(i-1) == error_judge(i+1)`：如果当前元素的前一个和后一个误差符号相同。
%         - `temp = (error(i-1)+error(i+1))/2;`：计算前一个和后一个误差的平均值作为调整值`temp`。
%     - `else`：否则
%         - `temp = 0;`：调整值`temp`为0。
%     - `wout(i) = wout(i) + alpha*temp;`：根据调整值`temp`和反馈增益`alpha`对当前信号元素`wout(i)`进行调整。
% 
% 4. **第二轮判决和误差反馈调整**：
%     - `Q = decision(wout, constellation);`：再次根据星座图对调整后的信号`wout`进行判决，得到新的判决结果`Q`。
%     - `error = wout - Q;`：计算新的误差`error`。
%     - `error_judge = sign(error);`：获取新的误差符号，得到新的误差判决向量`error_judge`。
%     - `for i = 2:length(Q) - 1`：再次遍历信号向量`Q`的中间元素。
%     - `if error_judge(i-1) == error_judge(i+1)`：如果当前元素的前一个和后一个误差符号相同。
%         - `if error_judge(i-1) == error_judge(i)`：如果当前元素的误差符号与前一个和后一个误差符号都相同。
%             - `temp = error(i-1)+error(i+1);`：计算前一个和后一个误差的和作为调整值`temp`。
%         - `else`：否则
%             - `temp = 0;`：调整值`temp`为0。
%     - `else`：如果当前元素的前一个和后一个误差符号不同。
%         - `judge = sign(error(i-1)+error(i+1));`：计算前一个和后一个误差之和的符号`judge`。
%         - `if error_judge(i) ~= judge`：如果当前误差符号与`judge`不同。
%             - `temp = 0;`：调整值`temp`为0。
%         - `else`：否则
%             - `temp = error(i-1)+error(i+1);`：计算前一个和后一个误差的和作为调整值`temp`。
%     - `wout(i) = wout(i) + alpha*temp/2;`：根据调整值`temp`和反馈增益`alpha`对当前信号元素`wout(i)`进行第二次调整。
% 
% 5. **函数返回**：
%     - 函数最终返回经过两轮误差反馈调整后的信号向量`wout`。
% 
% 综上所述，这段代码通过对信号进行判决、计算误差、根据误差符号进行反馈调整，实现了对输入信号的解码或均衡处理，以改善信号的质量或准确性。需要注意的是，代码中使用了未定义的`decision`函数，实际应用中需要根据具体的星座图和判决规则来实现该函数。 