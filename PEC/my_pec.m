function wout = my_pec(win, alpha, constellation)
% 示例：% 对于FFE的增强噪声进行消除
%     pec_out = my_pec(ffe_out, alpha, constellation);
Q = decision(win, constellation);
error = win - Q;
error_judge = sign(error);
wout = win;

for i = 2:length(Q)-1
    if error_judge(i-1) == error_judge(i+1)
        if error_judge(i-1) == error_judge(i)
            temp = error(i-1)+error(i+1); % +++/--- 判决有可能发生错误，更正。
        else
            temp = 0; % +-+/-+- 证明判决正确，不更正。
        end
    else
        temp = 0;
    end
    wout(i) = wout(i) + alpha*temp/2;
end

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
for i = 2:length(Q) - 1
    judge = sign(error(i-1)+error(i+1));
    if error_judge(i) ~= judge
        temp = 0;
    else
        temp = error(i-1)+error(i+1);
    end
    wout(i) = wout(i) + alpha*temp/2;
end

end



% 实现了一种基于误差反馈的信号处理算法，可能用于通信系统中的信号判决和纠错。
% 
% ### 函数定义和输入输出

% - 输入参数：
%     - `win`：输入信号向量。
%     - `alpha`：一个控制反馈调整强度的参数。
%     - `constellation`：星座图，用于对信号进行判决。
% - 输出参数：
%     - `wout`：经过处理后的输出信号向量。
% 
% ### 初始化和初步判决
% ```matlab
% Q = decision(win, constellation);
% error = win - Q;
% error_judge = sign(error);
% wout = win;
% ```
% - `Q = decision(win, constellation);`：调用`decision`函数，根据给定的星座图`constellation`对输入信号`win`进行判决，得到判决结果`Q`。这里`decision`函数应是自定义的，用于将信号映射到星座点上。
% - `error = win - Q;`：计算输入信号`win`与判决结果`Q`之间的误差。
% - `error_judge = sign(error);`：获取误差的符号，生成误差判决向量`error_judge`，用于后续判断误差的方向。
% - `wout = win;`：将输出信号`wout`初始化为输入信号`win`。
% 
% ### 第一轮误差反馈调整
% ```matlab
% for i = 2:length(Q)-1
%     if error_judge(i-1) == error_judge(i+1)
%         if error_judge(i-1) == error_judge(i)
%             temp = error(i-1)+error(i+1);
%         else
%             temp = 0;
%         end
%     else
%         temp = 0;
%     end
%     wout(i) = wout(i) + alpha*temp/2;
% end
% ```
% - 此循环遍历判决结果向量`Q`中除第一个和最后一个元素之外的所有元素。
% - 如果当前元素的前一个和后一个误差符号相同（`error_judge(i-1) == error_judge(i+1)`）：
%     - 并且当前元素的误差符号与前一个误差符号也相同（`error_judge(i-1) == error_judge(i)`），则计算前一个和后一个误差的和作为调整值`temp`。
%     - 否则，将调整值`temp`设为0。
% - 如果当前元素的前一个和后一个误差符号不同，则将调整值`temp`设为0。
% - 最后，根据调整值`temp`和参数`alpha`对当前信号元素`wout(i)`进行调整，调整量为`alpha*temp/2`。
% 
% ### 第二轮判决和误差反馈调整
% ```matlab
% Q = decision(wout, constellation);
% error = wout - Q;
% error_judge = sign(error);
% for i = 2:length(Q) - 1
%     judge = sign(error(i-1)+error(i+1));
%     if error_judge(i) ~= judge
%         temp = 0;
%     else
%         temp = error(i-1)+error(i+1);
%     end
%     wout(i) = wout(i) + alpha*temp/2;
% end
% ```
% - 再次调用`decision`函数，根据星座图对经过第一轮调整后的信号`wout`进行判决，得到新的判决结果`Q`。
% - 计算新的误差`error`和误差符号`error_judge`。
% - 此循环再次遍历判决结果向量`Q`中除第一个和最后一个元素之外的所有元素。
% - 计算当前元素的前一个和后一个误差之和的符号`judge`。
% - 如果当前误差符号`error_judge(i)`与`judge`不同，则将调整值`temp`设为0。
% - 否则，计算前一个和后一个误差的和作为调整值`temp`。
% - 最后，根据调整值`temp`和参数`alpha`对当前信号元素`wout(i)`进行第二次调整，调整量为`alpha*temp/2`。
% 
% ```
% 函数返回经过两轮误差反馈调整后的输出信号向量`wout`。
% 
% 总体而言，这段代码通过对信号进行判决、计算误差、根据误差符号进行反馈调整，尝试优化信号的判决结果，以提高信号的准确性或改善信号质量。需要注意的是，代码中依赖的`decision`函数并未给出具体实现，实际应用中需要根据具体的星座图和信号处理需求来实现该函数。  