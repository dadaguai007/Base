function wout = my_pnc_enhanced(win, alpha, constellation)

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

Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);
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
    % judge = sign(error(i-1)+error(i+1));
    % if error_judge(i) ~= judge
    %     temp = 0;
    % else
    %     temp = error(i-1)+error(i+1);
    % end
    wout(i) = wout(i) + alpha*temp/2;
end

end





% ### 第一轮误差反馈调整
% ```matlab
% for i = 2:length(Q)-1
%     if error_judge(i-1) == error_judge(i+1)
%         if error_judge(i-1) == error_judge(i)
%             temp = error(i-1)+error(i+1); % +++/--- 判决有可能发生错误，更正。
%         else
%             temp = 0; % +-+/-+- 证明判决正确，不更正。
%         end
%     else
%         % temp = 0;
%         judge = sign(error(i-1)+error(i+1));
%         if error_judge(i) ~= judge
%             temp = 0;
%         else
%             temp = error(i-1)+error(i+1);
%         end
%     end
%     wout(i) = wout(i) + alpha*temp/2;
% end
% ```
% - 此循环遍历判决结果向量`Q`中除第一个和最后一个元素之外的所有元素。
% - 如果当前元素的前一个和后一个误差符号相同（`error_judge(i-1) == error_judge(i+1)`）：
%     - 并且当前元素的误差符号与前一个误差符号也相同（`error_judge(i-1) == error_judge(i)`），则计算前一个和后一个误差的和作为调整值`temp`。这是因为连续相同符号的误差可能表示判决有误，需要进行更正。
%     - 否则，说明误差的变化模式可能是正确的，将调整值`temp`设为0。
% - 如果当前元素的前一个和后一个误差符号不同：
%     - 计算前一个和后一个误差之和的符号`judge`。
%     - 如果当前误差符号`error_judge(i)`与`judge`不同，说明当前误差可能是正常的，将调整值`temp`设为0。
%     - 否则，计算前一个和后一个误差的和作为调整值`temp`。
% - 最后，根据调整值`temp`和参数`alpha`对当前信号元素`wout(i)`进行调整，调整量为`alpha*temp/2`。
% 
% ### 第二轮误差反馈调整
% - 再次调用`decision`函数，根据星座图对经过第一轮调整后的信号`wout`进行判决，得到新的判决结果`Q`。
% - 计算新的误差`error`和误差符号`error_judge`。
% - 与第一轮调整类似，再次遍历判决结果向量`Q`中除第一个和最后一个元素之外的所有元素，进行相同的误差判断和调整操作。
% 

% 函数最终返回经过两轮误差反馈调整后的输出信号向量`wout`。
% 
% 综上所述，该函数通过两轮基于误差符号和误差和的符号判断的反馈调整，对输入信号进行优化处理，以提高信号判决的准确性，改善信号质量。需要注意的是，代码中依赖的`decision`函数需要根据具体的星座图和信号处理需求进行自定义实现。  