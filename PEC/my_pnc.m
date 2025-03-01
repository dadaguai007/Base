function wout = my_pnc(win, alpha, constellation, stage)
% pnc_out = my_pnc(ffe_out, alpha, constellation);
if nargin < 4
    stage = 2;
end
wout = win;
Q = decision(wout, constellation);
error = wout - Q;
error_judge = sign(error);

for i = 2:length(win)-1
    if error_judge(i-1) == error_judge(i+1)
        temp = (error(i-1)+error(i+1))/2;
    else
        temp = 0;
    end
    wout(i) = wout(i) + alpha*temp;
end

if stage >= 2
    Q = decision(wout, constellation);
    error = wout - Q;
    error_judge = sign(error);
    
    for i = 2:length(win)-1
        if error_judge(i-1) == error_judge(i+1)
            temp = error(i-1)+error(i+1);
            % temp = (error(i-1)+error(i+1))/2^2;
        else
            if abs(error(i-1)) >= abs(error(i+1))
                temp = error(i+1) / 2;
            else
                temp = error(i-1) / 2;
            end
        end
        wout(i) = wout(i) + alpha*temp;
    end 
end

if stage == 3
    Q = decision(wout, constellation);
    error = wout - Q;
    error_judge = sign(error);
    
    for i = 2:length(win)-1
        if error_judge(i-1) == error_judge(i+1)
            temp = error(i-1)+error(i+1);
            % temp = 0;
        else
            if abs(error(i-1)) >= abs(error(i+1))
                temp = error(i+1) / 2;
            else
                temp = error(i-1) / 2;
            end
        end
        wout(i) = wout(i) + alpha*temp;
    end
end

end


%函数实现了一种基于误差反馈的信号处理算法，可能用于通信系统中的信号判决和校正。
% 
% ### 函数定义与输入输出
% ```matlab
% function wout = my_pnc(win, alpha, constellation, stage)
% ```
% - 输入参数：
%     - `win`：输入信号向量。
%     - `alpha`：一个控制反馈调整强度的参数。
%     - `constellation`：星座图，用于对信号进行判决。
%     - `stage`：一个可选参数，用于指定处理的阶段，默认值为2。
% 
% - 输出参数：
%     - `wout`：经过处理后的输出信号向量。
% 
% 
% ### 初始化与初步处理
% ```matlab
% wout = win;
% Q = decision(wout, constellation);
% error = wout - Q;
% error_judge = sign(error);
% ```
% - 将输入信号`win`赋值给输出信号`wout`。
% - 调用`decision`函数（该函数需自定义实现），根据星座图`constellation`对`wout`进行判决，得到判决结果`Q`。
% - 计算信号`wout`与判决结果`Q`之间的误差`error`。
% - 获取误差的符号，得到误差判决向量`error_judge`。
% 
% ### 第一轮误差反馈调整
% ```matlab
% for i = 2:length(win)-1
%     if error_judge(i-1) == error_judge(i+1)
%         temp = (error(i-1)+error(i+1))/2;
%     else
%         temp = 0;
%     end
%     wout(i) = wout(i) + alpha*temp;
% ```
% - 遍历输入信号`win`中除第一个和最后一个元素之外的所有元素。
% - 如果当前元素的前一个和后一个误差符号相同，则计算前一个和后一个误差的平均值作为调整值`temp`。
% - 否则，将调整值`temp`设为0。
% - 根据调整值`temp`和参数`alpha`对当前信号元素`wout(i)`进行调整。
% 
% ### 第二轮误差反馈调整（`stage >= 2`）
% ```matlab
% if stage >= 2
%     Q = decision(wout, constellation);
%     error = wout - Q;
%     error_judge = sign(error);
%     
%     for i = 2:length(win)-1
%         if error_judge(i-1) == error_judge(i+1)
%             temp = error(i-1)+error(i+1);
%             % temp = (error(i-1)+error(i+1))/2^2;
%         else
%             if abs(error(i-1)) >= abs(error(i+1))
%                 temp = error(i+1) / 2;
%             else
%                 temp = error(i-1) / 2;
%             end
%         end
%         wout(i) = wout(i) + alpha*temp;
%     end 
% end
% ```
% - 如果`stage`大于或等于2，则进行第二轮处理。
% - 再次对经过第一轮调整后的信号`wout`进行判决，计算新的误差和误差符号。
% - 再次遍历信号元素：
%     - 如果当前元素的前一个和后一个误差符号相同，则将前一个和后一个误差相加作为调整值`temp`。
%     - 如果当前元素的前一个和后一个误差符号不同，则比较前一个和后一个误差的绝对值，
%       取绝对值较小的误差的一半作为调整值`temp`。
% - 根据调整值`temp`和参数`alpha`对当前信号元素`wout(i)`进行第二次调整。
% 
% ### 第三轮误差反馈调整（`stage == 3`）
% - 如果`stage`等于3，则进行第三轮处理。
% - 与第二轮处理类似，再次对信号进行判决、计算误差和误差符号，然后遍历信号元素进行调整。
% 
% ### 函数返回
% 函数返回经过处理后的输出信号向量`wout`。
% 
% 总体而言，该函数通过多轮基于误差反馈的调整，根据不同的`stage`设置对输入信号进行逐步优化，以改善信号的判决结果，提高信号质量。但需要注意的是，代码中依赖的`decision`函数需要根据具体的星座图和信号处理需求进行自定义实现。  