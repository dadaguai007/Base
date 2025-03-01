function [filtercoeff, wc] = PRS_poly(D, span, sps)
% 示例：% Linear encoding
% D = 0.8;
% lcoeff = PRS_poly(D, 2, 1);
% lcoeff = lcoeff / sum(lcoeff);
syms x;
taps = span*sps;
y = (1+x)^taps;
coeff = double(flip(coeffs(expand(y))));
clear x y;
filtercoeff = 1;
for i = 1:taps
    filtercoeff = [filtercoeff, D^i];
end
filtercoeff = filtercoeff .* coeff;
% filtercoeff = filtercoeff / sum(filtercoeff);
wc = getcutoff(filtercoeff, 3);
% fprintf('归一化截止频率为%f\n', wc);
end

% 用于生成特定的多项式滤波器系数，并计算其截止频率
% 
% ### 函数定义与输入输出

% - **输入参数**：
%     - `D`：一个参数，可能用于控制滤波器系数的某种特性。
%     - `span`：滤波器的跨度，用于确定多项式的阶数。
%     - `sps`：每个符号的采样点数，与滤波器的抽头数相关。
% - **输出参数**：
%     - `filtercoeff`：生成的滤波器系数向量。
%     - `wc`：计算得到的滤波器截止频率。
% 
% ### 符号运算与多项式系数生成
% ```matlab
% syms x;
% taps = span*sps;
% y = (1+x)^taps;
% coeff = double(flip(coeffs(expand(y))));
% clear x y;
% ```
% - `syms x;`：定义符号变量`x`，用于后续的符号运算。
% - `taps = span*sps;`：计算滤波器的抽头数，由`span`和`sps`相乘得到。
% - `y = (1+x)^taps;`：构建一个多项式`(1 + x) ^ taps`。
% - `coeff = double(flip(coeffs(expand(y))));`：
%     - `expand(y)`：展开多项式`y`。
%     - `coeffs(expand(y))`：获取展开后多项式的系数。
%     - `flip`：将系数向量翻转，因为`coeffs`函数返回的系数顺序可能与预期不符。
%     - `double`：将符号系数转换为数值类型。
% - `clear x y;`：清除不再需要的符号变量`x`和`y`，释放内存。
% 
% ### 生成滤波器系数向量
% ```matlab
% filtercoeff = 1;
% for i = 1:taps
%     filtercoeff = [filtercoeff, D^i];
% end
% filtercoeff = filtercoeff.* coeff;
% ```
% - `filtercoeff = 1;`：初始化滤波器系数向量`filtercoeff`，第一个元素为1。
% - `for i = 1:taps`：循环`i`从1到`taps`。
%     - `filtercoeff = [filtercoeff, D^i];`：将`D`的`i`次幂依次添加到`filtercoeff`向量中，形成一个包含`1, D, D^2,..., D^taps`的向量。
% - `filtercoeff = filtercoeff.* coeff;`：将上述生成的向量与之前得到的多项式系数向量`coeff`对应元素相乘，得到最终的滤波器系数向量。
% 
% ### 计算截止频率
%`wc = getcutoff(filtercoeff, 3);`：调用`getcutoff`函数（该函数需自定义实现），传入滤波器系数向量`filtercoeff`和一个参数`3`，计算得到滤波器的截止频率`wc`。

% 总体而言，这段代码通过符号运算生成多项式系数，
% 结合输入参数`D`构建滤波器系数向量，并通过调用`getcutoff`函数计算其截止频率。 