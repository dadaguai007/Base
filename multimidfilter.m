function y=multimidfilter(x,m,order)
% 中值滤波 操作
if nargin < 3
    order = 5; %5th order median filtering is used by default
end

a=x;
for k=1 : m
    b=medfilt1(a, order); 
    a=b;
end
y=b;
%用于执行多次中值滤波（median filtering）操作。函数的输入参数和输出参数如下：

% **输入参数:**
% - `x`：输入信号。
% - `m`：中值滤波的次数。
% - `order`：中值滤波的阶数（默认为5）。
% 
% **输出参数:**
% - `y`：经过多次中值滤波后的信号。
% 
% 函数的主要过程如下：

% 2. 使用 `medfilt1` 函数对输入信号 `x` 进行一次中值滤波，并将结果存储在 `b` 中。
% 3. 将中值滤波的结果 `b` 赋值给变量 `a`。
% 4. 重复步骤2和3，共执行 `m` 次中值滤波。
% 5. 将最终的中值滤波结果赋值给输出变量 `y`。
% 
% 这个函数的目的是通过多次中值滤波来平滑输入信号，以减少噪声或波动。中值滤波是一种非线性滤波方法，通过取窗口内值的中位数来平滑信号。