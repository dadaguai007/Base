%E3_2_QuantArith.m
x=[7/8 zeros(1,15)];
y=zeros(1,length(x));   %存放原始运算结果
B=2;                    %量化位数
Qy=zeros(1,length(x));   %存放量化运算结果
Qy2=zeros(1,length(x));  %存放量化运算结果
Qy4=zeros(1,length(x));  %存放量化运算结果
Qy6=zeros(1,length(x));  %存放量化运算结果

%系统系数
A=0.5;
b=[1];
a=[1,A];

        %差分方程的作用是根据前一个时刻的输出值和当前时刻的输入信号来计算当前时刻的输出值。
        % 这种差分方程可以用来实现信号的平滑、去噪、滤波等操作。


%未经过量化处理的运算
% -A * Qy(i-1) 表示前一个时刻的输出值经过系数 -A 的缩放，它在滞后滤波器中引入了一个反馈项，有助于滤波器的延迟和稳定性。
for i=1:length(x);
    if i==1
        y(i)=x(i);
    else
        y(i)=-A*y(i-1)+x(i);
    end
end

%经过量化处理的运算
for i=1:length(x);
    if i==1
        Qy(i)=x(i);
        Qy(i)=round(Qy(i)*(2^(B-1)))/2^(B-1);
    else
        Qy(i)=-A*Qy(i-1)+x(i);
        Qy(i)=round(Qy(i)*(2^(B-1)))/2^(B-1);
    end
end
Qy2=Qy;

% 模拟将连续信号（浮点数）映射到离散值（有限位数的整数）的过程，模拟量化器的基本原理
% Qy(i) 乘以 2^(B-1)，这个操作将信号映射到一个范围内，该范围由量化位数 B 决定。
% 范围的上限和下限取决于 B，通常是 [-2^(B-1), 2^(B-1)-1]。
% round 函数用于四舍五入操作，将浮点数转换为最接近的整数
%/ 2^(B-1)：这一部分是为了将经过乘法和四舍五入后的结果重新映射回原始范围，以获得最终的量化输出值
B=4;
%经过量化处理的运算
for i=1:length(x);
    if i==1
        Qy(i)=x(i);
        Qy(i)=round(Qy(i)*(2^(B-1)))/2^(B-1);
    else
        Qy(i)=-A*Qy(i-1)+x(i);
        Qy(i)=round(Qy(i)*(2^(B-1)))/2^(B-1);
    end
end
Qy4=Qy;

B=6;
%经过量化处理的运算
for i=1:length(x);
    if i==1
        Qy(i)=x(i);
        Qy(i)=round(Qy(i)*(2^(B-1)))/2^(B-1);
    else
        Qy(i)=-A*Qy(i-1)+x(i);
        Qy(i)=round(Qy(i)*(2^(B-1)))/2^(B-1);
    end
end
Qy6=Qy;

xa=0:1:length(x)-1;
plot(xa,y,'-',xa,Qy2,'--',xa,Qy4,'O',xa,Qy6,'+');
legend('原系统运算结果','2bit量化运算结果','4bit量化运算结果','6bit量化运算结果')
xlabel('运算次数');ylabel('滤波结果');

%量化位越小，振荡值越大，输出响应在几次运算后形成固定值的来回振荡过程

