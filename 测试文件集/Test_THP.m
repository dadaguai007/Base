
% 信道矩阵B，需要得到最优的信道估计效果
% 所以THP编码需要一个信道矩阵
B = [1, 2, 3; 4, 5, 6; 7, 8, 9];
% B=[2.5,1,-0.5];
% 进行QR分解
%进行共轭转置
B=B.';
[Q, R] = qr(B);
% 显示结果
disp('Q:');
disp(Q);
disp('R:');
disp(R);
%选取共轭转置
L=R.';
Q=Q.';
A=2;
% THP 编码 Tomlinson-Harashiam precoding
% N为信号长度，L为进行LQ分解得到的上矩阵
% 第一个编码信号保持不变，后续的码元进行DPC编码，最后进行THP的取模操作
% N 应该为信道的长度，是不是要将信号进行分组进行传输？

N = size(B,1);
coding=zeros(N,1);
send=zeros(N,1);
sigtx = randi([0, 1], 1, 1000);
for i=1:floor(length(sigtx)/N)
    symbol=sigtx((i-1)*N+1:i*N);
    for m = 1:N
        coding(m) = symbol(m);
        if m==1
            coding(m)=coding(m);
        else
            for n = 1:m-1
                coding(m) = coding(m) - coding(n) * L(m, n) / L(m, m);
            end
        end

        coding(m) = THP_modulo(coding(m), A); % Tomlinson-Harashiam precoding
    end

    Q_H=Q.';
    % 完成THP编码过程
    for m = 1:N
        for n = 1:N
            send(m) = send(m) + Q_H(m, n) * coding(n);
        end
    end
    sigTx((i-1)*N+1:i*N)=send;
end

