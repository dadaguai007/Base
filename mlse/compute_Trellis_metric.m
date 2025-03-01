% Computer Trellis vector 计算网格转移矩阵

function Trellis=compute_Trellis_metric(M)

% 对于不同的调制格式，转移矩阵不同
switch M

    case 2
        % 约束长度定为3，索引数为2^3,状态数为2^2
        Prev_State=[1,2;3,4;1,2;3,4];% 上一时刻状态
        Prev_State_trans=Prev_State.';
        Prev_Ip=[1,1;1,1;2,2;2,2]; % 上一时刻输入
        Prev_Ip_trans=Prev_Ip.';
        %Outputs_prev = [1,5;2,6;3,7;4,8];% 分支索引(索引对不上)
        Outputs_prev =[8,4;6,2;7,3;5,1];% 分支索引
        index_temp = [0;1*2;2*2;3*2]; %for linear indexing.
    case 4
        % 约束长度定为2，索引数为4^2, 状态数为4^1
        Prev_State = [1,2,3,4;1,2,3,4;1,2,3,4;1,2,3,4]; % 上一时刻状态
        Prev_State_trans = Prev_State.'; % 转置
%         Prev_Ip = [1,2,3,4;1,2,3,4;1,2,3,4;1,2,3,4]; 
        Prev_Ip_trans=[1,1,1,1;2,2,2,2;3,3,3,3;4,4,4,4]; % 上一时刻输入
        Outputs_prev = [1,5,9,13;2,6,10,14;3,7,11,15;4,8,12,16];% 分支索引
        index_temp = [0;1*4;2*4;3*4]; %for linear indexing
end

Trellis=struct();
Trellis.Prev_State=Prev_State;
Trellis.Outputs_prev=Outputs_prev;
Trellis.Prev_State_trans=Prev_State_trans;
Trellis.Prev_Ip_trans=Prev_Ip_trans;
Trellis.index_temp=index_temp;