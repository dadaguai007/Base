% Viterbi算法
function dec_ip=Viterbi_Decode(num_input,M,States_Number,decoding_delay,branch_metric,Trellis)
% 初始设置 num_input为输入数组的数量

% 网格参数初始化
Prev_State=Trellis.Prev_State;
Outputs_prev=Trellis.Outputs_prev;
Prev_State_trans=Trellis.Prev_State_trans;
Prev_Ip_trans=Trellis.Prev_Ip_trans;
index_temp=Trellis.index_temp;

ip=0; % initialization (原始输入)
dec_ip = zeros(1,num_input-decoding_delay); % decoded symbols (-decoding delay is because,ignoring last transient part)
survivor_node = zeros(States_Number,num_input); % survivor nodes 幸存路径
survivor_ip = zeros(States_Number,num_input); % survivor inputs  幸存输入
path_metric = zeros(States_Number,num_input+1); % path metrics
% 度量存储矩阵，进行判定

for  n=1:num_input
    % 选出度量值较小,  按行比较

    switch M
        % OOK
        case 2
            [path_metric(:,n+1),idx]=min([path_metric(Prev_State(:,1),n)+branch_metric(Outputs_prev(:,1),n) ...
                path_metric(Prev_State(:,2),n)+branch_metric(Outputs_prev(:,2),n)],[],2);
        case 4
            % PAM4
            [path_metric(:,n+1),idx]=min([path_metric(Prev_State(:,1),n)+branch_metric(Outputs_prev(:,1),n) ...
                path_metric(Prev_State(:,2),n)+branch_metric(Outputs_prev(:,2),n) ...
                path_metric(Prev_State(:,3),n)+branch_metric(Outputs_prev(:,3),n) ...
                path_metric(Prev_State(:,4),n)+branch_metric(Outputs_prev(:,4),n)],[],2);
    end
    % 对幸存路径和幸存输入进行赋值
    survivor_node(:,n)=Prev_State_trans(idx+index_temp);
    survivor_ip(:,n)=Prev_Ip_trans(idx+index_temp);
    % 进行回溯
    if (n>decoding_delay)
        % 找到最优路径的索引
        [~,trace_back] = min(path_metric(:,n+1));
        % 在幸存路径和幸存输入中进行索引
        for bk_cnt= 1 : decoding_delay+1
            ip = survivor_ip(trace_back,n+1-bk_cnt);
            trace_back = survivor_node(trace_back,n+1-bk_cnt);
        end
       % 输出
        dec_ip(sym_cnt-decoding_delay)=ip;
    end

end

end

% DetectedBits         = zeros(1,Nbits - decision_delay); 存储最终检测到的比特序列，长度为总比特数减去决策延迟
% DBits_idx            = 1; % 检测比特的索引
% pathmetric           = zeros(States_Number,1); 当前时刻各状态的路径度量
% pathmetric_n         = zeros(States_Number,1); 下一时刻各状态的路径度量
% SurvivorPath         = zeros(States_Number,decision_delay); 用于存储幸存路径信息，记录每个状态的前一个状态
% decision_delay    回溯长度
% for n=1:1:Nbits-pulse_length
% Branche_metrics_cum   = zeros(Branche_Number,1);
% 遍历所有分支，计算所有分支的分支度量
% for i = 1:Branche_Number
% （带入度量计算公式）
%     if(n~=1)
%         z = real(exp(-1i*Branches(i,1))*trapz(received_signal((n-1)*os:n*os).*conj(Branche_metrics(i,:)))*Ts);
%         Branche_metrics_cum(i) = z;
%     else
%         z = real(exp(-1i*Branches(i,1))*trapz(received_signal(1:1*os+1).*conj(Branche_metrics(i,:)))*Ts);
%         Branche_metrics_cum(i) = z;
%     end
% end
% 路径度量更新：
% for i=1:M_ary:总分支数
%     for k=0:M_ary-1
%         分支度量比较(k+1) = Branche_metrics_cum(分支索引) + 前驱状态度量;
%     end
%     [新路径度量, 最优分支索引] = max(分支度量比较);
%     SurvivorPath(目标状态, 当前时间步) = 前驱状态索引;（直接将前一个状态的索引记录到 SurvivorPath 矩阵的对应位置）
% end
% 更新前驱状态度量， 即前驱状态度量等于新路径度量
% 当时间步超过 decision_delay 时，滑动窗口更新 SurvivorPath，丢弃最早的路径记录

% 回溯操作
% if(n>decision_delay)
%-------------------- trace back unit----------------------------
%     [ ~，idxpath] = max(pathmetric);最大状态量
%     currState         = idxpath;
%     for jj = decision_delay:-1:1
%         prevState         = SurvivorPath(currState,jj);% 从 SurvivorPath 中回溯决策延迟步，找到最可能的路径。
%         if(jj>1)
%              currState         = prevState;
%         end
%     end
% 对当前最优路径进行回溯：
% 从当前最优状态 currState 出发，沿着 SurvivorPath 逆向查找前驱状态。
% 回溯深度为 decision_delay 步，确保路径收敛。
%     DetectedBits(DBits_idx) = Branches(states_sort(find(ismember(states_sort(:,1:2),[prevState currState],'rows')),3),2);
% 根据回溯路径中的分支索引（计算出分支索引），提取对应的输入符号
%     DBits_idx = DBits_idx+1;
% end