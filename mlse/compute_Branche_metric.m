% Computer Branch metric 计算分支度量

function [Branches,Branche_vector]=compute_Branche_metric(input,M,const,pulse_length,lecoff)

% 分支矩阵
Branche_vector =get_state(M, const, pulse_length);
% 度量存储矩阵
Branches                = zeros(M.^pulse_length,length(input)-pulse_length+1);

% 计算度量值
for i=1:length(Branche_vector)
    % 选出所有的组合向量
    vector=Branche_vector(:,i);
    % 计算对应的度量值
    Branches(i,:)=input(pulse_length:end)-lecoff*vector;
end

end



% 首先要明确所有的状态数和分支数

% % 所有状态数
% States_Number=M.^(pulse_length-1)*k;
%
% % 所有分支数
% Branche_Number=M.^(pulse_length)*k;


%  每个调制电平状态下状态数
% state_vector_number  = M.^(pulse_length-1);
%  每个调制电平状态下分支数  分支数应为状态数的M倍
% Branche_vector_number = M.^(pulse_length);

% 构建状态矩阵
% state_vector=get_state(M, const, pulse_length-1)
% states                  = zeros(States_Number,(pulse_length-1)+1);
% states(1:end,1)         = repelem(phase_states,state_vector_number);  % 第一列重复所有 调制格式 state_vector_number次
% states(1:end,2:end)     = repmat(state_vector,length(phase_states),1);

% 构建分支矩阵
% Branche_vector =get_state(M, const, pulse_length)
% Branches                = zeros(Branche_Number,pulse_length+1);
% Branches(1:end,1)       = repelem(phase_states,Branche_vector_number);   % 第一列重复所有 调制格式 Branche_vector_number次
% Branches(1:end,2:end)   = repmat(Branche_vector,length(phase_states),1);


% 构建状态转移网格（Trellis），即确定每个状态如何转移到下一个状态，以及对应的分支路径
% 保存起始状态索引、终止状态索引和分支索引
% states_transitions       = zeros(length(Branches),3);
% for i =1:length(Branches)
%states_transitions(i,1)=intersect(Branches(i,1:pulse_length),states,'rows');;% 起始索引
% states_transitions(i,2)=intersect([Next Branches(i,3:end)],states,'rows');;% 终止索引
% states_transitions(i,3)=i;% 分支索引
% end
% 排序（按照终止索引进行）
% [~,idx]              = sort(states_transitions(:,2), 'ascend');
% states_sort          = states_transitions(idx,:);                % Comparable Branches are sorted one after another.


