function [L_ext,decoded_bits]=forward_recursion_decoding(trellis, mode, msg_length, b_state_metrics,f_state_metrics, branch_probs, app, L_int,priors, L_ext, decoded_bits)


number_states = trellis.number_states;
number_inputs = trellis.number_inputs;


next_state_table = trellis.next_state_table;


%     # Forward Recursion
for time_index =2: msg_length+1   % 序号可以考虑进行更改
    % 初始化
    app(:) = 0;
    for current_state =1:number_states
        for current_input =1: number_inputs
            next_state = next_state_table(current_state, current_input);
            branch_prob = branch_probs(current_input, current_state, time_index-1);
            %                 # Compute the forward state metrics
            f_state_metrics(next_state, 1 )=   f_state_metrics(next_state, 1 )+ ...
                (f_state_metrics(current_state, 9) *branch_prob *priors(current_input, time_index-1));

            %             # Compute APP
            app(current_input) =app(current_input)+ ...
                (f_state_metrics(current_state, 1) *branch_prob *b_state_metrics(next_state, time_index));
        end
        % LLR 计算
        lappr = L_int(time_index-1) + log(app(2)/app(1));% LLR_bit1/LLR_bit0 ;
        L_ext(time_index-1) = lappr;
    end
    if strcmp(mode,'decode')
        if lappr > 0
            decoded_bits (time_index-1)= 1;
        else
            decoded_bits(time_index-1) = 0;
        end
    end
    %         # Normalization of the forward state metrics  归一化
    f_state_metrics(:,1) = f_state_metrics(:,1)/sum(f_state_metrics(:,1));

    f_state_metrics(:,0) = f_state_metrics(:,1);
    f_state_metrics(:,1) = 0;

end