function bits_de = conv_soft_e4_ampm(sym_in)

% Input : 1*5000
% Output : 1*100000

num_of_states = 8;                                                  
depth_of_trellis = length(sym_in);                                       

bits_de = zeros(1,depth_of_trellis*2);
survivor = NaN(num_of_states,depth_of_trellis + 1);   
branch = NaN(num_of_states,num_of_states);    
path = NaN(num_of_states,depth_of_trellis + 1);            
final_path_state = zeros(1,depth_of_trellis + 1);

% get from trellis 

next_state = [0	 2	1	3;
4	6	5	7;
1	3	0	2;
5	7	4	6;
2	0	3	1;
6	4	7	5;
3	1	2	0;
7	5	6	4];
output = [0	 1	2	3;
4	5	6	7;
0	1	2	3;
4	5	6	7;
0	1	2	3;
4	5	6	7;
0	1	2	3;
4	5	6	7];


%start from state 00
path(1,1) = 0;
path(2:num_of_states,1) = Inf;        

% for each time 
for t = 1:depth_of_trellis

    for state = 0:num_of_states-1
        % for each input bit 0,1,2,3
        for bit = 0:3
            % get next state and possible bit  
            nextstate = next_state(state+1,bit+1);
            possible_bit = output(state+1,bit+1);
            distance = euclidean_dist_ampm(possible_bit,sym_in(t));
            % save distance from state to next state
            branch(state+1,nextstate+1) = distance;
        end
    end

    for state = 0: num_of_states - 1
        % current path + next possible branch
        term = path(:,t) + branch(:, state+1);
        % minimum accumulated path 
        [val, ind] = min(term);
        % save the minimum path to path matrix
        path(state+1, t+1) = val;
        if ~isnan(val)
            % save the survived state 
            survivor(state+1, t+1) = ind - 1;
        end
    end

    branch(:)= NaN;
end


% trace back, the last state = 00
[~, ind] = min(path(:,end));
state = ind - 1;

 for i = depth_of_trellis+1 : -1 : 1
     final_path_state(i) = state;
     state = survivor(state+1, i);
 end

 for i = 1:depth_of_trellis
     possible_state = next_state(final_path_state(i) + 1, :);
     predict_state = final_path_state(i+1);
     % corresponds to next_state column, 
     % 
     bits = find( possible_state == predict_state );
     if bits == 4
         bits_de(2*i-1) = 1;
         bits_de(2*i) = 1;
     elseif bits == 3
         bits_de(2*i-1) = 1;
         bits_de(2*i) = 0;
     elseif bits == 2
         bits_de(2*i-1) = 0;
         bits_de(2*i) = 1;
     else
         bits_de(2*i-1) = 0;
         bits_de(2*i) = 0;
     end
 end

end
