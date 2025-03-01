function bits_de = conv_soft_e1_qpsk(sym_in,G)



% G3 = [1,1,0,0,1;1,1,0,1,1]; for 3nd case
% QPSK : constellation = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 - 1i)]/sqrt(2);
% 00: +1+1, 01: +1-1, 10: -1+1, 11: -1-1


[n,N] = size(G);
num_of_states = 2^(N-1);                                                  
depth_of_trellis = length(sym_in);                                       

bits_de = zeros(1,depth_of_trellis);
survivor = NaN(num_of_states,depth_of_trellis + 1);   
branch = NaN(num_of_states,num_of_states);    
path = NaN(num_of_states,depth_of_trellis + 1);            
final_path_state = zeros(1,depth_of_trellis + 1);

% get from trellis 
% 0: 00, 1:01, 2:10, 3:11 etc
next_state = [0,2;0,2;1,3;1,3];
output = [0,3;3,0;1,2;2,1];

% QPSK
% constellation = [1 + 1i, 1 - 1i, -1 + 1i , -1 - 1i]/sqrt(2);

%start from state 00
path(1,1) = 0;
path(2:num_of_states,1) = Inf;        

% for each time 
for t = 1:depth_of_trellis

    for state = 0:num_of_states-1
        % for each input bit 0,1
        for bit = 0:1
            % get next state and possible bit  
            nextstate = next_state(state+1,bit+1);
            possible_bit = output(state+1,bit+1);
            distance = euclidean_dist(possible_bit,sym_in(t));
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
     % 1st column -> bit 0, 2nd column -> bit 1
     bits = find( possible_state == predict_state ) - 1;
     bits_de(i) = bits;
 end

end
