function bits_de = conv_hard(bits_in,G)

% G = [1,0,1,1,1;1,0,1,1,0] for 2rd case

[n,N] = size(G);
num_of_states = 2^(N-1);                                                  
depth_of_trellis = length(bits_in)/n;                                       

bits_de = zeros(1,depth_of_trellis);
survivor = NaN(num_of_states,depth_of_trellis + 1);   
branch = NaN(num_of_states,num_of_states);    
path = NaN(num_of_states,depth_of_trellis + 1);            
final_path_state = zeros(1,depth_of_trellis + 1);


% change the received bit from binary to integer,00-0,01-1,10-2,11-3
y_hat = buffer(bits_in, 2)' ;
y_hat = reshape(bi2de(y_hat ,'left-msb'), 1, []) ;

% input = reshape(bits_in,n,depth_of_trellis);
% input = zeros(num_of_states);

% get from trellis
next_state = [0	8;
0	8;
1	9;
1	9;
2	10;
2	10;
3	11;
3	11;
4	12;
4	12;
5	13;
5	13;
6	14;
6	14;
7	15;
7	15];
output = [0	2;
3	1;
0	2;
3	1;
3	1;
0	2;
3	1;
0	2;
3	1;
0	2;
3	1;
0	2;
0	2;
3	1;
0	2;
3	1];

%start from state 00
path(1,1) = 0;
path(2:num_of_states,1) = Inf;        

% for each time 
for t = 1:depth_of_trellis

    % for each state 00(0),01(1),10(2),11(3)
    for state = 0:num_of_states-1
        % for each input bit 0,1
        for bit = 0:1
            % get next state and possible bit  
            nextstate = next_state(state+1,bit+1);
            possible_bit = output(state+1,bit+1);
            % calculate hamming distance  
            distance = hamming_dist(possible_bit,y_hat(t));
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


% trace back 
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

  
