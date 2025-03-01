% Trellis for M = 2, L=2
function [Prev_State,Prev_State_trans,Prev_Ip,Prev_Ip_trans,Outputs_prev]= Get_Trellis()
Prev_State = [1,2;3,4;1,2;3,4]; % PREVIOUS STATES
Prev_State_trans =  [1,3,1,3;2,4,2,4]; % TRANSPOSE OF P_State
Prev_Ip = [1,1;1,1;2,2;2,2]; % Previous inputs
Prev_Ip_trans = [1,1,2,2;1,1,2,2]; % tranpose of PREVIOUS INPUTS
Outputs_prev = [1,2;3,4;5,6;7,8];% GAMMA INDICES FOR THE ALPHA RECURSION (THINK OF BCJR ALGORITHM)


% Prev_State = [4,3;4,3;2,1;2,1]; % PREVIOUS STATES
% Prev_State_trans =  Prev_State.'; % TRANSPOSE OF P_State
% Prev_Ip = [2,2;1,1;2,2;1,1]; % Previous inputs
% Prev_Ip_trans = Prev_Ip.'; % tranpose of PREVIOUS INPUTS
% Outputs_prev = [1,5;2,6;3,7;4,8];
% Outputs_prev =[8,4;6,2;7,3;5,1];
% Outputs_prev =[4,8;2,6;3,7;1,5];
end