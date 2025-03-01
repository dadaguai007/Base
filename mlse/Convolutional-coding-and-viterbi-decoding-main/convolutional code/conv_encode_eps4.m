function [bits_en] = conv_encode_eps4(bits_in)

% Convolutional encoder, R=2/3
% Input : 1*100000
% Output: 1*150000

% trellis

u_1 = zeros(1, length(bits_in)/2);
u_2 = zeros(1, length(bits_in)/2);

c_1 = zeros(1, length(bits_in)/2);
c_2 = zeros(1, length(bits_in)/2);
c_3 = zeros(1, length(bits_in)/2);

D_0 = zeros(1, length(bits_in)/2);
D_1 = zeros(1, length(bits_in)/2);
D_2 = zeros(1, length(bits_in)/2);

for i = 1:length(bits_in)/2
    
    u_1(i) = bits_in(2*i-1);
    u_2(i) = bits_in(2*i);


    c_1(i) = D_2(i);
    c_2(i) = u_1(i);
    c_3(i) = u_2(i);


    D_0(i+1) = D_2(i); 
    D_1(i+1) = bitxor(D_0(i),u_2(i));
    D_2(i+1) = bitxor(D_1(i),u_1(i));


end

bits_en = [c_1;c_2;c_3];
bits_en = reshape(bits_en,1,[]);


end