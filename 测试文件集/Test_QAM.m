clc;clear;


M=16;
N=log2(M);
h_bit=N-1;
v_bit=floor(N/2)-1;

z = (0:M-1)';
y = qammod(z,M,'gray');

rng(1)
b=randi([0,1],2000,1);
M=16;
iq1 = qammod(b,M,'InputType','bit');
h = qamdemod(iq1,M,'OutputType','bit');




function iq = qam_mapping(b, M)
b=b.';
bits = log2(M);
    if mod(length(b), bits) ~= 0
        warning('b.size not divisible by %d, zeros appended', bits);
        b = [b, zeros(1, bits - mod(length(b), bits))];
    end

    % get chosen constellation map
    sym2iq_map = qammod((0:M-1)',M,'gray');

    % generate iq points based on chosen constellation map
    idx = zeros(1, length(b)/bits);
    for i = 1:bits
        %二进制左移
        idx = idx + bitshift(b(i:bits:end), i-1);
    end

    iq = sym2iq_map(idx + 1);
end



% h=qamdemod(iq,M,'OutputType','bit');




function b_out = qam_demapping(iq, order)
% order=16;
% bits = log2(order);
    sym2iq_map = qammod((0:order-1)',order,'gray');

    sym2iq_map_mod = (sym2iq_map - min(real(sym2iq_map)) * (1 + 1i))/2;
    demap_map = zeros(1, order);
    demap_map(bitshift((real(sym2iq_map_mod)), bits/2) + (imag(sym2iq_map_mod))+1) = 1:order;

    b_out = zeros(1, bits * length(iq));

    % modify and saturate
    iq_Mod = (iq - min(real(sym2iq_map)) * (1 + 1i))/2;
    iq_modreal = max(real(iq_Mod),0);
    iq_modimag = max(imag(iq_Mod),0);
    iq_modreal = min(iq_modreal , 2.^floor(bits/2)-1);
    iq_modimag = min(iq_modimag , 2.^floor(bits/2)-1);

    iq_idx = bitshift((round(iq_modreal)), bits/2) + (round(iq_modimag));
    sym = demap_map(iq_idx + 1);

    for j = 1:bits
        b_out(j:bits:end) = bitget(sym, j);
    end
end
