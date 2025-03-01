% Test PNC and PEC

if 1
% Encoding 经过编码
% Linear encoding
alpha = 0.5;

% 对信号进行编码，可有可无
% sig信号为初始的信号，后续要经过脉冲成型
lcoeff = [1, alpha];
% Flitering
prs_sig = filter(lcoeff, 1, sig);

% ffe 过后的数据，经过PNC或者PEC

pnc_out = my_pnc(ffe_out, alpha, constellation);

my_pnc_out = my_pnc_enhanced(ffe_out, alpha, constellation);
my_pnc_out = my_pnc_enhanced(my_pnc_out, alpha, constellation);
else
% PRS_PNC

% Linear encoding
D = 0.8;
lcoeff = PRS_poly(D, 2, 1);
lcoeff = lcoeff / sum(lcoeff);

% Flitering
prs_sig = filter(lcoeff, 1, sig);

pnc_out = my_pnc(ffe_out, alpha, constellation);
% out_ffe = mlseeq(ffe_out, lcoeff, constellation, 1000, 'rst');
end