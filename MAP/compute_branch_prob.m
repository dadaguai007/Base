function branch_prob=compute_branch_prob(code_bit_0, code_bit_1, rx_symbol_0, rx_symbol_1,noise_variance)

    code_symbol_0 = 2*code_bit_0 - 1;
    code_symbol_1 = 2*code_bit_1 - 1;

    x = rx_symbol_0 - code_symbol_0;
    y = rx_symbol_1 - code_symbol_1;

%     # Normalized branch transition probability
    branch_prob = exp(-(x*x + y*y)/(2*noise_variance));

end