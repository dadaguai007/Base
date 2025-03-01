function code = GrayCode(n)
    % Gray code generator.
    %
    % Parameters:
    %   n: length of the codeword in bits.
    %
    % Returns:
    %   code: array of binary strings of the gray code.

    code = cell(1, 2^n); % Initialize as zero vector

    for i = 0:(2^n - 1)
        % Generating the decimal
        % values of gray code then using
        % bitxor to convert them to binary form
        val = bitxor(i, bitshift(i, -1));

        % Converting to binary string
        s = dec2bin(val, n);
        code{i + 1} = s;
    end
%     code = cell2mat(code);
end
