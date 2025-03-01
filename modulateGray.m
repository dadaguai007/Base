function symbols = modulateGray(bits, M, constType)
% Modulate bit sequences to constellation symbol sequences (w/ Gray mapping).
%
% Parameters
% ----------
% bits : sequence of data bits.
% M : order of the modulation format.
% constType :  'qam', 'psk', 'pam' or 'ook'.

if M ~= 2 && constType == "ook"
    warning("OOK has only 2 symbols, but M != 2. Changing M to 2.")
    M = 2;
end
bitsSymb = log2(M);
const = GrayMapping(M, constType); % use the custom function defined before

symb = reshape(bits, [bitsSymb, length(bits) / bitsSymb])';
symbInd = bi2de(symb);

symbols = const(symbInd + 1); 
end