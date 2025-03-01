function bits = demodulateGray(symb, M, constType)
% Demodulate symbol sequences to bit sequences (w/ Gray mapping).
%
% Hard demodulation is based on minimum Euclidean distance.
%
% Parameters
% ----------
% symb : array of complex constellation symbols
%     sequence of constellation symbols to be demodulated.
% M : int
%     order of the modulation format.
% constType : string
%     'qam', 'psk', 'pam' or 'ook'.
%
% Returns
% -------
% array of ints
%     sequence of demodulated bits.

if M ~= 2 && constType == "ook"
    warning("OOK has only 2 symbols, but M != 2. Changing M to 2.")
    M = 2;
end
const = GrayMapping(M, constType); % use the custom function defined before

% get bit to symbol mapping
indMap = minEuclid(const, const); % use the custom function defined before
bitMap = de2bi(indMap, log2(M));
b = log2(M);
bitMap = reshape(bitMap, [M, b]);

% demodulate received symbol sequence
indrx = minEuclid(symb, const); % use the custom function defined before
bits = demap(indrx, bitMap); % use the custom function defined before
end