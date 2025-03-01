function w = getcutoff(ch, ampdB)
amp = 10^(ampdB/10);
[H, f] = freqz(ch, 1, 2^12);
amx = max(abs(H));
[~, idx] = min(abs(abs(H)-amx/sqrt(amp)));
w = f(idx) / pi;
end