function [wout, snr_dB] = quantizer(win, bits, ENoB, type, snr_dB)
[row, col] = size(win);
win = win(:);
if ~isreal(win)
    step = [(max(real(win))-min(real(win))) / round(2^bits-1), (max(imag(win))-min(imag(win))) / round(2^bits-1)];
    switch lower(type)
        case 'riser'
            wout = step(1)*(round(real(win)/step(1))+0.5) + 1j*step(2)*(round(imag(win)/step(2))+0.5);
        case 'tread'
            wout = step(1)*sign(real(win)).*(round(abs(real(win))/step(1)+0.5)-0.5) + 1j*step(2)*sign(imag(win)).*(round(abs(imag(win))/step(2)+0.5)-0.5);
        otherwise
            error('The quantizer only supports two quantization methods: riser and tread.');
    end
else
    step = (max(win)-min(win)) / round(2^bits-1);
    switch lower(type)
        case 'riser'
            wout = step * (round(win/step)+0.5);
        case 'tread'
            wout = step * sign(win).*(round(abs(win)/step+0.5)-0.5);
        otherwise
            error('The quantizer only supports two quantization methods: riser and tread.');
    end
end
dpwr = bandpower(wout - win);
ipwr = bandpower(win);
pn = ipwr / (3*2^(2*ENoB-1)-1) - dpwr;
pn(pn<0) = 0;
if ~isreal(win)
    noise = wgn(size(win, 1), 1, pn, 'linear', 'complex');
else
    noise = wgn(size(win, 1), 1, pn, 'linear', 'real');
end
if snr_dB == inf
    snr_dB = 10*log10(bandpower(wout)/bandpower(noise));
else
    snr = 10^(snr_dB/10);
    p = bandpower(wout);
    pn0 = p / (snr+1);
    ps = p - pn0;
    pn = pn0 + bandpower(noise);
    snr_dB = 10*log10(ps/pn);
end
wout = wout + noise;
wout = reshape(wout, row, col);
end