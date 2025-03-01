% Test PRS  部分响应信号序列
% Rs>>Fs
clc;clear;close all;
% 参数
sps = 10;
Rs = 510e9;
Fs = 256e9;
% 信号生成
M=4;
data_2bit=randi([0,1],log2(M),8000);
% 相当于四个电平
symbols = 2.^(0:log2(M)-1)*data_2bit;
symbols=symbols.';
% Mapeia bits para pulsos eletricos
symbTx = pammod(symbols,M,0,'gray');
% symbTx = pnorm(symbTx);

% Pulso
hsqrt = rcosdesign(0.01,256,sps,'sqrt');

if 1
    if 0
        % Encoding
        % Linear encoding
        lcoeff = PRS_poly(0.5, 4, 1);
        lcoeff = lcoeff / sum(lcoeff);
        % nonlinear encoding
        nlcoeff = exp(1j*gaussdesign(0.5, 2, 1));
        nlcoeff = nlcoeff / sum(nlcoeff);
        nlcoeff = [];
        step = 0.1;
        w = [lcoeff, step*nlcoeff];
        numtaps = [length(lcoeff), length(nlcoeff)];

        % Flitering
        L = length(symbTx) + 1 - 5;
        prs_sig = zeros(L, 1);
        block = zeros(length(w), 1);
        % temp_2 = zeros(2, 1);
        for i = 1:L
            block(1:numtaps(1)) = [symbTx(i); block(1:numtaps(1)-1)];
            % temp_2 = [sig(i); temp_2(1:end-1)];
            % block(sum(numtaps(1))+1:sum(numtaps(1:2))) = gen_vol_block(temp_2, 2);
            prs_sig(i) = w * block;
        end
    else
        D = 0.25;
        taps1 = 5;
        taps2 = 3;
        taps3 = 1;
        % Linear encoding
        lcoeff = PRS_poly(D, taps1-1,1);

        % Nonlinear encoding
        nlcoeff = randn(1, taps2*(taps2+1)/2 + taps3*(taps3+1)*(taps3+2)/6);
        nlcoeff = nlcoeff / sum(nlcoeff);
        w = [lcoeff, nlcoeff];

        % Flitering (PRS filter)
        step = 1e-2;
        prs_sig = zeros(8000, 1);
        for i = 1:8000+1-taps1
            block1 = symbTx(i:i+taps1-1);
            block2 = step * voltrans(block1(round((taps1-taps2)/2)+1 : end - fix((taps1-taps2)/2)), 2).';
            block3 = step.^2 * voltrans(block1(round((taps1-taps3)/2)+1 : end - fix((taps1-taps3)/2)), 3).';
            block = [block1; block2; block3];
            prs_sig(i+taps1-1) = w * block;
        end
    end

else

    D = 1;
    taps1 = 2;
    % Linear encoding
    lcoeff = PRS_poly(D, taps1-1, 1);
    lcoeff = lcoeff / sum(lcoeff);

    % Flitering
    prs_sig = filter(lcoeff, 1, symbTx);

end

% Upsampling
symbolsUp = upsample(prs_sig,sps);
% % pulse shaping
ps_sig=conv(symbolsUp,hsqrt,'same');

% Anti-aliasing filtering
aaf_sig = AAF(ps_sig, Fs, Rs*sps);

% Resampling
sig = resample(aaf_sig, Fs, Rs*sps);

figure;
plot(real(sig))