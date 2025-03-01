function Pb = theoryBER(M, EbN0dB, constType)
    % Theoretical (approx.) bit error probability for PAM/QAM/PSK in AWGN channel.

    % Input arguments:
    % M: Modulation order.
    % EbN0: Signal-to-noise ratio (SNR) per bit in dB.
    % constType: Modulation type: 'pam', 'qam' or 'psk'

    % Output:
    % Pb: Theoretical probability of bit error.
%     EbN0 = SNR/;

%snrdB  = EbN0dB + 10*log10(log2(M))

    EbN0lin = 10^(EbN0dB/10);
    k = log2(M);

    if strcmp(constType, 'qam')
        L = sqrt(M);
        Pb = 2 * (1 - 1/L) / log2(L) * qfunc(sqrt(3 * log2(L) / (L^2 - 1) * (2 * EbN0lin)));
    elseif strcmp(constType, 'psk')
        Ps = 2 * qfunc(sqrt(2 * k * EbN0lin) * sin(pi/M));
        Pb = Ps / k;
    elseif strcmp(constType, 'pam')
        Ps = (2 * (M - 1) / M) * qfunc(sqrt(6 * log2(M) / (M^2 - 1) * EbN0lin));
        Pb = Ps / k;
    end
end
