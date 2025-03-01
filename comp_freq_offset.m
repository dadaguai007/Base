function comp_signal = comp_freq_offset(sig, freq_offset, Fs)
    % Compensate for frequency offset in signal


    % Fix number of stuff
    % sig should be the 2×N
    npols = size(sig, 1);

    % Output Vector
    comp_signal = zeros(size(sig));
    Nfft=size(sig,2);

    t = (0:(Nfft-1)) * (1/Fs);

    for l = 1:npols
        lin_phase = 2 * pi * t * freq_offset(l);
        comp_signal(l, :) = sig(l, :) .* exp(-1i * lin_phase);
    end
    
end
