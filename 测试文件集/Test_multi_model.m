% WDM 的多个模式模拟
for indMode = 1:param.Nmodes
    %模式数
    fprintf('  mode #%d\t power: %.2f dBm\n', indMode, 10 * log10((Pch(indCh) / param.Nmodes) / 1e-3));

    % Generate random bits
    bitsTx = randi([0, 1], 1, param.Nbits);

    % Map bits to constellation symbols
    %             symbTx = modulateGray(bitsTx, param.M, param.constType);
    %You need to manually change the parameters yourself
    symbTx=qammod(bitsTx,param.M,'InputType','bit','UnitAveragePower',1);

    % Normalize symbols energy to 1
    %             symbTx = symbTx / sqrt(Es);

    symbTxWDM(:, indMode, indCh) = symbTx;

    % Upsampling
    symbolsUp = upsample(symbTx, param.SpS);

    % Pulse shaping
    sigTx = firFilter(pulse, symbolsUp);
    % Optical modulation
    if indMode == 1
        % Generate LO field with phase noise
        phi_pn_lo = phaseNoise(param.lw, length(sigTx), 1 / Fs);
        sigLO = exp(1i * phi_pn_lo);
    end
    sigTxCh = iqm(sigLO, 0.5 * sigTx);
    sigTxCh = sqrt(Pch(indCh) / param.Nmodes) * (sigTxCh/max(abs(sigTxCh).^2));
    %mode and channel sum
    sigTxWDM(:, indMode) = sigTxWDM(:, indMode) + sigTxCh.*exp(1i*2*pi*(freqGrid(indCh)/Fs)* t);

    Pmode = Pmode + signal_power(sigTxCh);
end