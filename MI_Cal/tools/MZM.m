function out = MZM(drive, laser, Vpi, Vdc, Vbias, extinctionRatio, insertionLoss)
% This function is a brief Mach-Zehnder modulator model.
%   Format:out = MZM(drive, laser, Vpi, Vdc, extinctionRatio, insertionLoss)

    eRlin = 10^(extinctionRatio/10);
    e = sqrt(1/eRlin);

    a = sqrt(0.5+e);

    alpha = 10^(insertionLoss/10);

    % Spectralratio = [0.5 0.5];
    Spectralratio = [1 1];
    Vrf_upper = sqrt(Spectralratio(1)) * drive;
    Vrf_lower = sqrt(Spectralratio(2)) * drive;

    laser_upper = laser * a;
    laser_lower = laser * sqrt(1-a^2);

    Vpi_rf = Vpi;
    Vpi_dc = Vdc;
    
    if isscalar(Vbias)
        Vdc_upper = Vbias;
        Vdc_lower = Vbias;
    elseif numel(Vbias) == 2
        Vdc_upper = Vbias(1);
        Vdc_lower = Vbias(2);
    elseif numel(Vbias) > 2
        slog('The length of Vpi array must be 1 or 2', 'ERR');
    end

    phi_upper = pi*Vrf_upper/Vpi_rf + pi*Vdc_upper/Vpi_dc;
    phi_lower = pi*Vrf_lower/Vpi_rf + pi*Vdc_lower/Vpi_dc;

    out = 1/sqrt(2)/sqrt(alpha) * (laser_upper.*exp(1j*phi_upper) + laser_lower.*exp(-1j*phi_lower));

end

