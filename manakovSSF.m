function [Ex, Ey] = manakovSSF(Ei, hz, Lspan, Ltotal, alpha, gamma, D, Fc, Fs)
Ex = Ei(:, 1);
Ey = Ei(:, 2);

c = 299792458; % speed of light (vacuum)
c_kms = c / 1e3;
lambda = c_kms / Fc;
alpha = alpha / (10 * log10(exp(1)));
beta2 = -(D * lambda^2) / (2 * pi * c_kms);


Nfft = length(Ex);
w = 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);

Nspans = floor(Ltotal / Lspan);
Nsteps = floor(Lspan / hz);

Ex = fft(Ex); % Pol. X
Ey = fft(Ey); % Pol. Y

linOperator = exp(-(alpha / 2) * (hz / 2) + 1i * (beta2 / 2) * (w.^2) * (hz / 2));

for spanN = 1:Nspans
    for stepN = 1:Nsteps

        % First linear step (frequency domain)
        Ex = Ex .* linOperator;
        Ey = Ey .* linOperator;

        % Nonlinear step (time domain)
        Ex = ifft(Ex);
        Ey = ifft(Ey);

        Ex = Ex .* exp(1i * gamma * 8/9 * (Ex .* conj(Ex) + Ey .* conj(Ey)) * hz);
        Ey = Ey .* exp(1i * gamma * 8/9 * (Ex .* conj(Ex) + Ey .* conj(Ey)) * hz);

        % Second linear step (frequency domain)
        Ex = fft(Ex);
        Ey = fft(Ey);
        Ex = Ex .* linOperator;
        Ey = Ey .* linOperator;
    end

    Ex = Ex .* exp(alpha * Nsteps * hz);
    Ey = Ey .* exp(alpha * Nsteps * hz);
end

Ex = ifft(Ex);
Ey = ifft(Ey);
end
