function Eo = edfa_lin(signal, gain, nf, fc, fs)
% simple linear EDFA model.
h = 6.62607004e-34;
nf_lin = 10^(nf/10);
gain_lin = 10^(gain/10);
nsp = (gain_lin * nf_lin - 1) / (2 * (gain_lin - 1));
s_ase = (gain_lin - 1) * nsp * h * fc;

p_noise = s_ase * fs;
mean_noise = 0;
noise_real = randn(size(signal));
noise_imag = randn(size(signal));
noise = (mean_noise + sqrt(p_noise) * (noise_real + 1i * noise_imag));

Eo = signal * gain_lin + noise;
end
