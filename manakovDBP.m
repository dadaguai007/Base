function [Ech, param] = manakovDBP(Ei, param)
% Default parameters
Ltotal = 400; %km
Lspan = 80;
hz= 0.5;
alpha=0.2;
D = 16;
gamma = 1.3;
Fc = 193.1e12;
amp='edfa';
maxIter = 10;
tol=1e-5;
nlprMethod='true';
maxNlinPhaseRot=2e-2;

% check input parameters

if isfield(param, 'Fs')
    Fs = param.Fs;
end
if isfield(param, 'Ltotal')
    Ltotal = param.Ltotal;
end
if isfield(param, 'Lspan')
    Lspan = param.Lspan;
end
if isfield(param, 'hz')
    hz = param.hz;
end
if isfield(param, 'alpha')
    alpha = param.alpha;
end
if isfield(param, 'D')
    D = param.D;
end
if isfield(param, 'gamma')
    gamma = param.gamma;
end
if isfield(param, 'Fc')
    Fc = param.Fc;
end

if isfield(param, 'amp')
    amp = param.amp;
end
if isfield(param, 'maxIter')
    %积分中的最大迭代次数
    maxIter = param.maxIter;
end

if isfield(param, 'tol')
    %积分中的收敛公差
    tol = param.tol;
end
if isfield(param, 'nlprMethod')
    % 自适应步长基于非线性相位旋转的开关
    nlprMethod = param.nlprMethod;
end
if isfield(param, 'maxNlinPhaseRot')
    %最大非线性相位旋转容差的阈值
    maxNlinPhaseRot = param.maxNlinPhaseRot;
end
%输出的光纤跨度的索引列表
saveSpanN=Ltotal/Lspan;


% channel parameters
c = 299792458; % speed of light (vacuum) in m/s
c_kms = c / 1e3;
lambda = c_kms / Fc;
alpha = alpha / (10 * log(exp(1)));
beta2 = -(D * lambda^2) / (2 * pi * c_kms);

% generate frequency axis
Nfft = length(Ei);
omega = 2 * pi * Fs * (-0.5:1/Nfft:0.5-1/Nfft);
Nspans = floor(Ltotal / Lspan);


Ech_x = Ei(:, 1:2:end).';
Ech_y = Ei(:, 2:2:end).';

% Define static part of the linear operator
% important
argLimOp = (alpha / 2) - 1i * (beta2 / 2) * (omega.^2);

if size(Ech_x, 1) > 1
    argLimOp = repmat(argLimOp, [size(Ech_x, 1), 1]);
else
    argLimOp = reshape(argLimOp, [1, numel(argLimOp)]);
end

if ~isempty(saveSpanN)
    Ech_spans = zeros([size(Ei, 1), size(Ei, 2) * saveSpanN]);
    indRecSpan = 1;
end

for spanN = 1:Nspans
    % Reverse amplification step
    if ismember(amp, {'edfa', 'ideal'})
        Ech_x = Ech_x .* exp(-alpha / 2 * Lspan);
        Ech_y = Ech_y .* exp(-alpha / 2 * Lspan);
    elseif strcmp(amp, 'None')
        Ech_x = Ech_x .* exp(0);
        Ech_y = Ech_y .* exp(0);
    end

    Ex_conv = Ech_x;
    Ey_conv = Ech_y;
    z_current = 0;

    % Reverse fiber propagation steps
    while z_current < Lspan
        Pch = Ech_x .* conj(Ech_x) + Ech_y .* conj(Ech_y);

        phiRot = nlinPhaseRot(Ex_conv, Ey_conv, Pch, gamma);

        if strcmp(nlprMethod, 'true')
            hz_ = maxNlinPhaseRot / max(phiRot);
            if Lspan - z_current < hz_
                hz_ = Lspan - z_current;
            end
        elseif Lspan - z_current < hz
            hz_ = Lspan - z_current;
        else
            hz_ = hz;
        end

        % Define the linear operator
        linOperator = exp(argLimOp * (hz_ / 2));

        % First linear step (frequency domain)
        Ex_hd = ifft(fft(Ech_x) .* linOperator);
        Ey_hd = ifft(fft(Ech_y) .* linOperator);

        % Nonlinear step (time domain)
        for nIter = 1:maxIter
            rotOperator = exp(-1i * phiRot * hz_);

            Ech_x_fd = Ex_hd .* rotOperator;
            Ech_y_fd = Ey_hd .* rotOperator;

            % Second linear step (frequency domain)
            Ech_x_fd = ifft(fft(Ech_x_fd) .* linOperator);
            Ech_y_fd = ifft(fft(Ech_y_fd) .* linOperator);

            % Check convergence of trapezoidal integration in phiRot
            lim = convergenceCondition(Ech_x_fd, Ech_y_fd, Ex_conv, Ey_conv);

            Ex_conv = Ech_x_fd;
            Ey_conv = Ech_y_fd;

            if lim < tol
                break;
            elseif nIter == maxIter
                warning(['Warning: target SSFM error tolerance was not achieved ' ...
                    'in ' num2str(maxIter) ' iterations']);
            end

            phiRot = nlinPhaseRot(Ex_conv, Ey_conv, Pch, gamma);
        end

        Ech_x = Ech_x_fd;
        Ech_y = Ech_y_fd;

        z_current = z_current + hz_;  % Update propagated distance
    end

    if ismember(spanN, saveSpanN)
        Ech_spans(:, (2*indRecSpan-1):(2*indRecSpan)) = [Ech_x.'; Ech_y.'];
        indRecSpan = indRecSpan + 1;
    end
end

if ~isempty(saveSpanN)
    Ech = Ech_spans;
else
    Ech_x = Ech_x.';
    Ech_y = Ech_y.';
    Ech = Ei;
    Ech(:, 1:2:end) = Ech_x;
    Ech(:, 2:2:end) = Ech_y;
end
end
