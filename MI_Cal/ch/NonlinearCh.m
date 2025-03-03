% class NonlinearCh, filename 'NonlinearCh.m'
%   This class is a brief nonlinear optical fiber channel model.
% 
% Function:
%   This class implements the optical channel propagation model.
%   The theoretical model of this class is based on the Non-Linear Schrödinger Equation(NLSE).
%   This class employs a split-step Fourier method(SSFM) for modeling dispersion and nonlinearity.
%
% Attention:
%   1.Pay attention to the units of input parameters.
%   2.Ensure that the input data has two columns. This class can only handle two polarization mode signals, meaning the first two columns of the signal.
%   3.Regarding the addition of polarization rotation:ignore.
%
% Param conventions:
%   1.Dispersion coeff: D = 2*pi*c/lambda^2.
%   2.Dispersion slope: S = (2*pi*c/lambda^2)^2*beta3 + (4*pi*c/lambda^3)*beta2.
%
% Examples:
%   @code
%   sig = rand(1e3, 2) + 1j*rand(1e3, 2);
%   param.sig = struct('Fs', 20e9, 'Fc', const.c/1550e-9, 'Rs', 20e9, 'PCol', '[pwr(20, {-16, 'dBm'}), pwr(25, {-18, 'dBm'})]');
%   fiberch = NonlinearCh();
%   chfield = fiberch.generate(sigfield);
%   @endcode
%
% from @OCG-LAB OPC Framework
%
% 简体中文版注释：
% 类 NonlinearCh, 文件名为 'NonlinearCh.m'
%   该类是一个简要的非线性光纤信道模型。
% 
% 功能:
%   该类实现了光信道传播模型。
%   该类的理论模型基于非线性薛定谔方程(NLSE)。
%   该类采用分步长傅里叶方法(SSFM)来建模色散和非线性。
%
% 注意事项:
%   1. 注意输入参数的单位。
%   2. 确保输入数据有两列。该类只能处理双偏振信号，即信号的前两列。
%   3. 关于偏振旋转的添加方式：忽略。
%
% 参数约定:
%   1. 色散系数: D = 2*pi*c/lambda^2。
%   2. 色散斜率: S = (2*pi*c/lambda^2)^2*beta3 + (4*pi*c/lambda^3)*beta2。
%
% 示例:
%   @code
%   sig = rand(1e3, 2) + 1j*rand(1e3, 2);
%   param.sig = struct('Fs', 20e9, 'Fc', const.c/1550e-9, 'Rs', 20e9, 'PCol', '[pwr(20, {-16, 'dBm'}), pwr(25, {-18, 'dBm'})]');
%   fiberch = NonlinearCh();
%   chfield = fiberch.generate(sigfield);
%   @endcode
%
% 来自 @OCG-LAB OPC 框架

classdef NonlinearCh < unit

    properties
        nInputs = 1;
        nOutputs = 1;
        spans;                          % fiber + amplifier, default:1
        L;                              % length of each span (km), default:80
        stepsize;                       % stepsize of SSFM (km), default:10
        alpha;                          % fiber attenuation (km^-1), default:0.2
        gamma;                          % nonlinear coefficient (W^-1*km^-1), default:1.2
        D;                              % dispersion coefficient (ps/nm/km), default:17
        S;                              % dispersion slope (ps/nm^2/km), default:0
        idealEDFAEnabled                % EDFA without noise
        EDFAEnabled;                    % Whether to enable EDFA, default:1
        EDFAgain;                       % Gain of EDFA, default:16
        EDFAgaintype;                   % Gaintype of EDFA, default:'dB'
        EDFANF;                         % noise figure of EDFA, default:3
        dispersionEnabled;              % add dispersion, default:1
        dispersionCompensationEnabled;  % add dispersion compensation, default:0
    end

    properties(Hidden=true, SetAccess = public)
        beta2cal;
        beta3cal;
    end

    methods
        function obj = NonlinearCh(param)
            if isfield(param, 'spans') && isscalar(param.spans), obj.spans = param.spans; else, obj.spans = 1; end
            if isfield(param, 'L') && isscalar(param.L), obj.L = param.L ; else, obj.L = 80; end
            if isfield(param, 'stepsize') && isscalar(param.stepsize), obj.stepsize = param.stepsize ; else, obj.stepsize = 10; end
            if isfield(param, 'alpha') && isscalar(param.alpha), obj.alpha = param.alpha ; else, obj.alpha = 0.2; end
            if isfield(param, 'gamma') && isscalar(param.gamma), obj.gamma = param.gamma ; else, obj.gamma = 1.2; end
            if isfield(param, 'D') && isscalar(param.D), obj.D = param.D ; else, obj.D = 17; end
            if isfield(param, 'S') && isscalar(param.S), obj.S = param.S ; else, obj.S = 0; end
            if isfield(param, 'idealEDFAEnabled') && isscalar(param.idealEDFAEnabled), obj.idealEDFAEnabled = param.idealEDFAEnabled; else, obj.idealEDFAEnabled = 0; end
            if isfield(param, 'EDFAEnabled') && isscalar(param.EDFAEnabled), obj.EDFAEnabled = param.EDFAEnabled ; else, obj.EDFAEnabled = 1; end
            if isfield(param, 'EDFAgain') && isscalar(param.EDFAgain), obj.EDFAgain = param.EDFAgain ; else, obj.EDFAgain = 16; end
            if isfield(param, 'EDFAgaintype') && isscalar(param.EDFAgaintype), obj.EDFAgaintype = param.EDFAgaintype ; else, obj.EDFAgaintype = 'dB'; end
            if isfield(param, 'EDFANF') && isscalar(param.EDFANF), obj.EDFANF = param.EDFANF ; else, obj.EDFANF = -15.98; end
            if isfield(param, 'dispersionCompensationEnabled') && isscalar(param.dispersionCompensationEnabled), obj.dispersionCompensationEnabled = param.dispersionCompensationEnabled ; else, obj.dispersionCompensationEnabled = 0; end
        end

        function out = generate(obj, in)
            % Calculate key parameters
            ckms       = const.c*1e-3;
            lambda     = ckms / in.Fc;
            beta2      = -obj.D * lambda^2 / (2*pi*ckms);   
            % beta3      = beta2.^2./obj.D.*(obj.S./obj.D + 2/lambda);
            beta3      = lambda^2/(2*pi*ckms)^2*(lambda^2*obj.S+2*lambda*obj.D);
            % beta2 = 0;
            % beta3 = 0;

            % ll = const.c / in.Fc;
            % lambda1 = ll;
            % beta31      = lambda1^2/(2*pi*const.c)^2*(lambda1^2*0+2*lambda1*16e-6);

            % beta3      = 0;
            beta       = [zeros(2,1); beta2; beta3]; 
            nz         = ceil(obj.L/obj.stepsize);        
            dz         = obj.L / nz;                          
            alphalin   = obj.alpha/(10*log10(exp(1)));

            obj.beta2cal = beta2;
            obj.beta3cal = beta3;

            slog(['Nonlinear channel params are initialized:\n' ...
            '\t\tAverage step size: %d km\n' ...
            '\t\tNumber of fiber spans: %d.'], dz, obj.spans);
            
            % Multi-span transmission
            for idx = 1:obj.spans
                sig = in.get;
                slog('Span %d input power: %1.2f dBm.', idx, in.P.Ptot);
                Pin = mean(pwr.meanpwr(sig));
                % SSFM adding dispersion and nonlinearity.
                sig = ssprop(sig, in.Fs, dz, nz, alphalin, beta, obj.gamma);
                Pout = mean(pwr.meanpwr(sig));
                PColout = in.PCol .* (Pout./Pin);
                inparams = in.params;
                inparams.PCol = PColout;
                in = signal_interface(sig, inparams);
                slog('Span %d output power: %1.2f dBm.', idx, in.P.Ptot);
                if obj.EDFAEnabled
                    if obj.idealEDFAEnabled
                        param.gain.gain     = obj.EDFAgain;
                        param.gain.gainType = obj.EDFAgaintype;
                        gain = Gain(param.gain);
                        in = gain.generate(in);
                        slog('Span %d EDFA output power: %1.2f dBm.', idx, in.P.Ptot);
                    else
                        param.edfa.gain     = obj.EDFAgain;
                        param.edfa.gaintype = obj.EDFAgaintype;
                        param.edfa.NF       = obj.EDFANF;
                        edfa = EDFA(param.edfa);
                        in = edfa.generate(in);
                        slog('Span %d EDFA output power: %1.2f dBm.', idx, in.P.Ptot);
                    end
                end
    
                %  Dispersion compensation fiber(DCF).
                if obj.dispersionCompensationEnabled
                    for n = 1:in.N
                        sig(:, n) = ssprop(sig(:, n), in.Ts, dz, nz, 0, -beta, 0);
                    end
                    Power = pwr.meanpwr(sig);
                    slog('Span %d CD output power: %1.2f dBm.', idx, 10*log10((Power/1e-3)));
                    
                    in = signal_interface(sig, in.params);
                end
            end
            out = in;
        end
    end
end