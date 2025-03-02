% class EDFA, filename 'EDFA.m'
%   This class is a brief implementation of an EDFA model.
% 
% Function:
%   It can amplify the input signal based on the gain parameter, and add ASE noise to the amplified signal according to the noise figure.
%
% Examples:
%   @code
%   sig = rand(1e3, 2);
%   param.sig = struct('Fs', 20e9, 'Fc', const.c/1550e-9, 'Rs', 20e9, 'PCol', '[pwr(20, {-16, 'dBm}), pwr(25, {-18, 'dBm'})]);
%   sigfield = signal_interface(sig, param.sig)
%   param.edfa.gain = 16;
%   param.edfa.NF = 5;
%   param.edfa.gaintype = 'dB';
%   edfa = EDFA(param.edfa);
%   edfafield = edfa.generate(sigfield);
%   @endcode
%
% from @OCG-LAB OPC Framework
%
% 简体中文版注释：
% 类 EDFA，文件名为 'EDFA.m'
%   该类是一个简要的EDFA模型实现。
% 
% 功能：
%   它可以根据增益参数对输入的信号进行放大，并且根据噪声系数为放大后的信号添加ASE噪声。
%
% 示例：
%   @code
%   sig = rand(1e3, 2);
%   param.sig = struct('Fs', 20e9, 'Fc', const.c/1550e-9, 'Rs', 20e9, 'PCol', '[pwr(20, {-16, 'dBm'}), pwr(25, {-18, 'dBm'})]');
%   sigfield = signal_interface(sig, param.sig)
%   param.edfa.gain = 16;
%   param.edfa.NF = 5;
%   param.edfa.gaintype = 'dB';
%   edfa = EDFA(param.edfa);
%   edfafield = edfa.generate(sigfield);
%   @endcode
%
% 来自 @OCG-LAB OPC 框架

classdef EDFA < unit
    
    properties
        nInputs = 1;
        nOutputs = 1;
        gain;       % Gain
        gaintype;   % gaintype, default:'dB'
        NF;         % noise figure, default:3
    end

    methods
        function obj = EDFA(param)
            if isfield(param, 'gain') && isscalar(param.gain), obj.gain = param.gain; else, slog('Please input "gain"(dB).', 'ERR'); end
            if isfield(param, 'NF') && isscalar(param.NF), obj.NF = param.NF; else, obj.NF = 3; end
            if isfield(param, 'gaintype') && ischar(param.gaintype), obj.gaintype = param.gaintype; else, obj.gaintype = 'dB'; end
        end

        function out = generate(obj, in)
            NFlin = 10^(obj.NF/10);
            switch obj.gaintype
                case 'dB'
                    Glin = 10^(obj.gain/10);
                case 'lin'
                    Glin = obj.gain;
                otherwise
                    slog('Please input the data type of gain "gaintype", it should be "lin" or "dB".', 'ERR');
            end
            nsp = (Glin*NFlin-1) / (2*(Glin-1));
            np = (Glin-1) .* nsp * const.h * in.Fc;
            P = max(0, np*in.Fs);
            noise = wgn(in.L, in.N, P, 'linear', 'complex');
            PCol_noise = repmat(pwr(-inf, {P, 'W'}), 1, in.N);
            PCol_out = in.PCol * Glin + PCol_noise;

            Sout = sqrt(Glin) * get(in) + noise;
            out = signal_interface(Sout, struct('Fs',in.Fs,'Rs',in.Rs, 'PCol', PCol_out, 'Fc', in.Fc));
        end

    end

end