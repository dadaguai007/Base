% class Gain, filename 'Gain.m'
%   This class is a brief implementation of an Gain model.
% 
% Function:
%   It can amplify the input signal based on the gain parameter without noise.
%
% Examples:
%   @code
%   sig = rand(1e3, 2);
%   param.sig = struct('Fs', 20e9, 'Fc', const.c/1550e-9, 'Rs', 20e9, 'PCol', '[pwr(20, {-16, 'dBm}), pwr(25, {-18, 'dBm'})]);
%   sigfield = signal_interface(sig, param.sig)
%   param.gain.gain = 16;
%   param.gain.gaintype = 'dB';
%   gain = Gain(param.gain);
%   edfafield = gain.generate(sigfield);
%   @endcode
%
% from @OCG-LAB OPC Framework
%
% 简体中文版注释：
% 类 Gain，文件名为 'Gain.m'
%   该类是一个简要的Gain模型实现。
% 
% 功能：
%   它可以根据增益参数对输入的信号进行无噪声放大。
%
% 示例：
%   @code
%   sig = rand(1e3, 2);
%   param.sig = struct('Fs', 20e9, 'Fc', const.c/1550e-9, 'Rs', 20e9, 'PCol', '[pwr(20, {-16, 'dBm'}), pwr(25, {-18, 'dBm'})]');
%   sigfield = signal_interface(sig, param.sig)
%   param.gain.gain = 16;
%   param.gain.gaintype = 'dB';
%   gain = Gain(param.gain);
%   edfafield = gain.generate(sigfield);
%   @endcode
%
% 来自 @OCG-LAB OPC 框架
classdef Gain < unit 

    properties
        nInputs = 1;
        nOutputs = 1;
        gain;           % Gain
        gaintype;       % gaintype, default:'dB'
    end

    methods
        function obj = Gain(param)
            if isfield(param, 'gain') && isscalar(param.gain), obj.gain = param.gain; else, slog('Please input "gain"', 'ERR'); end
            if isfield(param, 'gaintype') && ischar(param.gaintype), obj.gaintype = param.gaintype; else, obj.gaintype = 'dB'; end
        end

        function out = generate(obj, in)
            insig = in.get;
            switch obj.gaintype
                case 'dB'
                    Glin = 10^(obj.gain/10);
                case 'lin'
                    Glin = obj.gain;
                    obj.gain = 10 * log10(obj.gain);
                otherwise
                    slog('Please input the data type of gain "gaintype", it should be "lin" or "dB".', 'ERR');
            end
            insig = sqrt(Glin) * insig;
            out = signal_interface(insig, struct('Fs',in.Fs,'Fc',in.Fc,'Rs',in.Rs,'P', pwr(in.P.SNR, in.P.Ptot + obj.gain)));
        end
    end

end