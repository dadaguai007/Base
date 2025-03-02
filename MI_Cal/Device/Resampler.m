classdef Resampler < unit

    properties
        nInputs = 1;
        nOutputs = 1;
        resampleMethod;
        outputSamplingRate;
    end

    methods
        function obj = Resampler(param)
            if isfield(param, 'resampleMethod') && ischar(param.resampleMethod), obj.resampleMethod = param.resampleMethod; else, obj.resampleMethod = 'matlab'; end 
            if isfield(param, 'outputSamplingRate') && isscalar(param.outputSamplingRate), obj.outputSamplingRate = param.outputSamplingRate; else, slog('Please input "outputSamplingRate".', 'ERR'); end 
        end

        function out = generate(obj, in)
            out = in.get;
            switch lower(obj.resampleMethod)
                case 'matlab'
                    out = resample(out, obj.outputSamplingRate, in.Fs);
                case 'spline'
                    for idx = 1:size(out, 2)
                        temp = splineResampler(out(:, idx), obj.outputSamplingRate, in.Fs);
                        out(:, idx) = temp;
                    end
                otherwise
                    slog('This class only supports "matlab" resampler and "spline" resampler', 'ERR');
            end
            inparams = in.params;
            inparams.Fs = obj.outputSamplingRate;
            out = signal_interface(out, inparams);
        end

    end

end