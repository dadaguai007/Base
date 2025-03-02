

classdef Laser < unit 

    properties
        nInputs;
        nOutputs = 1;
        carrierFreq;
        symbolRate;
        samplingRate;
        Lnoise;
        outputpwr;
        phaseNoiseModelobj;
        phaseNoiseLimitpwr;
        frequencyModulationNoise;
        phaseModulationNoise;
        frequencyNoise;
        phaseNoise;
        L;
        idealEnabled;
    end

    methods
        function obj = Laser(param)
            if isfield(param, 'carrierFreq') && isscalar(param.carrierFreq), obj.carrierFreq = param.carrierFreq; else, obj.carrierFreq = const.c / 1550e-9; end 
            if isfield(param, 'symbolRate') && isscalar(param.symbolRate), obj.symbolRate = param.symbolRate; else, obj.symbolRate = 1; end
            if isfield(param, 'outputpwr') && isa(param.outputpwr, 'pwr'), obj.outputpwr = param.outputpwr; else, obj.outputpwr = pwr(inf, 0); end
            if isfield(param, 'idealEnabled') && isscalar(param.idealEnabled), obj.idealEnabled = param.idealEnabled; else, obj.idealEnabled = 1; end
            phaseNoiseModelparams = paramsref('PhaseNoise', param);
            phaseNoiseModelparams.type = 'linear';
            obj.phaseNoiseModelobj = PhaseNoise(phaseNoiseModelparams);
            
            if isfield(param, 'phaseNoiseLimitpwr') && isscalar(param.phaseNoiseLimitpwr), obj.phaseNoiseLimitpwr = param.phaseNoiseLimitpwr; else, obj.phaseNoiseLimitpwr = 0; end
            
            judge = paramdetect(param, 'scale', 'samplingRate', 'Lnoise');
            if prod(judge)
                obj.samplingRate = param.samplingRate;
                obj.Lnoise = param.Lnoise;
                obj.nInputs = 0;
                slog('It is not allowed to input laser signal.', 'BSCIF');
            else
                obj.nInputs = 1;
                slog('Laser parameters will be copied from input signal.', 'BSCIF');
            end
            
            if isfield(param, 'L') && isscalar(param.L), obj.L = param.L; else, obj.L = 2^11; end
                

        end

        function out = generate(obj, varargin)
            if obj.nInputs == 1 && isempty(varargin)
                slog('Missing input signal from where to copy parameters', 'ERR');
            elseif obj.nInputs == 0 && ~isempty(varargin)
                slog('Too many input arguments. Laser parameters are already set.', 'ERR');
            end

            if obj.nInputs == 1
                obj.samplingRate = varargin{1}.samplingRate;
                obj.symbolRate = varargin{1}.symbolRate;
                obj.Lnoise = varargin{1}.Lnoise;
                obj.phaseNoiseModelobj.fMin = obj.samplingRate / 2;
                obj.phaseNoiseModelobj.fMax = obj.samplingRate / obj.phaseNoiseModelobj.Lpsd;
            end

            if ~obj.idealEnabled
                [obj.frequencyNoise, obj.phaseNoise] = calculate(obj);
            else
                obj.frequencyNoise = zeros(obj.Lnoise,1);
                obj.phaseNoise = zeros(obj.Lnoise,1);
            end

%             t = (1/obj.samplingRate : 1/obj.samplingRate : obj.Lnoise / obj.samplingRate).';
            laserwave = sqrt(obj.outputpwr.Ptot('W')) .* exp(1j*obj.phaseNoise(1:obj.Lnoise));
            out = signal_interface(laserwave(:), struct('Fs', obj.samplingRate, 'Rs', obj.symbolRate, 'P', obj.outputpwr, 'Fc', obj.carrierFreq));

        end

        function [frequencyNoise, phaseNoise] = calculate(obj)
            obj.Lnoise = (obj.Lnoise + 2*obj.phaseNoiseModelobj.Lpsd) / 2;
            [obj.frequencyModulationNoise, obj.phaseModulationNoise] = obj.phaseNoiseModelobj.generatePSD();
            f = [-flipud(obj.phaseNoiseModelobj.FMfreq); 0; obj.phaseNoiseModelobj.FMfreq];
            t = linspace(-obj.phaseNoiseModelobj.Lpsd, obj.phaseNoiseModelobj.Lpsd, 2*obj.phaseNoiseModelobj.Lpsd+1) / obj.samplingRate;
            Hfn = obj.frequencyModulationNoise(:);
            if length(Hfn) > 1
                Hfn = sqrt([flipud(Hfn(:)); 0; Hfn(:)] / 2);
                hfn = real(obj.idft(f, Hfn, t));
            else
                hfn = obj.samplingRate * sqrt(Hfn/2);
            end
            ase = awgn(zeros(2*obj.Lnoise, 1), 0) / sqrt(obj.samplingRate);
            frequencyNoise = conv(ase, hfn, 'valid');
            phaseNoise = cumsum(2*pi*frequencyNoise/obj.samplingRate);
            if obj.phaseNoiseLimitpwr
                limit = [-1 1] * obj.phaseNoiseLimitpwr;
                phaseNoise = phaseNoise - mean(phaseNoise);
                phaseNoise(phaseNoise>max(limit)) = -phaseNoise(phaseNoise>max(limit)) + 2*max(limit);
                phaseNoise(phaseNoise<min(limit)) = -phaseNoise(phaseNoise<min(limit)) + 2*min(limit);
            end
            obj.Lnoise = 2 * (obj.Lnoise - obj.phaseNoiseModelobj.Lpsd);
        end

    end

    methods(Static)
        function y = idft(f,x,t)
            shape = size(t);
            f = f(:);
            x = x(:);
            t = t(:);

            df = diff(f);
            df(length(df)/2) = 0;
            df=([0;df]+[df;0]) / 2;

            y(1:length(t), 1) = 0;
            for idx = 1:length(t)
                y(idx) = exp(2*pi*1j*t(idx)*f') * (x.*df);
            end
            y = reshape(y, shape);
        end
    end

end