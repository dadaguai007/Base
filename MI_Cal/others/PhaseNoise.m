

classdef PhaseNoise

    properties
        Lpsd;
        Llw;
        fMin;
        fMax;
        linewidth;
        linewidth3dB;
        LFLW1GHz;
        HFLW;
        fr;
        K;
        alpha;
        type;
        fcutoff;

    end

    methods
        function obj = PhaseNoise(param)
            if isfield(param, 'Llw') && isscalar(param.Llw), obj.Llw = 2^14; end
            if isfield(param, 'fMin') && isscalar(param.fMin), obj.fMin = param.fMin; else, obj.fMin = 1e3; end
            if isfield(param, 'fMax') && isscalar(param.fMax), obj.fMax = param.fMax; else, obj.fMax = 50e9; end 
            if isfield(param, 'type') && ischar(param.type), obj.type = param.type; else, obj.type = 'linear'; end
            if isfield(param, 'fcutoff') && isscalar(param.fcutoff), obj.fcutoff = param.fcutoff; else, obj.fcutoff = 0; end 
            
            if isfield(param, 'linewidth')
                slog('Using Lorentzian mode.', 'BSCIF');
                if isscalar(param.linewidth)
                    obj.linewidth = param.linewidth;
                else
                    obj.linewidth = 100e3;
                end
                obj.Lpsd = 1;
                obj.linewidth3dB = 2 * obj.linewidth;
            else
                if isfield(param, 'Lpsd') && isscalar(param.Lpsd)
                    obj.Lpsd = param.Lpsd; 
                else
                    obj.Lpsd = 2^11;
                end
                judge = paramdetect(param, 'scale', 'LFLW1GHz', 'HFLW', 'fr', 'K', 'alpha');
                if prod(judge)
                    obj.LFLW1GHz = param.LFLW1GHz;
                    obj.HFLW     = param.HFLW;
                    obj.fr       = param.fr;
                    obj.K        = param.K;
                    obj.alpha    = param.alpha;
                else
                    for idx = 1:length(judge)
                        if ~judge(idx)
                            slog('Please input "%s".', 'ERR', paramstruct(idx));
                        end
                    end
                end
                if obj.LFLW1GHz > 50
                    obj.linewidth3dB = 1e5 * sqrt(obj.LFLW1GHz);
                else
                    obj.linewidth3dB = 2 * obj.HFLW;
                end
            end

            if isfield(param, 'linewidth3dB') && isscalar(param.linewidth3dB)
                obj.linewidth3dB = param.linewidth3dB;
            end

        end

        function FMfreq = FMfreq(obj)
            switch obj.type
                case 'linear'
                    FMfreq = linspace(obj.fMin, obj.fMax, obj.Lpsd);
                case 'log'
                    FMfreq = logspace(log10(obj.fMin), log10(obj.fMax), obj.Lpsd);
                otherwise
                    slog('Incorrect "type" input, only support two kinds, "linear" and "log".', 'ERR');
            end
        end

        function [frequencyModulationNoise, phaseModulationNoise] = generatePSD(obj)
            FMfreq = obj.FMfreq();
            if ~isempty(obj.linewidth)
                frequencyModulationNoise = obj.linewidth / pi * ones(size(FMfreq));
                % frequencyModulationNoise = 2*pi*obj.linewidth*randn(obj.Llw, 1);
            else
                frequencyModulationNoise = obj.LFLW1GHz / pi * 1e9 ./ FMfreq + obj.HFLW / pi * (1/(1+obj.alpha^2)) + ...
                obj.HFLW / pi * (obj.alpha^2/(1+obj.alpha^2)) * obj.fr^4 ./ ((obj.fr^2-FMfreq.^2).^2+(obj.K/2/pi)^2*obj.fr^4*FMfreq.^2);
            end
            if obj.fcutoff ~= 0 && obj.fcutoff >= obj.fMin
                frequencyModulationNoise = frequencyModulationNoise .* (FMfreq < obj.fcutoff);
            end
            phaseModulationNoise = (2*pi)^2 * frequencyModulationNoise ./ ((2*pi)^2*FMfreq.^2);
        end

        function effectiveLinewidth = geteffLinewidth(symbolRate, FMnoise, FMfreq)
            effectiveLinewidth = 2*pi*(1/symbolRate)*trapz(FMfreq, FMnoise.*sinc(FMfreq/symbolRate).^2);
        end
    end

end