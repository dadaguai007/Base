classdef setOSNR < unit 

    properties
        OSNR;
        Nb = 0.1;  % Noise bandwidth (nm)
        nInputs = 1;
        nOutputs = 1;
    end

    methods

        function obj = setOSNR(param)
            if isfield(param, 'OSNR') && isscalar(param.OSNR), obj.OSNR = param.OSNR; else, slog('Please input "OSNR".', 'ERR'); end
            if isfield(param, 'Nb') && isscalar(param.Nb), obj.Nb = param.Nb; else, obj.Nb = 0.1; end
        end

        function out = generate(obj, in)
            lambda = const.c / in.Fc;
            NbHz = const.c/(lambda-.5*1e-9*obj.Nb)-const.c/(lambda+.5*1e-9*obj.Nb);

            Ps = in.P.Ps('W');
            Pn = in.P.Pn('W');
            currentOSNR = pwr.getOSNR(in);
            if obj.OSNR > currentOSNR
                slog('The input optical signal-to-noise ratio (OSNR) is already greater than the signal optical signal-to-noise ratio, so it cannot be set.', 'ERR');
            end

            OSNRlin = 10^(obj.OSNR/10);
            N0 = Ps / OSNRlin / NbHz;
            Pn = N0*in.Fs - Pn;
            PnPerCol = Pn / in.N;

            p = struct( ...
                'Fs', in.Fs, ...
                'Fc', in.Fc, ...
                'Rs', in.Rs, ...
                'P', pwr(-inf, {Pn, 'W'}) ...
                );
            
            noise = wgn(in.L, in.N, PnPerCol, 'linear', 'complex');
            noise = signal_interface(noise, p);

            out = in + noise;
        end

    end

end