
classdef setSNR < unit 

    properties
        nInputs = 1;
        nOutputs = 1;
        SNR;
        M;
    end

    methods
        function obj = setSNR(param)
            if isfield(param, 'SNR') && isscalar(param.SNR), obj.SNR = param.SNR; else, slog('Please input SNR', 'ERR'); end
            if isfield(param, 'M') && isscalar(param.M), obj.M = param.M; else, slog('Please input M', 'ERR'); end 
        end

        function out = generate(obj, in)
            SNRo = obj.SNR;
            if obj.SNR > in.P.SNR
                slog('SNR can only be decreased', 'ERR');
            elseif in.P.SNR < inf
                slog('Assume the input signal might be noisy (Use SNR information)');
                obj.SNR = 10*log10(in.P.SNR('lin')*(10^(obj.SNR/10))/(in.P.SNR('lin')-10^(obj.SNR/10)));
            end
            in = in.fun1(@(x) obj.addNoise(x, in.Rs, in.Fs));
            PCol = in.PCol;
            for idx = 1:in.N
                Pn = in.PCol(idx).Ptot - obj.SNR;
                PCol(idx) = pwr(SNRo, {in.PCol(idx).Ptot('W')+1e-3*10^(Pn/10), 'W'});
            end
            out = set(in, 'PCol', PCol);
        end

        function out = addNoise(obj, in, Rs, Fs)
            EbNo = 10^(obj.SNR/10);
            Es = sum(abs(in).^2) / length(in) / Rs;
            Eb = Es / log2(obj.M);
            No = Eb / EbNo;
            pn = No * Fs / 2;
            if isreal(in)
                noise = sqrt(pn) * randn(1, length(in))';
            else
                noise = sqrt(pn) * (randn(1, length(in)) + 1j*randn(1, length(in)))';
            end
            out = in + noise;
        end
    end

end