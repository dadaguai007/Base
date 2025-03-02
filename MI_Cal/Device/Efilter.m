
classdef Efilter < unit

    properties
        nInputs = 1;
        nOutputs = 1;
        filterType;
        filterOrder;
        filterBandwidth;
        outputVoltage;
        levelDC;
    end

    methods

        function obj = Efilter(param)
            if isfield(param, 'filterType') && ischar(param.filterType), obj.filterType = param.filterType; else, obj.filterType = 'gauss'; end
            if isfield(param, 'filterOrder') && isscalar(param.filterOrder), obj.filterOrder = param.filterOrder; else, obj.filterOrder = 0; end
            if isfield(param, 'filterBandwidth') && isscalar(param.filterBandwidth), obj.filterBandwidth = param.filterBandwidth; else, slog('Please input the bandwidth of electrical filter "filterBandwidth"', 'ERR'); end
            if isfield(param, 'outputVoltage') && isscalar(param.outputVoltage), obj.outputVoltage = param.outputVoltage; else, obj.outputVoltage = 1; end
            if isfield(param, 'levelDC') && isscalar(param.levelDC), obj.levelDC = param.levelDC; else, obj.levelDC = 0; end
        end

        function out = generate(obj, in)
            if obj.filterOrder
                if strcmp(obj.filterType, 'window')
                    slog('Rectanguler filter has no order, please do not input filterOrder', 'ERR');
                end
            else
                if ~strcmp(obj.filterType, 'window')
                    slog('Please input filterOrder', 'ERR');
                end
            end
            switch lower(obj.filterType)
                case 'window'
                    filtercoeftime = (2*obj.filterBandwidth/(in.Fs))*sinc(2*in.SPS*(obj.filterBandwidth/(in.Fs))*linspace(-floor(in.L/2-1)/(in.SPS), floor(in.L/2-1)/(in.SPS), floor(in.L/2-1)*2+1)).';
                    coefshift = -floor(length(filtercoeftime)/2);
                    filtercoeftime(in.L) = 0;
                    filtercoeffreq = fftshift(fft(circshift(filtercoeftime(:), [coefshift 0])));
                case 'gauss'
                    filtercoeffreq = exp(-log(sqrt(2))*((linspace(-0.5,0.5,in.L)/(obj.filterBandwidth/(in.Fs))).').^(2*obj.filterOrder));
                case 'bessel'
                    [b,a] = besself(obj.filterOrder, obj.filterBandwidth);
                    besselcoef = polyval(b, in.Fs*2j*pi*linspace(-0.5,0.5,in.L))./polyval(a, in.Fs*2j*pi*linspace(-0.5,0.5,in.L));
                    filtercoeffreq = exp(1j*angle(besselcoef(:)));
            end
            sig = in.get;
            for idx = 1:in.N
                sig(:, idx) = fftshift(fft(sig(:, idx)));
                sig(:, idx) = ifft(ifftshift(sig(:, idx).*filtercoeffreq(:)));
            end
            Pout = pwr.meanpwr(sig);
            for idx = 1:length(in.PCol)
                PCol_new(idx) = pwr(in.PCol(idx).SNR, {Pout(idx), 'W'});
            end
            in = signal_interface(sig, struct('Rs', in.Rs, 'Fs', in.Fs, 'Fc', in.Fc, 'PCol', PCol_new));
            sig = in.getRaw;
            for idx = 1:size(sig,2)
                if obj.outputVoltage
                    maxVoltage = max(max(abs([real(sig) imag(sig)])));
                    sig(:,idx) = (real(sig(:,idx)) + 1j*imag(sig(:,idx)))*obj.outputVoltage/maxVoltage;
                else
                    sig(:,idx) = (real(sig(:,idx)) + 1j*imag(sig(:,idx)));
                end
                sig(:,idx) = sig(:,idx) + obj.levelDC + 1j*obj.levelDC;
            end

            P = pwr.meanpwr(sig);
            for idx = length(P):-1:1
                Power{idx} = pwr(in.PCol(idx).SNR, {P(idx), 'w'});
            end
            out = signal_interface(sig, struct('Rs', in.Rs, 'Fs', in.Fs, 'Fc', in.Fc, 'PCol', [Power{:}]));

        end
    end
end