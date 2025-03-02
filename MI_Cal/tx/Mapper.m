classdef Mapper < unit
    
    properties
        nInputs = 1;
        nOutputs = 1;
        symbolRate;
        M;
        modulationFormat;
        codeFormat;
        N;
    end

    methods
        function obj = Mapper(param)
            if isfield(param, 'symbolRate') && isscalar(param.symbolRate), obj.symbolRate = param.symbolRate; else, obj.symbolRate = 1; end
            if isfield(param, 'M') && isscalar(param.M), obj.M = param.M; else, error('Please input modulation order "M".'); end
            if isfield(param, 'modulationFormat') && ischar(param.modulationFormat), obj.modulationFormat = lower(param.modulationFormat); else, error('Please input "modulationFormat".'); end
            if isfield(param, 'codeFormat') && ischar(param.codeFormat), obj.codeFormat = lower(param.codeFormat); else, obj.codeFormat = 'bin'; end
            if isfield(param, 'N') && isscalar(param.N), obj.N = param.N; else, obj.N = 1; end

            if sum(mod(log2(obj.M),1) ~= 0)
                slog('Modulation Order, M, should be a power of 2.', 'ERR');
            end
        end

        function out = generate(obj, in)
            insigRaw = in.getRaw;
            h = nextpow2(obj.M);
%                 xpolmap = binmap(insigRaw(:, 1:h));
%                 ypolmap = binmap(insigRaw(:, h+1:end));
%                 out = [xpolmap ypolmap];
%             else
%                 out = binmap(insigRaw);
%             end
            out = zeros(in.L, obj.N);
            for i = 1:obj.N
                out(:, i) = binmap(insigRaw(:, 1+(i-1)*h:i*h));
            end
            switch obj.codeFormat
                case 'gray'
                    out = bin2gray(out);
                case 'bin'
                otherwise
                    slog('This framework only supports "gray" and "bin" code format so far.', 'ERR');
            end
            constellation = constellationref(obj.modulationFormat, obj.M);
            out = constellation(out+1);
            out = signal_interface(out, struct('Fs', obj.symbolRate, 'Rs', obj.symbolRate, 'Fc', in.Fc));
            out = out.set('P', pwr(inf, {1, 'W'}));
        end
    end

end