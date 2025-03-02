classdef PulseShaper < unit

	properties (Access = public)
		%> Number of inputs
		nInputs = 1;
		%> Number of outputs
        nOutputs = 1;
        %> Output samples per symbol
        samplesPerSymbol;
        %> Output sample rate
        symbolRate = 1;
        %> shaper type;
		filterType;
		%> filter bandwidth limit
		bt;
		%> filter span
		span;
    end

	methods

		function obj = PulseShaper(param)
            if isfield(param, 'bt') && isscalar(param.bt), obj.bt = param.bt; else, slog('Please input bandwidth limitation "bt"', 'ERR'); end
			if isfield(param, 'span') && isscalar(param.span), obj.span = param.span; else, slog('Please input the number of symbols "span"', 'ERR'); end
			if isfield(param, 'filterType') && ischar(param.filterType), obj.filterType = param.filterType; else, obj.filterType = 'gauss'; end
			if isfield(param, 'samplesPerSymbol') && isscalar(param.samplesPerSymbol), obj.samplesPerSymbol = param.samplesPerSymbol; else, slog('Please input samples per symbol "samplesPerSymbol".', 'ERR'); end
			if isfield(param, 'symbolRate') && isscalar(param.symbolRate), obj.symbolRate = param.symbolRate; else, obj.symbolRate = 1; end
		end

		function out = generate(obj, in)
			insig = in.get;
			switch lower(obj.filterType)
				case 'gauss'
					insig = upsample(insig, obj.samplesPerSymbol);
					ch = gaussdesign(obj.bt, obj.span, obj.samplesPerSymbol);
					out = conv(insig, ch, 'same');
				case 'rrc'
					insig = upsample(insig, obj.samplesPerSymbol);
					ch = rcosdesign(obj.bt, obj.span, obj.samplesPerSymbol, 'normal');
                    out = conv(insig.', ch, 'same').';
				case 'square'
					insig = repelem(insig, obj.samplesPerSymbol, 1);
					ch = fir1(obj.samplesPerSymbol*obj.span, obj.bt, 'low');
					out = conv(insig, ch, 'same');
				otherwise
					slog('No filter type is found', 'ERR');
            end
			slog('Pulse shaping is done, shape format: %s, filter length: %d', obj.filterType, obj.span * obj.samplesPerSymbol);
			outparams.Rs = max(obj.symbolRate, in.Rs);
			outparams.Fs = outparams.Rs * obj.samplesPerSymbol;
            for i = 1:in.N
                outparams.PCol(i) = pwr(inf, {pwr.meanpwr(out(:, i)), 'W'});
            end
			out = signal_interface(out, outparams);
		end

	end

end