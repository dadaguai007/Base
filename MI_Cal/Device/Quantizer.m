classdef Quantizer < unit

	properties
		nInputs = 1;
		nOutputs = 1;
		bitResolution;
		ENoB;
		quantizationType % riser or tread
	end

	methods

		function obj = Quantizer(varargin)
			if nargin == 1 && ~isempty(fieldnames(varargin{1}))
                param = varargin{1}; 
				obj.bitResolution = paramdefault(param, 'bitResolution', 8);
				obj.ENoB = paramdefault(param, 'ENoB', 16);
				obj.quantizationType = paramdefault(param, 'quantizationType', 'riser');
			elseif nargin == 0 || isempty(fieldnames(varargin{1}))
				obj.bitResolution = 8;
				obj.ENoB = 16;
				obj.quantizationType = 'riser';
			else
				slog('Too many input arguments', 'ERR');
			end
		end

		function out = quantization(obj, in)
			if ~isreal(in)
				step = [(max(real(in))-min(real(in))) / round(2^obj.bitResolution-1), (max(imag(in))-min(imag(in))) / round(2^obj.bitResolution-1)];
				switch lower(obj.quantizationType)
					case 'riser'
						out = step(1)*(round(real(in)/step(1))+0.5) + 1j*step(2)*(round(imag(in)/step(2))+0.5);
					case 'tread'
						out = step(1)*sign(real(in)).*(round(abs(real(in))/step(1)+0.5)-0.5) + 1j*step(2)*sign(imag(in)).*(round(abs(imag(in))/step(2)+0.5)-0.5);
					otherwise
						slog('The quantizer only supports two quantization methods: riser and tread.', 'ERR');
				end
			else
				step = (max(in)-min(in)) / round(2^obj.bitResolution-1);
				switch lower(obj.quantizationType)
					case 'riser'
						out = step * (round(in/step)+0.5);
					case 'tread'
						out = step * sign(in).*(round(abs(in)/step+0.5)-0.5);
					otherwise
						slog('The quantizer only supports two quantization methods: riser and tread.', 'ERR');
				end
			end
		end

		function out = addnoise(obj, raw, in)
			rawsig = raw.get;
			insig = in.get;
			outsig = insig;
			p = in.params;
			distortionpwr = pwr.meanpwr(rawsig-insig);
			inputpwr = pwr.meanpwr(rawsig);
			noisepwr = inputpwr / (3*2^(2*obj.ENoB-1)-1) - distortionpwr;
			% distortionpwr = mean(2^(obj.bitResolution-1)*rawsig(:, idx)-insig(:, idx).^2);
			% inputpwr = mean(2^(obj.bitResolution-1)*rawsig(:,idx).^2);
			% noisepwr = inputpwr / (3*2^(2*obj.ENoB-1)-1) - distortionpwr;
			noisepwr(noisepwr<0) = 0;
			outsig = outsig + sqrt(noisepwr).*randn(size(insig));
			for idx = 1:in.N
				p.PCol(idx) = pwr(10*log10(inputpwr(idx)./(noisepwr(idx)+distortionpwr(idx))), {inputpwr(idx), 'W'});
			end
			out = signal_interface(outsig, p);
			out.set('P', out.setPtot);
		end

		function out = generate(obj, in)
			out = in.fun1(@(x) obj.quantization(x));
			out = obj.addnoise(in, out);
		end
    end
end