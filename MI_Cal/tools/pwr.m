classdef pwr

    properties(SetAccess=immutable, Hidden=true)
        P_dBW;
        SNR_dB;
    end

    methods(Access=private, Hidden=true, Static=true)
        function res = validpwr(P)
            res = false;
            if isscalar(P) && isreal(P), res = true; end
        end

        function P_format = outputP(P_dBW,varargin)
            type = defaultunit('dBm',varargin{:});
            switch lower(type)
                case 'w'
                    P_format = dB2lin(P_dBW);
                case 'mw'
                    P_format = dB2lin(P_dBW)*1e3;
                case 'uw'
                    P_format = dB2lin(P_dBW)*1e6;
                case 'dbw'
                    P_format = P_dBW;
                case 'dbm'
                    P_format = P_dBW+30;
                case 'dbu'
                    P_format = P_dBW+60;
                otherwise
                    error('Power unit can be W, mW, uW, dBW, dBm or dBu.');
            end
        end

        function SNR_format = outputSNR(SNR_dB,varargin)
            type = defaultunit('dB',varargin{:});
            switch lower(type)
                case 'lin'
                    SNR_format = dB2lin(SNR_dB);
                case 'db'
                    SNR_format = SNR_dB;
                otherwise
                    error('SNR unit can be lin for linear or dB.');
            end
        end

        function P_dB = inputP(P)
            if pwr.validpwr(P)
                type = 'dBm';
            elseif iscell(P) && numel(P)==2 && pwr.validpwr(P{1}) && ischar(P{2})
                type = lower(P{2});
                P = P{1};
            else
                error('Power must be specified as a real scalar (in dBm) or a cell array must be constructed as {value,unit}');
            end
            switch lower(type)
                case 'w'
                    P_dB = lin2dB(P);
                case 'mw'
                    P_dB = lin2dB(P/1e3);
                case 'uw'
                    P_dB = lin2dB(P/1e6);
                case 'dbw'
                    P_dB = P;
                case 'dbm'
                    P_dB = P-30;
                case 'dbu'
                    P_dB = P-60;
                otherwise
                    error('Power unit can be W, mW, uW, dBW, dBm or dBu.');
            end
        end
        
        function SNR_dB = inputSNR(SNR)
            if pwr.validpwr(SNR)
                type = 'dB';
            elseif iscell(SNR) && numel(SNR)==2 && pwr.validpwr(SNR{1}) && ischar(SNR{2})
                type = lower(SNR{2});
                SNR = SNR{1};
            else
                error('SNR must be specified as a real scalar (in dB) or a cell array must be constructed as {value,unit}.');
            end
            switch lower(type)
                case 'lin'
                    SNR_dB = lin2dB(SNR);
                case 'db'
                    SNR_dB = SNR;
                otherwise
                    error('SNR unit can be lin for linear or dB.');
            end
        end
    end

    methods (Static=true)
        function [y,scfactor] = normpwr(x,varargin)
            [type,P,unit] = defaultargs({'average',1,'linear'},varargin);

            switch lower(unit)
                case {'linear','lin'}
                    P_ = P;

                case {'db'}
                    P_ = 10^(P/10);

                case {'dbm'}
                    P_ = 1e-3*10^(P/10);

                otherwise
                    error('Unit must be ''linear'', ''dB'' or ''dBm''.');
            end
    
            % TODO Handling real and complex cases
            [Pav,~,Prange] = pwr.meanpwr(x(:));
            switch lower(type)
                case {'average','avg'}
                    scfactor = Pav/P_;

                case {'maximum','max'}
                    scfactor = Prange(2,:)/P_;

                otherwise
                    error('Normalization type must be ''average'' or ''maximum''.');
            end

%             y = bsxfun(@mrdivide,x,sqrt(scfactor));
            y = x/sqrt(scfactor);
        end

        function [P,E,Prange,Ppeak] = meanpwr(x)
            
            if isvector(x)
                L = numel(x);
            else
                L = size(x,1);
            end
            %absxsq = abs(x).^2;
            absxsq = x.*conj(x);        %faster
            Ppeak = max(absxsq);
            E = sum(absxsq);
%             if isreal(x)
%                 E = E/2;
%             end
            P = E/L;
            Prange = findrange(absxsq);
        end

        function OSNR = getOSNR(sig, varargin)
            if nargin < 3
                NBW = 0.1;
            else
                NBW = varargin{2};
            end
            if sig.Fc == 0
                slog('The signal carrier frequency is required to compute the noise bandwidth', 'ERR');
            end
            lambda = const.c/sig.Fc;
            NBW_Hz = const.c./(lambda-.5*1e-9*NBW)-const.c./(lambda+.5*1e-9*NBW);
            OSNR_dB = sig.P.SNR+10*log10(sig.Fs/NBW_Hz);
            if nargin > 1
                OSNR = pwr.outputSNR(OSNR_dB,varargin(1));
            else
                OSNR = pwr.outputSNR(OSNR_dB,{});
            end
        end
    end
    
    methods 
        function obj = pwr(SNR,P)
            if nargin<1
                error('Signal-to-noise ratio must be specified.');
            elseif nargin<2
                P = 0;
            end
            obj.SNR_dB = pwr.inputSNR(SNR);
            obj.P_dBW = pwr.inputP(P);
        end

        function obj = plus_1elem(obj1,obj2)
            Ps = obj1.Ps('W')+obj2.Ps('W');
            Pn = obj1.Pn('W')+obj2.Pn('W');
            %zero total power is a special case
            if isinf(log10(Ps))&&isinf(log10(Pn))
                obj = pwr(0, -Inf);
            else
                obj = pwr({Ps/Pn,'lin'},{Ps+Pn,'W'});
            end
        end

        function obj = plus(obj1, obj2)
            if numel(obj1)~=numel(obj2)
                error('When adding two power object arrays, number of elements in array 1 must equal number of elements in array 2');
            end
            obj = plus_1elem(obj1(1), obj2(1));
            for jj=2:numel(obj1)
                obj(jj) = plus_1elem(obj1(jj), obj2(jj));
            end
        end

        function obj = minus(obj1,obj2)
            warning('Subtraction is equivalent to addition -- no correlation.')
            obj = plus(obj1,obj2);
        end

        function obj = mtimes(in1,in2)
            if isa(in1, 'pwr')&&~isa(in2, 'pwr')
                Ps = dB2lin([in1.P_dBW])*in2;
                SNR = [in1.SNR_dB];
            elseif ~isa(in1, 'pwr')&&isa(in2, 'pwr')
                Ps = in1*dB2lin([in2.P_dBW]);
                SNR = [in2.SNR_dB];
            end
            for jj=1:numel(Ps)
                obj(jj) = pwr(SNR(jj), {Ps(jj), 'W'});
            end
        end

        function obj = times(in1, in2)
            if numel(in1)==numel(in2)
                for jj=1:numel(in1), obj(jj) = in1(jj)*in2(jj); end
            elseif (numel(in1)==1)||(numel(in2)==1)
                if numel(in1)==1
                    for jj=1:numel(in2), obj(jj) = in1*in2(jj); end
                else
                    for jj=1:numel(in1), obj(jj) = in1(jj)*in2; end
                end
            else
                error('Matrix dimensions must agree');
            end
        end

        function obj = mrdivide(in1, C)
            obj = mtimes(in1, 1/C);
        end

        function P = P(obj,varargin)
            warning('Please use Ptot instead of P function of pwr.');
            P = Ptot(obj,varargin{:});
        end

        function SNR = SNR(obj,varargin)
            SNR = pwr.outputSNR(obj.SNR_dB,varargin);
        end

        function Ptot = Ptot(obj,varargin)
            Ptot = pwr.outputP(obj.P_dBW,varargin);
        end

        function Ps = Ps(obj,varargin)
            if isinf(obj.SNR('lin')) % If SNR is infinite, there is no noise, so Ps = P;
                Ps = obj.Ptot('dBW');
            else
                Ps = obj.Ptot('dBW')+obj.SNR('dB')-10*log10(obj.SNR('lin')+1);
            end
            Ps = pwr.outputP(Ps,varargin);
        end

        function Pn = Pn(obj,varargin)
            Pn = pwr.outputP(obj.Ptot('dBW')-10*log10(obj.SNR('lin')+1),varargin);
        end
    end

end