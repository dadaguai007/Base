% class signal_interface, filename 'signal_interface.m'
%   It is a brief signal interface class definition
%
% Basic information:
%   signal_interface class contains some main parameters of signal.
%   1.waveform, the sig.
%   2.carrier frequency, Fc. 
%   3.sampling rate, Fs.
%   4.symbol rate, Rs.
%   5.power, P or PCol.
%   
%   signal_interface class is a data class specification defined in this framework for passing data between units.
%   
% Function:
%   1.Stores signals along with their key data and calculates key parameters.
%   2.Acts as a signal carrier for passing between classes.
%   3.Perform operations on signals.
%
% Attention:
%   1.Fc can be 0, which represents the signal as a logic or electrical signal with no carrier frequency.
%   2.Fs must have a definite positive value. When it is 1, it indicates that the signal is a logic signal.
%   3.When Rs is not input, it defaults to 1, indicating that the signal is considered analog.
%
% Examples:
%   @code
%   sig = rand(1e3, 1);
%   param.sig.Fs = 10e9;
%   param.sig.Fc = 0;
%   param.sig.Rs = 5e9;
%   param.sig.P = pwr(20, 3);
%   sigfield = signal_interface(sig, param.sig);
%   @endcode
%
% from @OCG-LAB OPC Framework
%
% 简体中文版注释：
% signal_interface类，文件名'signal_interface.m'
% 这是一个简要的信号接口类定义
%
% 基本信息：
% signal_interface类包含信号的一些主要参数。
% 1.波形，即sig。
% 2.载波频率，即Fc。
% 3.采样率，即Fs。
% 4.符号速率，即Rs。
% 5.功率，即P或PCol。
%
% signal_interface类是在该框架中定义的数据类规范，用于在单元之间传递数据。
%
% 功能：
% 1.存储信号及其关键数据，并计算关键参数。
% 2.作为信号载体在类之间传递。
% 3.对信号执行操作。
%
% 注意事项：
% 1.Fc可以为0，表示信号为逻辑或电信号，没有载波频率。
% 2.Fs必须具有明确的正值。当它为1时，表示信号为逻辑信号。
% 3.当未输入Rs时，默认为1，表示信号被视为模拟信号。
%
% 示例：
% @code
% sig = rand(1e3, 1);
% param.sig.Fs = 10e9;
% param.sig.Fc = 0;
% param.sig.Rs = 5e9;
% param.sig.P = pwr(20, 3);
% sigfield = signal_interface(sig, param.sig);
% @endcode
%
% 来自@OCG-LAB OPC框架

classdef signal_interface
    properties(SetAccess=protected)
        Fc = 0;     % Carrier frequency (Hz)
        Fs;         % Sampling rate (Hz)
        Rs;         % symbol rate (Baud)
        P;          % Total signal power (pwr object)
        PCol;       % Signal power per column (pwr object)
    end

    properties(SetAccess=protected, Hidden=true)
        sig;        % waveform
    end

    methods
        function obj = signal_interface(signal, param)
            if isfield(param,'Fc') && isscalar(param.Fc), obj.Fc = param.Fc; else, obj.Fc = 1; end
            if isfield(param,'Fs') && isscalar(param.Fs), obj.Fs = param.Fs; else, slog('Sampling frequency must be specified','ERR'); end
            if isfield(param,'Rs') && isscalar(param.Rs), obj.Rs = param.Rs; else, obj.Rs = 1; end
            if isvector(signal)
                signal = signal(:);
            end
            obj.sig = signal;
            if isfield(param, 'PCol')
                if ~isa(param.PCol, 'pwr')
                    slog('Power should be a "pwr" object.', 'ERR')
                else
                    if numel(param.PCol) ~= size(obj.sig, 2)
                        slog('The number of signal inputs are not matched with the length of PCol', 'ERR');
                    end
                    obj.PCol = param.PCol;
                    obj.P = setPtot(obj);
                end
            elseif isfield(param, 'P')
                if ~isa(param.P, 'pwr')
                    slog('Power should be a "pwr" object.', ERR);
                end
                if length(param.P) > 1
                    slog('The size of "P" should not be larger than 1, please use "PCol" instead.', 'ERR');
                end
                obj.P = param.P;
                avgpwr = pwr.meanpwr(signal);
                pwrfraction = avgpwr / sum(avgpwr);
                obj.PCol = repmat(obj.P, obj.N, 1) .* pwrfraction;
            else
                slog('None "pwr" object is detected, power will be calculated automatically.', 'BSCIF');
                slog('SNR will be set in dB and Ptot will be set in dBm', 'BSCIF');
                avgpwr = pwr.meanpwr(signal);
                obj.PCol = pwr(inf, {avgpwr(1), 'W'});
                for idx = 2:obj.N
                    obj.PCol(idx) = pwr(inf, {avgpwr(idx), 'W'});
                end
                obj.P = setPtot(obj);
            end
        end

        function p = params(obj)
            p = struct( ...
                'Fs', obj.Fs, ...
                'Rs', obj.Rs, ...
                'Fc', obj.Fc, ...
                'P', obj.P, ...
                'PCol', obj.PCol ...
            );
        end

        function L = L(obj)
            L = size(obj.sig,1);
        end

        function N = N(obj)
            N = size(obj.sig,2);
        end

        function SPS = SPS(obj)
            % Samples per symbol
            SPS = obj.Fs/obj.Rs;
        end

        function Ts = Ts(obj)
            % Sampling time
            Ts = 1/obj.Fs;
        end

        function Pout = setPtot(obj)
            % Calculate total signal power
            Pout = obj.PCol(1);
            for idx = 2:numel(obj.PCol)
                Pout = Pout + obj.PCol(idx);
            end
        end

        function s = get(obj)
            % get waveform
            s = obj.sig;
            Pin = pwr.meanpwr(s);
            for idx = 1:obj.N
                Pout(idx) = obj.PCol(idx).Ptot('W');
            end
            if Pin ~= 0
                s = bsxfun(@times, s, sqrt(Pout./Pin));
            end
        end

        function s = getRaw(obj)
            % get raw waveform
            s = obj.sig;
        end

        function s = getNormalized(obj)
            % get normalized waveform
            s = obj.sig;
            s = s/sqrt(mean(pwr.meanpwr(s)));
            s = set(obj, s);
        end
        
        function obj = fun1(obj,fun)
            % Overloading fun1 can be used for operations on signal_interface objects.
            s = cellfun(fun,mat2cell(get(obj),obj.L,ones(1,obj.N)),'UniformOutput',false);
            if ~all(cellfun(@(c)isequal(size(c),size(s{1})),s))
                slog('Outputs after fun1 have different lengths. Cannot concatenate.','ERR');
            end
            obj.sig = cell2mat(s);
            Pout = pwr.meanpwr(obj.getRaw);
            for i=1:length(obj.PCol)
                PCol_new(i) = pwr(obj.PCol(i).SNR, {Pout(i), 'W'});
            end
            obj.PCol = PCol_new;
            obj.P = setPtot(obj);
        end

        function obj = mtimes(obj, M)
            % Overloading '×' can be used for operations on signal_interface objects.
            if isscalar(M)
                M = M * eye(obj.N);
            end
            if ~isequal(size(M),([obj.N obj.N]))
                slog('Jones matrix must be a NxN square matrix.','ERR');
            end

            Fnew=zeros(size(get(obj)));
            Pscale=zeros(obj.N, 1);
            for i=1:obj.N
                Fnew(:,i) = get(obj) * M(i,:).';
                Pscale(i) = norm(M(:,i), 2)^2;
            end
            obj.sig = Fnew;
            %power scaling
            obj.PCol = obj.PCol.*Pscale;
            obj.P = setPtot(obj);
        end

        function obj = plus(obj1,obj2)
            % Overloading '+' can be used for operations on signal_interface objects.
            param = struct('Fs',obj1.Fs);
            N = obj1.N;
            L = obj1.L;
            if param.Fs~=obj2.Fs || N~=obj2.N || L~=obj2.L
                slog('Sampling rate, and signal sizes of both signals must be equal.','ERR');
            end
            param.Rs = obj1.Rs;
            if param.Rs~=obj2.Rs
                slog('Assuming symbol rate of the first signal.', 'WRN');
            end
            param.PCol = obj1.PCol+obj2.PCol;   %TEMPORARY - for SNR only; power should be calculated numerically b/c of coherence issues
            
            %frequency offset
            Fc1 = obj1.Fc;
            Fc2 = obj2.Fc;
            param.Fc = (Fc1+Fc2)/2;
            s1 = bsxfun(@times,getScaled(obj1),exp(2j*pi*(Fc1-param.Fc)/param.Fs*(0:L-1)'));
            s2 = bsxfun(@times,getScaled(obj2),exp(2j*pi*(Fc2-param.Fc)/param.Fs*(0:L-1)'));
            
            % waveform addition
            s = s1+s2;
            
            %power tracking
            avpower = pwr.meanpwr(s);
            for jj = 1:N
                param.PCol(jj) = pwr(param.PCol(jj).SNR, {avpower(jj), 'W'});
            end
            
            obj = signal_interface(s,param);
            obj.P = setPtot(obj);
        end

        function s = getScaled(obj)
            % get power scaling waveform
            s = obj.sig;
            Pin = pwr.meanpwr(s);
            for jj=1:obj.N
                Pout(jj) = obj.PCol(jj).Ptot('W');
            end
            if Pin ~=0
                s = bsxfun(@times, s, sqrt(Pout./Pin));
            else
                if Pout>0
                    slog('Cannot scale waveform with 0 power to %.2e W','WRN', Pout)
                end
            end
        end

        function obj = set(obj,varargin)
            rescaleP = 0;
            rescalePCol = 0;
            if numel(varargin)==1
                obj.sig = varargin{1}; % TODO input checking --- can only be a matrix
                if numel(obj.PCol)~=size(obj.sig, 2)
                    slog('The number of "pwr" objects in PCol must match the number of columns of the signal.', 'WRN');
                    slog('Resetting quasi-arbitrarily.  This warning will become an error in the future.', 'WRN');
                    if numel(obj.PCol)>size(obj.sig, 2)
                    	obj.PCol = obj.PCol(1:size(obj.sig, 2));
                    else
                        obj.PCol = repmat(obj.P/obj.N, obj.N, 1);
                    end
                    obj.P = setPtot(obj);
                end
            elseif rem(nargin,2)
                for i=1:(nargin-1)/2
                    obj.(varargin{2*i-1}) = varargin{2*i};
                    if strcmp(varargin{2*i-1}, 'P'), rescalePCol=1; end
                    if strcmp(varargin{2*i-1}, 'PCol')
                        if numel(varargin{2*i}) ~= size(obj.sig, 2)
                            slog('The number of "pwr" objects in PCol must match the number of columns of the signal.', 'ERR');
                        end
                        rescaleP=1;
                    end
                end
            else
                slog('Bad key-value pairs.','ERR');
            end
            if rescaleP, obj.P = setPtot(obj); end
            if rescalePCol, obj.PCol = repmat(obj.P/obj.N, obj.N, 1); end
            if rescalePCol && rescaleP
                slog('Both power per column and total power were specified.  Power per column takes precedence')
            end
        end

        function obj = combine(obj1, obj2)
            % combine two signal_interface obj into one signal_interface obj
            slog('To combine signal1 and signal2, assuming that "Fc", "Fs", "Rs" of these two signal are equal.', 'BSCIF');
            objparams = obj1.params;
            if obj1.N ~= obj2.N
                slog('The modes of the two signals must be equal to combine.', 'ERR');
            end
            objparams.L = obj1.L + obj2.L;
            s = [obj1.get; obj2.get];
            objparams.PCol = (obj1.PCol.Ptot('W')*obj1.L + obj2.PCol.Ptot('W')*obj2.L) / objparams.L;
            obj = signal_interface(s, objparams);
        end

        function obj = demuxcomplex(obj1)
            % Separate the real and imaginary parts of the complex signal.
            Sin = obj1.get;
            Sout = zeros(size(Sin, 1), 2*size(Sin, 2));
            Sout(:, 1:2:end) = real(Sin);
            Sout(:, 2:2:end) = imag(Sin);
            avgpwr = pwr.meanpwr(Sout);
            p = params(obj1);
            for idx = 1:2*obj1.N
                p.PCol(idx) = pwr(obj1.PCol(round(idx/2)).SNR, {avgpwr(idx), 'w'});
            end
            obj = signal_interface(Sout, p);
        end

        function obj = cut(obj1, L1, L2)
            % Cut the signal
            s = obj1.get;
            s = s(L1:L2, :);
            obj = signal_interface(s, obj1.params);
        end

        function obj = reserve(obj, idx)
            s = obj.get;
            sout = s(:, idx);
            p = params(obj);
            p = rmfield(p, 'PCol');
            avgpwr = pwr.meanpwr(sout);
            p.P = pwr(obj.PCol(idx).SNR, {avgpwr, 'W'});
            obj = signal_interface(sout, p);
        end

        function obj = combineinto1(varargin)
            slog('To combine signal1 and signal2, assuming that "Fc", "Fs", "Rs" of these two signal are equal.', 'BSCIF');
            obj = varargin{1};
            for idx = 2:length(varargin)
                s = [obj.get, varargin{idx}.get];
                p = obj.params;
                p.PCol = [obj.PCol, varargin{idx}.PCol];
                obj = signal_interface(s, p);
            end
            obj.set('P', setPtot(obj));
        end

    end

end