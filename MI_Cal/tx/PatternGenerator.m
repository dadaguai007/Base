% PRBS或rand信号生成器
classdef PatternGenerator < unit

    properties
        nInputs = 0;
        nOutputs = 1;
        prbsOrder;
        M;
        N;
        sigLength;
        seed;
        type;
        zeroPercent;
    end
    
    methods

        function obj = PatternGenerator(param)
            if isfield(param, 'prbsOrder') && isscalar(param.prbsOrder), obj.prbsOrder = param.prbsOrder; else, obj.prbsOrder = 15; end
            if isfield(param, 'M') && isscalar(param.M), obj.M = param.M; else, slog('Please input the modulation order "M".', 'ERR'); end
            if isfield(param, 'N') && isscalar(param.N), obj.N = param.N; else, obj.N = 1; end
            if isfield(param, 'sigLength') && isscalar(param.sigLength), obj.sigLength = param.sigLength; else, obj.sigLength = 1e4; end
            if isfield(param, 'type') && ischar(param.type), obj.type = param.type; else, obj.type = 'prbs'; end
            if isfield(param, 'zeroPercent') && isscalar(param.zeroPercent), obj.zeroPercent = param.zeroPercent; else, obj.zeroPercent = 0.5; end
            chnum = log2(obj.M) * obj.N;

            % 产生种子
            if ~isfield(param, 'seed')
                slog('Seeds will be randomly generated.', 'BSCIF')
                obj.seed = randi([0,1], chnum, obj.prbsOrder);
%                 obj.seed = zeros(chnum, obj.prbsOrder);
            end
            
        end

        function out = generate(obj)
            chnum = log2(obj.M) * obj.N;
            if strcmp(obj.type, 'prbs')
                data = obj.prbs_generator(obj.prbsOrder, obj.seed, obj.sigLength);
            else
                data = obj.random_generator(obj.zeroPercent, chnum, obj.sigLength);
            end
            out = signal_interface(double(data), struct('Rs', 1, 'Fs', 1, 'P', pwr(inf, 0), 'Fc', 0));
        end
    end

    methods(Static)

        function out = prbs_generator(order, seed, sigLength)
            chnum = size(seed, 1);
            G = zeros(1, order); 
            switch(order)
                case 7
                    G(6:7) = 1;
                case 9
                    G(5) = 1;
                    G(9) = 1;
                case 11
                    G(11) = 1;
                    G(9) = 1;
                case 15
                    G(14:15) = 1;
                case 20
                    G(3) = 1;
                    G(20) = 1;
                case 23
                    G(18) = 1;
                    G(23) = 1;
                case 31
                    G(28) = 1;
                    G(31) = 1;
                otherwise
                    G(order-1:order) = 1;
            end
            pow = nextpow2(sigLength+1);
            index = find(G==1);
            leng = length(index);
            out = zeros(sigLength,chnum);
            for idx = 1:chnum
                reg = seed(idx,:);
                if pow <= order
                    for i = 1:sigLength
                        temp = reg(index(1));
                        for j = 2:leng
                            temp = xor(temp,reg(index(j)));
                        end
                        reg = [temp,reg(1:1:order-1)];
                        out(i, idx) = reg(order);
                    end
                else
                    T = zeros(2^order-1,1);
                    T(1) = reg(1);
                    for i = 1:2^order-1
                        temp = reg(index(1));
                        for j = 2:leng
                            temp = xor(temp,reg(index(j)));
                        end
                        reg = [temp,reg(1:1:order-1)];
                        T(i) = reg(order);
                    end
                    times = floor(sigLength/(2^order-1));
                    out(:,idx) = [repmat(T,times,1);T(1:sigLength-(2^order-1)*times)];
                end
            end
        end

        function out = random_generator(zeroPercent, chnum, lengthSequence)
            out = rand(lengthSequence, chnum) > zeroPercent;
        end
    end

end