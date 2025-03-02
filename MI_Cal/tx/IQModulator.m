classdef IQModulator < unit 

    properties
        nInputs = 2;
        nOutputs = 1;
        Vdc;
        Vpi;
        Vbias;
        Phase;
        IQGainImbalance;
        scaledVdirveEnabled;
        extinctionRatio;
        insertionLoss;
    end

    properties(SetAccess=private)
        nModes;
    end

    methods
        
        function obj = IQModulator(varargin)
            if nargin == 1
                param = varargin{1};
                obj.Vpi = paramdefault(param, 'Vpi', 5);
                obj.Vdc = paramdefault(param, 'Vdc', 5);
                obj.Vbias = paramdefault(param, 'Vbias', -1);
                obj.Phase = paramdefault(param, 'Phase', pi/2);
                obj.IQGainImbalance = paramdefault(param, 'IQGainImbalance', 0);
                obj.insertionLoss = paramdefault(param, 'insertionLoss', 6);
                obj.extinctionRatio = paramdefault(param, 'extinctionRatio', 35);
            elseif nargin == 0
                obj.Vpi = 5;
                obj.Vdc = 5;
                obj.Phase = pi/2;
                obj.IQGainImbalance = 0;
                obj.insertionLoss = 6;
                obj.extinctionRatio = 35;
            else
                slog('Too many input arguments', 'ERR');
            end
        end

        function out = generate(obj, drive, laser)
            
            if ~isreal(drive.getRaw)
                drive = demuxcomplex(drive);
            end

            obj.nModes = drive.N / 2;
            inputRs = drive.Rs;

            if mod(obj.nModes, 1)
                slog('If the input drive is a real number sequence, then N must be a multiple of 2.', 'ERR');
            end

            p = params(drive);
            pCol = p.PCol;
            for idx = 1:obj.nModes
                index = 2*(idx-1) + 1;
                pCol(index) = p.PCol(index) * 10^(-obj.IQGainImbalance/20);
                pCol(index+1) = p.PCol(idx+1) * 10^(obj.IQGainImbalance/20);
            end

            drive = set(drive, 'PCol', pCol);

            if laser.N > 1
                slog('The laser can have only one polarization.', 'ERR');
            end

            if laser.L > drive.L
                outL = drive.L;
                slog('Taking output length from drive', 'BSCIF');
            else
                outL = laser.L;
                slog('Taking output length from laser', 'BSCIF');
            end

            drivesig = drive.get;
            lasersig = laser.get;

            out = zeros(outL, obj.nModes);
            for idx = 1:obj.nModes
                index = 2*(idx-1) + 1;
                outI = MZM(drivesig(1:outL, index), lasersig(1:outL, 1)/sqrt(2), obj.Vpi, obj.Vdc, obj.Vbias, obj.extinctionRatio, obj.insertionLoss);
                outQ = MZM(drivesig(1:outL, index+1), exp(1j*obj.Phase)*lasersig(1:outL, 1)/sqrt(2), obj.Vpi, obj.Vdc, obj.Vbias, obj.extinctionRatio, obj.insertionLoss);
                out(:, idx) = outI/sqrt(2) + outQ/sqrt(2);
            end
            clear drive;

            outparams = laser.params;
            outparams = rmfield(outparams, 'PCol');
            outparams.Rs = inputRs;
            outparams.P = pwr(inf, {pwr.meanpwr(out) / obj.nModes, 'W'});
            avgpwr = 10^(outparams.P.Ptot/10) * 1e-3 * pwr.meanpwr(out) / mean(pwr.meanpwr(out));
            for idx = 1:size(out, 2)
                outparams.PCol(idx) = pwr(inf, {avgpwr(idx), 'W'});
            end
            outparams = rmfield(outparams, 'P');
            out = signal_interface(out, outparams);
        end
    end
end