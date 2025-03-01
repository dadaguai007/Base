function [Eo,sampDelay] = decimate(Ei, param)
    % Decimate signal.

    % Input:
    %   Ei: Input signal.
    %   param: - param.SpS_in : samples per symbol of the input signal.
    %    - param.SpS_out: samples per symbol of the output signal.
    %
    % Output:
    %   Eo: Decimated signal.

    if isfield(param, 'SpS_in')
        SpS_in = param.SpS_in;
    end
    if isfield(param, 'SpS_out')
        SpS_out = param.SpS_out;
    end


    [M,N]=size(Ei);
    if M < N
        error('the Ei should be the N × models')
    end
    % Ei should be N×(1,2)
    decFactor = floor(SpS_in / SpS_out);

    % Simple timing recovery
    sampDelay = zeros(1, size(Ei, 2));
    Nmodels =size(Ei, 2);
    % Finds best sampling instant (maximum variance sampling time)
    for k = 1:Nmodels
        a = reshape(Ei(:, k), [], 1);
        varVector = var(reshape(a,SpS_in,[]).');
        [~, idx] = max(varVector);
%         idx = mod(idx-1,decFactor)+1;
        sampDelay(k) = idx;
    end
    
    if sampDelay(k)==SpS_in
       sampDelay(k)=1;
    end
    % Downsampling
    Eo = zeros(size(Ei(1:decFactor:end, :)));
    for k = 1:Nmodels
        Ei(:, k) = circshift(Ei(:, k), -sampDelay(k));
        Eo(:, k) = Ei(1:decFactor:end, k);
    end
%     for k = 1:Nmodels
%         Eo(:, k) = Ei(idx:decFactor:end, k);
%     end


end
