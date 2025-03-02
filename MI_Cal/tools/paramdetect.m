function judge = paramdetect(param, type, varargin)
    judge = zeros(1, length(varargin));
    switch lower(type)
        case 'scale'
            for idx = 1:length(varargin)
                judge(idx) = isfield(param, varargin{idx}) * isscalar(param.(varargin{idx}));
            end
        case 'char'
            for idx = 1:length(varargin)
                judge(idx) = isfield(param, varargin{idx}) * ischar(param.(varargin{idx}));
            end
        case 'vector'
            for idx = 1:length(varargin)
                judge(idx) = isfield(param, varargin{idx}) * isvector(param.(varargin{idx}));
            end
        otherwise
            slog('This function only supports "scale", "char", "vector" param detection.', 'ERR');
    end
end