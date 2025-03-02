function init
persistent flag;
if isempty(flag), flag = false; end

if ~flag
    addpath(genpath(fileparts(mfilename('fullpath'))));
    flag = 1;
end