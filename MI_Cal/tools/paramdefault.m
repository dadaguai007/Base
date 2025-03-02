% This is a function that helps to default the input parameters
% 

function out = paramdefault(param, ref, value)

	judge = isfield(param, ref);
	if ~iscell(ref)
		ref = {ref};
	end
	ref = ref(judge);
	if any(judge)
		out = param.(ref{end});
	else
		out = value;
	end

end