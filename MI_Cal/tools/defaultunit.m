% file defaultunit.m
% This function judges if the unit is correct or not.

function varargout = defaultunit(defaults, in)

inleng = length(in);

if ~iscell(defaults)
	defaults = {defaults};
end

varargout = defaults;
varargout(1:inleng) = in(1:inleng);