function ind = minEuclid(symb, const)
% Find minimum Euclidean distance.
%
% Find closest constellation symbol w.r.t the Euclidean distance in the
% complex plane.
%
% Parameters
% ----------
% symb : np.array
%     Received constellation symbols.
% const : np.array
%     Reference constellation.
%
% Returns
% -------
% array of int
%     indexes of the closest constellation symbols.

ind = zeros(size(symb));
for ii = 1:length(symb)
    ind(ii) = find(abs(symb(ii) - const) == min(abs(symb(ii) - const)), 1); % find the first index of the minimum distance
end
end