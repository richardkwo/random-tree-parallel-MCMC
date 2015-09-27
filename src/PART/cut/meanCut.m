function [index, value] = meanCut(x, l, r, varargin)
assert(x(1)>=l & x(end)<=r, 'x not in range of [l,r]');
value = mean(x);
[~, index] = min(abs(x-value));
end
