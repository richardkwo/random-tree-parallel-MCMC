function [index, value] = kdCut(x, l, r, varargin)
assert(x(1)>=l & x(end)<=r, 'x not in range of [l,r]');
index = ceil(length(x)/2);
value = median(x);
end
