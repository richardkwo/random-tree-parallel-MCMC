function [index, value] = midPointCut(x, l, r, varargin)
assert(x(1)>=l & x(end)<=r, 'x not in range of [l,r]');
value = (max(x)+min(x))/2;
[~, index] = min(abs(x-value));
end
