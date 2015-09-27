function [index, value] = fastMLCut(x, l, r, varargin)
% [index, value] = fastMLCut(x, l, r, varargin)
% get the cut that maximizes the empirical likelihood via a line search
% 
% x: 1-D array
% l: left boundary
% r: right boundary
% optional: 'plot' (default no plotting)
% 
% index: the index of x that is the rightmost point lying left of cutting point
% value: the value of the cutting point
x = x(:)';
N = length(x);
[y, I] = sort(x); % increasing order
jvec = 1:(N-1); % cuts are numbered 1 ... N-1
assert(y(1)>=l & y(end)<=r, 'x not in range of [l,r]');
if N>500
    index_vec = randsample(N, 500, true)';
    cuts = y(index_vec);
    goodness_of_cuts = index_vec .* log(index_vec ./ (cuts-l)) + (N-index_vec) .* log((N-index_vec) ./ (r-cuts));
    [~, j] = max(goodness_of_cuts);
    index = I(index_vec(j));
    value = x(index);
else
    cuts = [(y(1)+y(2))/2, y(2:end-1)];
    % goodness of cut by directly computing the empirical log liklihood
    goodness_of_cuts = jvec .* log(jvec ./ (cuts-l)) + (N-jvec) .* log((N-jvec) ./ (r-cuts));
    % smooth with slide windowed averaging
%     goodness_of_cuts = conv(goodness_of_cuts, ones(1,5)/5, 'same'); 
    % not allow cutting too small area
%     imbalance = (cuts - l) ./ (r - cuts);
%     J = imbalance<1;
%     imbalance(J) = 1 ./ imbalance(J);
%     goodness_of_cuts(imbalance>10) = -Inf;

    [~, j] = max(goodness_of_cuts);
    index = I(j);
    value = x(index);
end
end