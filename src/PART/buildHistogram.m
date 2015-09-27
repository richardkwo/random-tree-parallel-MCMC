function cuts = buildHistogram(x, l, r, cutFun, varargin)
% generic 1-D histogram building with a supplied cutting function
% cuts = buildHistogram(x, l, r, cutFun)
% 
% x: data
% l: left boundary
% r: right boundary
% cutFun: [index value] = cutFun(x, l, r) 
% optional:
% 'minNodeSize': 10 (default)
% 'minDepth': 3 (default)
% 
% cuts: a set of cutting points

% cuts = histogramMaxLCut(x, l, r, varargin)
% 
assert(all(x>=l) && all(x<=r), 'not all x in the range [l,r]');
options = struct('minNodeSize', 10, 'minDepth', 3);

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = pair{1}; %# make case insensitive
   if any(strcmp(inpName,optionNames))
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

cuts = [l r];
depth = 0; 
while true
    new_cuts = [];
    for k=1:length(cuts)-1
        x_sub = x(x>=cuts(k) & x<cuts(k+1));
        if length(x_sub)<=options.minNodeSize || abs(cuts(k+1)-cuts(k))<0.01
            continue
        end
        [~, cut] = cutFun(x_sub, cuts(k), cuts(k+1));
%         fprintf('[%f, %f] -> %f\n', cuts(k), cuts(k+1), cut);
        n_left = sum(x_sub<=cut);
        n_right = sum(x_sub>cut);
        if n_left>options.minNodeSize && n_right>options.minNodeSize
            % accept this good cut
            new_cuts = [new_cuts, cut];
        end
    end
    
    if isempty(new_cuts)
        if depth < options.minDepth
            warning('Depth = %d < minimum depth = %d', depth, options.minDepth);
        end
        break
    else
        depth = depth + 1;
        cuts = sort([cuts, new_cuts]);
    end 
end
if length(cuts)>2
    cuts = cuts(2:end-1);
else
    cuts = [];
end

end