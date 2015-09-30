function option = part_options(varargin)
% option = part_options('param_1', 'value_1', ..., 'param_k', 'value_k')
% configure an option for running PART algorithm
%  
% Parameters:
% min_cut_length: \delta_a, minimum side length of a block. (default: 0.001)
% min_fraction_block: \delta_\rho, minimum fraction of samples contained by a block (default: 0.01)
% ntree: number of trees (default: 16)
% resample_N: number of resamples drawn from an aggregated pair of sub-chains (default: 10000)
% cut_type: partition rule, must be 'kd' or 'ml' (default: 'kd')
% parallel: whether trees are built in parallel with MATLAB's parallel computing toolbox, true or false (default: true)
% verbose: level of verbose (default: 1)
% local_gaussian_smoothing: whether using local Gaussian smoothing applied to each block, true or false (default: true)
% match: better matching sub-chains into pairs, true or false (default: true)
% halving: whether halving min_fraction_block after each stage (default: false)
% resample_data_only: whether only resampling original data points (default: false)
% 
% Example:
% option = part_options('min_cut_length', 1e-2, 'ntree', 40)

option.min_cut_length = 0.001;
option.min_fraction_block = 0.01;
option.ntree = 16;
option.resample_N = 10000;
option.cut_type = 'kd';
option.verbose = 1;
option.parallel = true;
option.local_gaussian_smoothing = true;
option.match = true;
option.halving = false;
option.resample_data_only = false;


% read the acceptable names
optionNames = fieldnames(option);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
  error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
  inpName = pair{1}; 
  if any(strcmp(inpName, optionNames))
      %# overwrite options. If you want you can test for the right class here
      %# Also, if you find out that there is an option you keep getting wrong,
      %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      option.(inpName) = pair{2};
      switch inpName
      	case 'min_cut_length'
      		assert(pair{2}>0, 'minimum length of block must be positive');
      	case 'min_fraction_block'
      		assert(pair{2}>0 && pair{2}<1, 'minimum fraction of samples contained in a block must be >0 and <1');
      	case 'cut_type'
      		assert(strcmp(pair{2}, 'kd') || strcmp(pair{2}, 'ml'), 'cut_type must be kd or ml');
      	case 'resample_N'
      		assert(pair{2}>=1000, 'resample_N too small');
      	case 'resample_data_only'
      		warning('resample_data_only = true only for debugging');
      	case 'ntree'
      		assert(pair{2}>0, 'ntree must be a positive integer');
      end
  else
      error('%s is not a recognized parameter name',inpName)
  end
end

fprintf('Options configured:\n');
disp(option);

end