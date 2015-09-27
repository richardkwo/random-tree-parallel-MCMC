function [trees, sampler, sampler_prob] = MultiStageMCMC(MCdraws, option)
% MultiStageMCMC is a function that takes an input of raw subset posteriors
% and combine the posterior densities via a multistage manner, i.e., a
% pairwise combining until a final one density is reached. The function
% calls OneStageMCMC for combining each pair of posterior samples
%
% Calls:
% [trees, sampler, sampler_prob] = MultiStageMCMC(MCdraws, option)
%
% Arguments:
%
% MCdraws: a 1 x m cell containing subset posteriors from m subsets
% option: option configured by part_options(...). Use `option = part_options()` for default settings. 
%
% Outputs:
%
% trees: A list containing multiple partition trees, each corresponding to a density estimation of the combined posterior
% sampler: A list containing multiple flattened partition trees
% sampler_prob: A list containing the probabilities for each flattened partition tree
%
% See also:
% part_options, buildForest, OneStageMCMC
   
    
    %Setting the default value for parameters
    if ~isfield(option,'area'), area = []; else area = option.area;end
    min_fraction_block = option.min_fraction_block;
    if option.halving
        fprintf('\nPairwise aggregation with min fraction = %f, halving enabled. \n', min_fraction_block);
    else
        fprintf('\nPairwise aggregation with min fraction = %f, halving disabled. \n', min_fraction_block);
    end
    
    %Initialization
    m = length(MCdraws); %number of subsets
    d = size(MCdraws{1}, 2); %number of dimensions
    N = zeros(1, m); % number of posterior samples of each subsets
    
    tic;
    %To avoid duplicate computing, we give the default area from the very top
    if isempty(area)
        area = zeros(d, 2);
        area(:,1) = min(MCdraws{1});
        area(:,2) = max(MCdraws{1});
        for i = 1:m
            N(i) = size(MCdraws{i},1); %update posterior sample sizes
            area(:,1) = min([area(:,1)';MCdraws{i}]); %update boundaries
            area(:,2) = max([area(:,2)';MCdraws{i}]);
        end
        
        %enlarge the area by a factor of 1.01
        for j = 1:d
            if area(j,1)>0
                area(j,1) = area(j,1)/1.01;
            else
                area(j,1) = area(j,1)*1.01;
            end
            if area(j,2)>0
                area(j,2) = area(j,2)*1.01;
            else
                area(j,2) = area(j,2)/1.01;
            end
        end
    else
        for i = 1:m
            N(i) = size(MCdraws{i},1);
        end
    end
    
    %inner_option is used by OneStageMCMC for combining paired subsets. We
    %inhibit the information display from there
    inner_option = option;
    inner_option.verbose = 0;
    number_of_stages = floor(log(m)/log(2));
    if option.halving
        current_min_fraction_block = min_fraction_block * (2^number_of_stages);
    else
        current_min_fraction_block = min_fraction_block;
    end
    inner_option.min_number_points_block = max(ceil(min(N) * current_min_fraction_block), 3);
    
    %Now we start combining posterior densities in a pairwise manner
    
    %We use multiple tree at each stage

    MC_left = m; %Number of remaining subsets
    Internal_list = MCdraws; %Remaining subset posteriors.

    %setting up waitbar
    if option.verbose>1
        h = waitbar(0,'Start building multistage tree ensemble...');
    end

    %Initialize stage number (for waitbar use)
    stage = 1;

    %Starting pairwise combining!
    while MC_left > 3
        %update waitbar information
        if option.verbose>1
            waitbar(0.01, h, ['Building multistage tree ensemble for stage ', num2str(stage), ' ...']);
        end

        MC_match = cell(1,floor(MC_left/2)); %Restoring all pairs

        %Initialize pairs as 1-2, 3-4, 5-6, ...
        for MC = 1:floor(MC_left/2)
            MC_match{MC} = [2*MC-1 2*MC];
        end

        %If the number of subset is odd, we add the last subset to the
        %previous pair. So like 1-2, 3-4, 5-6-7...
        if MC_left > floor(MC_left/2)*2;
            MC_match{floor(MC_left/2)} = [MC_match{floor(MC_left/2)} MC_left];
        end

        %If we want to do a better match
        if option.match
            %updating waitbar
            if option.verbose>1
                waitbar(0.01,h,'Starting matching subsets...');
            end

            %Compute means for all subset posteriors
            MC_means = zeros(d,MC_left);
            for MC = 1:MC_left
                MC_means(:,MC) = mean(Internal_list{MC})';
            end

            %match the subset with another subset in median distance
            MC_used = []; %store subset that has been used
            MC_left_set = 1:MC_left; %initialize subsets haven't been paired up
            pair = 1; %current number of pairs + 1
            for pair_first = 1:MC_left
                %update waitbar
                if option.verbose>1
                    waitbar(1-length(MC_used)/MC_left,h,'Matching subsets...');
                end

                if ~ismember(pair_first, MC_used) %check if the subset has been used
                    %compute distances
                    dis = MC_means(:,MC_left_set) -...
                        repmat(MC_means(:,pair_first),1,length(MC_left_set));
                    dis = sum(dis.^2,1);

                    %sort distances
                    [~,ix] = sort(dis);

                    %The subset to pair is the one possesses the median
                    pair_second = MC_left_set(ix(ceil(length(MC_left_set)/2)));

                    %argmenting the used subsets and remaining subsets
                    MC_used = [MC_used, pair_first, pair_second];
                    MC_left_set = setdiff(1:MC_left, MC_used);

                    %saving the pair into MC_match
                    MC_match{pair} = [pair_first pair_second];

                    %Terminating for odd or even number of subsets
                    if pair == floor(MC_left/2) && MC_left > 2*pair
                        last = setdiff(1:MC_left, MC_used);
                        MC_match{pair} = [MC_match{pair}, last];
                        break;
                    elseif pair == floor(MC_left/2)
                        break;
                    end
                    pair = pair+1;
                end
            end
        end

        MC_left = floor(MC_left/2); %Number of subsets left after the current round
        list_temp = cell(1, MC_left); %Storing the combined posterior samples for each pair

        %Starting merging
        for MC_pair = 1:MC_left
            %setting up waitbar
            if option.verbose>1
                waitbar(MC_pair/MC_left,h,['Pairwise combining subsets for stage ', num2str(stage), '...']);
            end

            %Call OneStageMCMC to combine the pair of subsets
            [~, sampler, prob, RawMCMC] = OneStageMCMC(Internal_list(MC_match{MC_pair}), inner_option);
            if option.resample_data_only
                inner_option.mark = RawMCMC.mark;
                inner_option.list = RawMCMC.list;
            else
                inner_option.list = [];
                inner_option.mark = [];
            end

            %Resample from the combined density
            inner_N = max([max(N), option.resample_N]);
            list_temp{MC_pair} = treeSampling(sampler, prob, inner_N, inner_option);
        end
        
        Internal_list = list_temp; %update the posterior samples for the remaining subsets
        stage = stage + 1; %moving forward!
        if option.halving
            % if "halve", min block size is halved after each stage
            current_min_fraction_block = current_min_fraction_block / 2;
        end
        inner_option.min_number_points_block = max(ceil(current_min_fraction_block * inner_N), 3);
    end

    %close waitbar
    if option.verbose>1
        close(h);
    end

    %For final stage, we call OneStageMCMC the last time to combine the
    %final two or three subsets.
    option.min_number_points_block = inner_option.min_number_points_block;
    option.min_cut_length = inner_option.min_cut_length;
    [trees, sampler, sampler_prob] = OneStageMCMC(Internal_list,option);
    
    fprintf('Pairwise aggregation: finsihed in %f seconds \n', toc);
    
end