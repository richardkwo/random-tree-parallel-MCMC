function output_chains = thinning(chains, N, varargin)
% Call: 
% output_chain = thinning(chain, N)
% {chain_1, chain_2,..} = thinning({chain_1, chain_2,..}, N)
% chains = thinning(chains, N, 'burn', burn_in)
% 
% thin a chain (or a cell of chains) to the desired length N, 
% optionally after discarding a burn-in period at head

burn_in = 0;
if ~isempty(varargin) && strcmp(varargin{1},'burn')
    burn_in = varargin{2};
end

if isa(chains, 'double')
    chains = chains(burn_in+1:end, :);
    if size(chains,1)>N
        idx = floor(linspace(1,size(chains,1),N));
        output_chains = chains(idx, :);
    else
        output_chains = chains;
    end
else
    assert(isa(chains, 'cell'), 'wrong chain');
    output_chains = cell(size(chains));
    for m=1:numel(chains)
        chains{m} = chains{m}(burn_in+1:end, :);
        if size(chains{m},1)>N
            idx = floor(linspace(1,size(chains{m},1),N));
            output_chains{m} = chains{m}(idx, :);
        else
            output_chains{m} = chains{m};
        end
    end
end

end

