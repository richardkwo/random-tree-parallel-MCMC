function plot_pairwise_compare(chains, labels, varargin)
% comparing posterior by pairwise plotting samples
% Call
% plot_pairwise_compare({chain_1, chain_2, ...}, {name_1, name_2, ...})
% plot_pairwise_compare({chain_1, chain_2, ...}, {name_1, name_2, ...}, [p1
% p2]) if only a few dims are compared

k = length(chains);
assert(k==length(labels));
[~, p] = size(chains{1});
if isempty(varargin)
    dims = 1:p;
else
    dims = varargin{1};
    assert(all(dims>=1) && all(dims<=p));
end

a = length(dims);
for u=1:length(dims)
    i = dims(u);
    for v=1:length(dims);
        j = dims(v);
        h = subplot(a,a,length(dims)*(u-1)+v);
        if i==j
            text(0.5, 0.5, ['parameter ',num2str(j)], 'Parent', h);
        else
            hold on;
            for c=k:-1:1
                if size(chains{c},1)>1000
                    x = randsample(size(chains{c},1), 1000);
                end
                if c==1
                    plot(chains{c}(x,i), chains{c}(x,j), 'ko', 'MarkerSize', 2, 'DisplayName', labels{c});
                else 
                    plot(chains{c}(x,i), chains{c}(x,j), 'o', 'MarkerSize', 2, 'DisplayName', labels{c})
                end
            end
            hold off;
            if i==1 && j==2
                legend('-DynamicLegend', 'Location', 'best'); 
            end
            xlabel(['parameter ',num2str(i)]);
            ylabel(['parameter ',num2str(j)]);
        end
    end
end
end

