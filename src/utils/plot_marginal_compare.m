function plot_marginal_compare(chains, labels, varargin)
% comparing posterior by plotting marginal densities
% Call
% plot_marginal_compare({chain_1, chain_2, ...}, {name_1, name_2, ...})
% plot_marginal_compare({chain_1, chain_2, ...}, {name_1, name_2, ...}, [p1
% p2]) if only a few dims are compared
% plot_marginal_compare({chain_1, chain_2, ...}, {name_1, name_2, ...}, [..true values of theta])
% if also want to annotate the true values

k = length(chains);
assert(k==length(labels));
[~, p] = size(chains{1});
true_values = [];
if isempty(varargin)
    dims = 1:p;
else
    dims = varargin{1};
    assert(all(dims>=1) && all(dims<=p));
    if length(varargin)>1
        true_values = varargin{2};
        assert(length(true_values)==p);
    end
end

lineStyles = {'-','--',':','-.'};
a = floor(sqrt(length(dims)));
b = ceil(length(dims)/a);
for jj=1:length(dims)
    % plot dim j
    j = dims(jj);
    cnt = 0;
    subplot(a,b,jj);
    hold on;
    for c=1:k
        [f, x] = ksdensity(chains{c}(1:10:end,j));
        if cnt==0
            % plot truth in red
            plot(x, f, '-r', 'DisplayName', labels{c}, 'LineWidth', 2);
        else
            plot(x, f, 'DisplayName', labels{c}, 'LineWidth', 2, 'LineStyle', lineStyles{mod(cnt, length(lineStyles))+1});
        end
        cnt = cnt + 1;
    end
    if ~isempty(true_values)
        % also show the true values
        line([true_values(j), true_values(j)], ylim, 'Color','red', 'LineStyle', '--');
    end

    if j==1, legend('-DynamicLegend', 'Location', 'best'); end
    title(['parameter ',num2str(j)]);
    hold off;
end
end

