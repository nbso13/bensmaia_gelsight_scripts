function bee_swarm(group_data, x, group_color, max_points)
    ks_underlay = false;
    group_mean = nanmean(group_data, 'all');
   
    if length(group_data) > max_points
        shuffle_idx = randperm(length(group_data));
        group_data = group_data(shuffle_idx(1:max_points));
    end
   
    [B,I] = sort(group_data(:));
    group_data = group_data(I);
   
    n_bins = ceil(1 + log2(numel(group_data)));
    group_hist = hist(group_data(:),n_bins);
    s = 0.4;
    group_hist_prop = (group_hist / max(group_hist)) * (s*1);
    r = max(group_hist_prop);
   
    if ks_underlay
        [freq, freq_x] = ksdensity(group_data(:), 'Bandwidth', 0.8); %
        freq = rescale(freq) * s;
        freq_idx = freq > 0.1;
        freq = freq(freq_idx); freq_x = freq_x(freq_idx);
        x_vec = [x-freq, fliplr(x+freq), x-freq(1)];
        y_vec = [freq_x, fliplr(freq_x), freq_x(1)];
        %patch(x_vec, y_vec, group_color)
        plot(x_vec, y_vec, 'Color', group_color)
    end
    
    for i = 1:length(group_hist)
        if group_hist(i) == 0
            n_bins = n_bins - 1;
            group_hist = [group_hist(1:i-1), group_hist(i+1:end)];
        end
        if length(group_hist) == i
            break
        end
    end
   
    dp_ind = 1;
    for bin = 1:n_bins
        x_range = x + linspace(-group_hist_prop(bin), group_hist_prop(bin), group_hist(bin));
        x_vals = x_range(randperm(length(x_range)));
        y_vals = group_data(dp_ind:dp_ind+group_hist(bin)-1);
        scatter(x_vals, y_vals, 10, 'MarkerFaceColor', group_color, 'MarkerEdgeColor', group_color,'MarkerFaceAlpha', 0.5);
        dp_ind = dp_ind + group_hist(bin);
    end
%     plot([x-r, x+r], [group_mean, group_mean], 'Color' , [0.2 0.2 0.2], 'LineWidth', 2)

end