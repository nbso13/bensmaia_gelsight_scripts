function [] = errorbars_group(bar_data, error_data, significance)
    % Find the number of groups and the number of bars in each group
    [ngroups, nbars] = size(bar_data);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    
    % for text
    x = (1:ngroups) - groupwidth/2 + (2*1.5-1) * groupwidth / (2*nbars); %only works for 2 bar groups
    y = max(bar_data, [], 2);
    for i = 1:length(x)
        if significance(i) < 0.05
            text(x(i), y(i)+ y(i)/6, "*", 'FontSize',14, 'FontWeight', 'bold')
        end
    end
    
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, bar_data(:,i), error_data(:,i), 'k', 'linestyle', 'none', 'HandleVisibility','off');
    end
end
