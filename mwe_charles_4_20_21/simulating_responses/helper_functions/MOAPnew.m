function [rmses, rs, mean_ratios] = MOAPnew(activities)
%motherOfAllPlots new takes in activity statistics from real afferents and
%simulated afferents from touchsim and gelsight and compares the simulation
%to the real values. difference is that now we plot on three plots not six.

num_textures = length(activities.names);
rmses = zeros(3, 2);
rmse_errors = rmses;
rs = zeros(3, 2);
RUPs = rs;
RLOs = rs;
mean_ratios = zeros(3, 2);
sd_ratios = mean_ratios;
all_ratios = {};

aff_class = ["PCs", "RAs", "SA1s"];
figure;
% scatter_list = [];
for i = 1:3 %for each afferent
    subplot(1, 3, i)
    hold on;
    x = 1:max([max(activities.real(:,i)), ...
        max(activities.gel(:, i)), max(activities.ts(:, i))]);
    y = x;
    plot(x,y, 'k', 'LineStyle','--');
    title(strcat(aff_class(i)));
    for j = 1:num_textures
        scatter(activities.real(j,i),activities.gel(j,i), [], 'r', 'o', 'filled');
        scatter(activities.real(j,i), activities.ts(j,i), [], 'b', '+');
    end
     if (i == 1)
        strs = {'unity line', 'gel profiles', 'touchsim profiles' };
        colors = [ 0,0,0; 1, 0, 0; 0, 0, 1];
        leg = [color_legend(strs', colors)];
    	leg = legend(leg);
        leg.Box = 0;
        ylabel("Simulated Activity (Hz)");
     end
    pbaspect([1 1 1])
    
    all_ratios{i, 1} = activities.gel(:,i)./activities.real(:,i);
    mean_ratios(i, 1) = mean(all_ratios{i, 1});
    sd_ratios(i, 1) = std(all_ratios{i, 1})/sqrt(length(all_ratios{i, 1}));
    [r, p, RLO, RUP] = corrcoef(activities.real(:,i),activities.gel(:,i));
    rs(i, 1) = r(1,2);
    RUPs(i,1) = RUP(1,2)- r(1,2);
    RLOs(i,1) = -RLO(1,2)+r(1,2);
    rmses(i, 1) = rmse_calc(activities.real(:,i),activities.gel(:,i));
    xlabel("Real Recorded Activity (Hz)");
    ax = gca;
    ax.FontSize = 12;
    ax.FontWeight = 'bold';
end

for i = 1:3
    all_ratios{i, 2} = activities.ts(:,i)./activities.real(:,i);
    mean_ratios(i, 2) = mean(all_ratios{i, 2});
    sd_ratios(i, 2) = std(all_ratios{i, 2})/sqrt(length(all_ratios{i, 2}));
    [r, p, RLO, RUP] = corrcoef(activities.real(:,i),activities.ts(:,i));
    rs(i, 2) = r(1,2);
    RUPs(i,2) = RUP(1,2)-r(1,2);
    RLOs(i,2) = -RLO(1,2)+r(1,2);
    rmses(i, 2) = rmse_calc(activities.real(:,i),activities.ts(:,i));
end

%% bar plot and p_value MEAN RATIO and STANDARD ERROR MEAN
figure;

aff_names = ["PCs", "RAs", "SAs"];
significance = zeros(1,3);
for i = 1:3
    [~, sig] = ttest(all_ratios{i, 1}, all_ratios{i, 2});
    significance(i) = sig;
end
b = bar(mean_ratios);
b(1).FaceColor = [0, 0.45, 0.74]; b(2).FaceColor = [1, 0, 0];
title("Mean Ratio of Simulated to Recorded Firing Rates Across Textures")
xticklabels(aff_names);
ylabel("Mean Ratio (+/- SEM)")
% ylim([-1 1]);
%errorbars
hold on
errorbars_group(mean_ratios, sd_ratios, significance);
yline(1, 'HandleVisibility','off');
strs = {'Gel', 'Raw'}';
colors = [[0, 0.45, 0.74];[1, 0, 0]];
leg = color_legend(strs, colors);
leg = legend(leg);
leg.Box = 0;
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';



%% bar plot and p_value CORRELATION COEFF and 95 CI
figure;

aff_names = ["PCs", "RAs", "SAs"];
b = bar(rs);
b(1).FaceColor = [0, 0.45, 0.74]; b(2).FaceColor = [1, 0, 0];
title("Correlation Coefficient Between Simulated and Recorded Firing Rates Across Textures")
xticklabels(aff_names);
ylabel("Correlation (+/- 95% CI)")
% ylim([-1 1]);
%errorbars
hold on
errorbars_group_both(rs, RLOs, RUPs, [1, 1, 1]);
strs = {'Gel', 'Raw'}';
colors = [[0, 0.45, 0.74];[1, 0, 0]];
leg = color_legend(strs, colors);
leg = legend(leg);
leg.Box = 0;
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';

%% bar plot and p_value RMSE
figure;

aff_names = ["PCs", "RAs", "SAs"];
b = bar(rmses);
b(1).FaceColor = [0, 0.45, 0.74]; b(2).FaceColor = [1, 0, 0];
title("RMSEs between Simulated and Recorded Firing Rates Across Textures")
xticklabels(aff_names);
ylabel("RMSE")
% ylim([-1 1]);
%errorbars
hold on
errorbars_group(rmses, rmse_errors, [1, 1, 1]);
strs = {'Gel', 'Raw'}';
colors = [[0, 0.45, 0.74];[1, 0, 0]];
leg = color_legend(strs, colors);
leg = legend(leg);
leg.Box = 0;
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';


end

function [rmse] = rmse_calc(real_rates, sim_rates)
    rmse = sqrt(mean((real_rates-sim_rates).^2));
end

function [] = errorbars_group_both(bar_data, error_data_down, error_data_up,  significance)
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
        errorbar(x, bar_data(:,i), error_data_down(:,i), error_data_up(:,i), 'k', 'linestyle', 'none', 'HandleVisibility','off');
    end
end
