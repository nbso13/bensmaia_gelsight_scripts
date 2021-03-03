%% mother of all plots
% Feb 24

close all
clear;

load('activities')

%% Find best amplitude for each texture
num_textures = length(activities.names);
num_amps = length(activities.amps);
mean_percent_errors_gel = zeros(num_amps, num_textures); %percent difference of real value to sim value
mean_percent_errors_ts = zeros(num_amps, num_textures);
for i = 1:num_textures
    real_vals = activities.real(i, 1:3);
    gel_responses = activities.gels{i};
    ts_responses = activities.ts{i};
    for j = 1:num_amps
        gel_percent_diffs = zeros(1,3);
        ts_percent_diffs = zeros(1,3);
        for k = 1:3
            gel_percent_diffs(k) = 100*abs(real_vals(k)- gel_responses(j,k))/real_vals(k);
            ts_percent_diffs(k) = 100*abs(real_vals(k)- ts_responses(j,k))/real_vals(k);
        end
        mean_percent_errors_gel(j,i) = mean(gel_percent_diffs);
        mean_percent_errors_ts(j,i) = mean(ts_percent_diffs);
    end
end

% get the best overall amplitude across g
% get index of minimum of each column
[~, I_gel] = min(mean_percent_errors_gel);
[~, I_ts] = min(mean_percent_errors_ts);
%get associated amplitude
best_amps_gel = activities.amps(I_gel);
best_amps_ts = activities.amps(I_ts);

% grab activities

gel_optimal_activities = zeros(num_textures, 6);
ts_optimal_activities = zeros(num_textures, 6);

for i = 1:num_textures
    gel_activity = activities.gels{i};
    gel_optimal_activities(i, :) = gel_activity(I_gel(i), :);
    ts_activity = activities.ts{i};
    ts_optimal_activities(i, :) = ts_activity(I_ts(i), :);
end


aff_class = ["PCs", "RAs", "SA1s"];
figure;
scatter_list = [];
for i = 1:3
    subplot(2, 3, i)
    hold on;
    title(strcat(aff_class(i), ", Gel"));
    for j = 1:num_textures
        scatter_list(j) = scatter(activities.real(j,i),gel_optimal_activities(j,i), 'filled');
        er = errorbar(activities.real(j,i),gel_optimal_activities(j,i), activities.real(j,i+3), 'horizontal');
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er = errorbar(activities.real(j,i),gel_optimal_activities(j,i), gel_optimal_activities(j,i+3));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    x = 1:max(max(activities.real(:,i)), max(gel_optimal_activities(:, i)));
    y = x;
    plot(x,y);
    xlabel("Real Recorded Activity (Hz)");
    ylabel("Gel Simulated Activity (Hz)");
    legend(scatter_list, activities.names, 'Interpreter', 'none')
end

scatter_list = [];
for i = 1:3
    subplot(2, 3, i+3)
    hold on;
    title(strcat(aff_class(i), ", TouchSim"));
    for j = 1:num_textures
        scatter_list(j) = scatter(activities.real(j,i),ts_optimal_activities(j,i), 'filled');
        er = errorbar(activities.real(j,i),ts_optimal_activities(j,i), activities.real(j,i+3), 'horizontal');
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er = errorbar(activities.real(j,i),ts_optimal_activities(j,i), ts_optimal_activities(j,i+3));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    x = 1:max(max(activities.real(:,i)), max(ts_optimal_activities(:, i)));
    y = x;
    plot(x,y);
    xlabel("Real Recorded Activity (Hz)");
    ylabel("TouchSim Simulated Activity (Hz)");
    legend(scatter_list, activities.names, 'Interpreter', 'none')
end

        
        

