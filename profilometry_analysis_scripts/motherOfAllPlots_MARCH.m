%% mother of all plots
% Feb 24

close all
clear;

load('activities_max_indentation_final')

%% Find best amplitude for each texture
num_textures = length(activities.names);

aff_class = ["PCs", "RAs", "SA1s"];
figure;
scatter_list = [];
for i = 1:3 %for each afferent
    subplot(2, 3, i)
    hold on;
    title(strcat(aff_class(i), ", Gel"));
    for j = 1:num_textures
        scatter_list(j) = scatter(activities.real(j,i),activities.gel(j,i), 'filled');
        er = errorbar(activities.real(j,i), activities.gel(j,i), activities.real(j,i+3), 'horizontal'); %real
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er = errorbar(activities.real(j,i),activities.gel(j,i), activities.gel(j,i+3));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    x = 1:max(max(activities.real(:,i)), max(activities.gel(:, i)));
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
        scatter_list(j) = scatter(activities.real(j,i), activities.ts(j,i), 'filled');
        er = errorbar(activities.real(j,i),activities.ts(j,i), activities.real(j,i+3), 'horizontal');
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er = errorbar(activities.real(j,i),activities.ts(j,i), activities.ts(j,i+3));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    x = 1:max(max(activities.real(:,i)), max(activities.ts(:, i)));
    y = x;
    plot(x,y);
    xlabel("Real Recorded Activity (Hz)");
    ylabel("TouchSim Simulated Activity (Hz)");
    legend(scatter_list, activities.names, 'Interpreter', 'none')
end

        
        

