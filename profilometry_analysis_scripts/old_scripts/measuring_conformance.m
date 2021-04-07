%% Measuring conformance
% do to this, you need to visualize the profile and write down the ranges
% for min and max by eye

% need to clean this up.


clear
close all
load('colorscheme.mat')
gel_constant = 1.49;
force = 200; %g

num_widths = 7; %0.35,0.5,0.75,1,1.2,1.5,2 

%nov 15 "35" gel was 2.81% or 34.65
%nov 3 "36" gel was 2.73% or 35.4
%jan 6 "36" gel 1 was 2.72% 35.7
%jan 6 "36" gel 2 was 2.94% or 33
%jan 14 gels 3, 4, and 5 were 2.80% or 34.7
%Jan 29 gels 7 8 9 10 11 12 were 2.7 (first three) 2.65 (second three)

craig_data = load('craig_data_10-Feb-2021');
craig_data = craig_data.craig_data;

% file_names_to_add =  {'210209_craig_gel_7_processed', '210209_craig_gel_8_processed', ...
%     '210209_craig_gel_10_processed', '210209_craig_gel_11_processed', '210209_craig_gel_12_processed'};
% for i = 1:length(file_names_to_add)
%     craig_data.file_names{end+1} = file_names_to_add{i};
% end
% 
% gels_added = length(file_names_to_add);

% craig_data.area_vecs{end+1, 1} = [14,15; 13,14; 11,12; 9,11; 7,9; 4,6; 0.5,4]; %array, 7x2, beginning with 13, 14; i.e., range capturing minima
% craig_data.area_vecs{end, 2} = [15, 15.5; 13.5,14; 12, 12.5; 10, 10.5; 8, 8.5; 6,6.5; 3,4]; %array, 7x2, beginning with 13, 14; i.e., range capturing maxima
% craig_data.area_vecs{end, 3} = 'h';

% craig_data.ratio(end:end+gels_added) = [2.7, 2.7, 2.65, 2.65, 2.65, nan];
% craig_data.days_old(end:end+gels_added) = [11, 11, 11, 11, 11, nan];
% craig_data.gel_id_num(end:end+gels_added) = [7, 8, 10, 11, 12, nan];

num_gels = length(craig_data.file_names);
gels = zeros(num_gels, num_widths);
most_recent = 21; %num most recent gels we want to look at
%% set empirically chosen parameters %CRAIG STIMS 

% each file name is a gel struct of a different gel with the craig
% stimulus.
 %each row is for a craig stim reading 32 is 1, 35 is 2, 36 is 3, 1st entry is min ranges, 
%2nd entry is max ranges. Starts with smallest gap (0.35)

%% run main for loop visualizing cross sections
plot_flag = 0;
for i = num_gels-most_recent+1: num_gels
    cd ../../mat_files
    load(craig_data.file_names{i})
    gel.profile(gel.profile<0) = 0; % by histogram inspection
    gel.profile = gel.profile.*gel_constant;
    cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
    gels(i, :) = find_grating_differences(gel, craig_data.area_vecs{i,3},...
                                               craig_data.area_vecs{i,1}, ...
                                               craig_data.area_vecs{i,2}, ...
                                               plot_flag)';
                                           
end
% 
% c = date;
% % save new craig_data
% save(strcat("craig_data_", c), 'craig_data')

%% plot gel data
plot_flag = 0;
human_x = [0.35, 0.5, 0.75, 1, 1.2, 1.5]'; % in mm, from gibson and craig '06, fig 5
human_y = [0.06, 0.1, 0.17, 0.28, 0.366, 0.4]';
if plot_flag
    figure;
    hold on
    dot_size = 15;
end
x_norm = [0.35, 0.5, 0.75, 1, 1.2, 1.5, 2]'; %stimulus widths
x_36 = [0.5, 0.75, 1, 1.2, 1.5]';
x_gel_2_0115 = [0.5, 0.75, 1, 1.2, 1.5, 2]';
scatters = {};
fits = {};
mean_percent_error = zeros(1,num_gels);
for i = num_gels-most_recent+1:num_gels
    if i == 3 %third one has no first or last entry (im dumb)
        x = x_36;
        y = gels(i, :)';
        y = y(2:end-1); 
        mean_percent_error(i) = meanPercentError(y, human_y(2:end));
    elseif i == 3 %has no first  entry (im dumb)
        x = x_gel_2_0115;
        y = gels(i, :)';
        y = y(2:end); 
        mean_percent_error(i) = meanPercentError(y(1:end-1), human_y(2:end));
    else
        x = x_norm;
        y = gels(i, :)';
        mean_percent_error(i) = meanPercentError(y(1:end-1), human_y);
    end
    %add 0 to x and why because we know that 0 width should have 0
    %protrusion
    
    x = [ x];
    y = [ y];
    fit_ob = fitlm(x,y);
    fits{i} = fit_ob;
    
    if plot_flag
        scatters{i} = scatter(x,y, dot_size, colorscheme(i-num_gels+most_recent+1, :));
        fit_plot = fit(x,y, 'poly1');
        ax{i} = plot(fit_plot);
        ax_handles = ax{i};
        set(ax_handles,'color',colorscheme(i-num_gels+most_recent+1, :));
    end
end

%% add gibson and craig data
x = human_x; % in mm, from gibson and craig '06, fig 5
y = human_y;
i = i+1;
fit_ob = fitlm(x,y);
fits{i} = fit_ob;

if plot_flag
    scatters{i} = scatter(x,y, dot_size, colorscheme(1, :));
    fit_plot = fit(x,y, 'poly1');
    ax{i} = plot(fit_plot);
    ax_handles = ax{i};
    set(ax_handles,'color',colorscheme(1, :));
    xlabel("Grating Width (mm)")
    ylabel("Conformance Depth (mm)")
    title("Width-Conformance relationship at 200g force")
end

%% stats on results - showing change over time for 33
num_gels  = num_gels+1;

rmses = zeros(num_gels,1);
SEs = zeros(num_gels,1);
slopes = zeros(num_gels,1);
intercepts = slopes;
for i = num_gels-most_recent+1:num_gels
    lin_fit = fits{i};
    slopes(i) = lin_fit.Coefficients.Estimate(2);
    SEs(i) = lin_fit.Coefficients.SE(2);
    rmses(i) = lin_fit.RMSE;
    intercepts(i) = lin_fit.Coefficients.Estimate(1);
end


gel_id_targets = [ 7, 8, 10, 11, 12];
% gel_id_targets = [gel_id_targets, 35, 36];
gel_stats = {}; 
% first row is slopes second is standard errors on slopes, third is days
% old reading (x coord)
figure;
hold on
for i = 1:length(gel_id_targets)
gel_stats{1, i} = slopes(craig_data.gel_id_num == gel_id_targets(i));
gel_stats{2, i} = SEs(craig_data.gel_id_num == gel_id_targets(i));
gel_stats{3, i} = craig_data.days_old(craig_data.gel_id_num == gel_id_targets(i));
gel_stats{4, i} = intercepts(craig_data.gel_id_num == gel_id_targets(i));
gel_stats{5, i} = mean_percent_error(craig_data.gel_id_num == gel_id_targets(i));
errorbar(gel_stats{3, i}, gel_stats{1, i}, gel_stats{2, i});
end

yline(slopes(end))
xlabel("days old")
ylabel("slope index")
legend( "7", "8", "10", "11", "12", "human")
title("Slopes")

figure;
hold on
for i = 1:length(gel_id_targets)
    plot(gel_stats{3,i}, gel_stats{4,i})
end
yline(intercepts(end)) %human
legend( "7", "8", "10", "11", "12", "human")
title("intercepts")

figure;
hold on
for i = 1:length(gel_id_targets)
    plot(gel_stats{3,i}, gel_stats{5,i})
end
legend("7", "8", "10", "11", "12")
title("mean percent error")
