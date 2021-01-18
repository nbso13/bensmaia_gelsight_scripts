%% Measuring conformance
% do to this, you need to visualize the profile and write down the ranges
% for min and max by eye
clear
close all
load('colorscheme.mat')
gel_constant = 1.49;
force = 200; %g

num_widths = 7; %0.35,0.5,0.75,1,1.2,1.5,2 

plot_gels = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; %determines which gels are visualized (see order below in filenames)
num_gels = length(plot_gels);
gels = zeros(num_gels, num_widths);
%% set empirically chosen parameters %CRAIG STIMS 
%nov 15 "35" gel was 2.89% or 34.65
%nov 3 "36" gel was 2.82% or 35.4
%jan 6 "36" gel 1 was 2.80% 35.7
%jan 6 "36" gel 2 was 3.03% or 33

days_old = [nan, 1, 8, 3, 30, 2, 2, 5, 5, 6, 6, nan];
ratio = [nan, 34.65, 35.4, 34.65, 35.4, 35.7, 33, 35.7, 33, 35.7, 33, nan];
file_names = {'201111_craig_stim_32_gel_processed', ...
    '201116_craig_stim_35_gel_processed', '201111_craig_stim_36_gel_processed', ...
    '201118_craig_35_gel_processed', '201204_craig_36_gel_processed', ...
    '210108_craig_36_gel_1_processed', '210108_craig_36_gel_2_processed', ...
    '210111_craig_36_gel_1_processed', '210111_craig_36_gel_2_processed' ...
    '210112_craig_36_gel_1_processed', '210112_craig_36_gel_2_processed'};
% each file name is a gel struct of a different gel with the craig
% stimulus.
area_vecs = {}; %each row is for a craig stim reading 32 is 1, 35 is 2, 36 is 3, 1st entry is min ranges, 
%2nd entry is max ranges. Starts with smallest gap (0.35)
area_vecs{1,1} = [13, 14; 12, 13; 10, 11.5; 8, 10; 6, 8; 3, 6; 0.1, 3];
area_vecs{1,2} = [14, 16; 12.5, 13.5; 11, 12; 9, 10; 7, 8; 5, 6; 2, 4];
area_vecs{1,3} = 'v';
area_vecs{2,1} = [0.01, 1; 1, 3; 3, 4; 4, 6; 6, 8.5; 8.5, 11; 11, 13.5];
area_vecs{2,2} = [1, 1.4; 2, 3; 4, 5; 6, 7; 8, 10; 10, 12; 12, 13.8];
area_vecs{2,3} = 'h';
area_vecs{3,1} = [0.1, 0.2; 10, 11; 8, 10; 6, 8; 3, 5; 0.01, 3; 0.1, 0.2]; %ALERT: 0.35 and 2mm gaps not measured here. NAN in data.
area_vecs{3,2} = [0.1, 0.2; 11, 11.5; 9, 10; 7, 8; 5, 6; 2, 4; 0.1, 0.2];       
area_vecs{3,3} = 'h';
area_vecs{4,1} = [13, 15; 12, 13; 10, 12; 8, 10; 6, 8; 3, 6; 1, 3]; 
area_vecs{4,2} = [14, 15; 12.5, 13; 11, 12; 9, 10; 7, 8; 5, 6; 2, 3.5];       
area_vecs{4,3} = 'h';
area_vecs{5,1} = [13, 14; 11, 13; 10, 11;  7, 9; 5, 7; 2, 5; 0.5, 2 ]; 
area_vecs{5,2} = [14, 14.5; 12, 13.5; 11, 12; 9, 10; 7, 8;  5, 6;  1, 3.5 ];       
area_vecs{5,3} = 'h';
area_vecs{6,1} = [13, 14; 11, 13; 10, 11;  7, 9; 5, 7; 2, 5; 0.5, 2 ]; 
area_vecs{6,2} = [14, 14.5; 12, 13.5; 11, 12; 9, 10; 7, 8;  5, 6;  1, 3.5 ];       
area_vecs{6,3} = 'h';
area_vecs{7,1} = [13, 14; 11, 13; 10, 11;  7, 10; 6, 8; 4, 6; 0.5, 2 ]; 
area_vecs{7,2} = [13.9, 14.5; 12.5, 13.5; 11, 11.3; 9, 10; 7, 8;  5, 6;  1, 3.5 ];       
area_vecs{7,3} = 'h';
area_vecs{8,1} = [14, 14.2; 12, 13; 10, 12;  8, 10; 6, 8; 4, 6; 0.5, 4 ]; 
area_vecs{8,2} = [14.25, 14.5; 12.9, 13.5; 11, 12; 9, 10; 7, 8;  5, 6;  2, 3.5 ];       
area_vecs{8,3} = 'h';
area_vecs{9,1} = [13.96, 14.5; 12.5, 13; 10.7, 11.5;  9, 9.4; 6.7, 7.3; 4, 5; 0.5, 3 ]; 
area_vecs{9,2} = [14.25, 14.45; 12.9, 13.1; 11, 12; 9.6, 9.9; 7.6, 7.9;  5, 5.7;  2.8, 3.2 ];       
area_vecs{9,3} = 'h';
area_vecs{10,1} = [13.7, 14; 12, 13; 10, 12;  8, 10; 6, 8; 4, 5; 0.5, 3 ]; 
area_vecs{10,2} = [14, 14.2; 12.7, 13; 11, 11.5; 9, 10; 7, 8;  5, 6;  2, 3.5 ];       
area_vecs{10,3} = 'h';
area_vecs{11,1} = [13.96, 14.5; 12.5, 13; 10.7, 11.5;  9, 9.4; 6.7, 7.3; 4, 5; 0.5, 3 ]; 
area_vecs{11,2} = [14.25, 14.45; 12.9, 13.1; 11, 12; 9.6, 9.9; 7.6, 7.9;  5, 5.7;  2.8, 3.2 ];       
area_vecs{11,3} = 'h';
%% run main for loop
plot_flag = 1;
for i = 1: num_gels
    cd ../../mat_files
    load(file_names{i})
    gel.profile(gel.profile<0) = 0; % by histogram inspection
    gel.profile = gel.profile.*gel_constant;
    cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
    gels(i, :) = find_grating_differences(gel, area_vecs{i,3},...
                                               area_vecs{i,1}, ...
                                               area_vecs{i,2}, ...
                                               plot_flag)';
end

%% plot gel data
figure;
hold on
dot_size = 15;
x_norm = [0.35, 0.5, 0.75, 1, 1.2, 1.5, 2]'; %stimulus widths
x_36 = [0.5, 0.75, 1, 1.2, 1.5]';
m = zeros(num_gels+1, 1); %slopes
n = m; %intercepts
ax = {}; %trendlines
scatters = {};
fits = {};
for i = 1:num_gels
    
    if ~plot_gels(i) %if we're not to print a gel, skip
        continue
    end
    if i == 3 %third one has no first or last entry (im dumb)
        x = x_36;
        y = gels(i, :)';
        y = y(2:end-1); 
    else
        x = x_norm;
        y = gels(i, :)';
        x = [0; 0; x];
        y= [0; 0; y];
    end
    %add 0 to x and why because we know that 0 width should have 0
    %protrusion
    
    scatters{i} = scatter(x,y, dot_size, colorscheme(i, :));
    fit_ob = fitlm(x,y);
    fits{i} = fit_ob;
    fit_plot = fit(x,y, 'poly1');
    ax{i} = plot(fit_plot);
    ax_handles = ax{i};
    set(ax_handles,'color',colorscheme(i, :));
end

%% add gibson and craig data
x = [0.35, 0.5, 0.75, 1, 1.2, 1.5]'; % in mm, from gibson and craig '06, fig 5
y = [0.06, 0.1, 0.17, 0.28, 0.366, 0.4]';
i = i+1;
scatters{i} = scatter(x,y, dot_size, colorscheme(i, :));
fit_ob = fitlm(x,y);
fits{i} = fit_ob;
fit_plot = fit(x,y, 'poly1');
ax{i} = plot(fit_plot);
ax_handles = ax{i};
set(ax_handles,'color',colorscheme(i, :));

% %automating plotting legend
% plot_legends = zeros(length(plot_gels)*2, 1)';
% for i = 1:length(plot_gels)
%     plot_legends(i*2) = plot_gels(i);
%     plot_legends(i*2 - 1) = plot_gels(i);
% end
% legend_labels = {'1:32 gel', '1:32 gel trendline', '1:35 gel', '1:35 gel trendline', '1:36 gel', '1:36 gel trendline', ...
%     '1:35_201118 gel', '1:35_201118 trendline', '1:36 gel dec', '1:36 dec trendline'};
% legend_labels = legend_labels(plot_legends);
% legend_labels = [legend_labels, 'Gibson & Craig, 2006 (human)', 'G&B trendline'];
% 
% 
% 
% legend(legend_labels);
xlabel("Grating Width (mm)")
ylabel("Conformance Depth (mm)")
title("Width-Conformance relationship at 200g force")
legend( "2.72%, 6 days old", "trend", "2.94%, 6 days old", "trend", "finger", "trend");
%% stats on results - showing change over time for 33
rmses = zeros(num_gels+1,1);
SEs = zeros(num_gels+1,1);
slopes = zeros(num_gels+1,1);
for i = 1:num_gels+1
    if i == 3
        continue
    end
    lin_fit = fits{i};
    slopes(i) = lin_fit.Coefficients.Estimate(2);
    SEs(i) = lin_fit.Coefficients.SE(2);
    rmses(i) = lin_fit.RMSE;
end
x = [2, 5, 6];
slope33 = slopes(ratio == 33);
SE33 = SEs(ratio == 33);
slope35 = slopes(ratio == 35.7);
SE35 = SEs(ratio == 35.7);
figure;
hold on
errorbar(x, slope33, SE33)
errorbar(x, slope35, SE33)
yline(slopes(end))
x = [1, 3];
slope34 = slopes(ratio == 34.65);
SE34 = SEs(ratio == 34.65);
errorbar(x, slope34, SE34);
% 
% 
% %% OLD MEASURE CONF
% clear
% close all
% gel_constant = 1.56;
% force = 200; %g
% num_gels = 2; %1:32, 1:36 
% num_widths = 2; %1mm, 2mm
% gels = zeros(num_gels, num_widths);
% 
% %% set empirically chosen parameters 
% file_names = {"Gel_grating_201021", ...
%     "Gel_grating_201019", 'craig_stim_32_gel_processed'; ...
%     '201104_small_grating_gel_processed', ...
%     '201104_large_grating_gel_processed', 'craig_stim_36_gel_processed'};
% %each row of cells contains gratings for one gel. for each grating, there
% %is a file.
% area_vecs = {[0, 8, 6.52, 6.66; 0, 8, 6.035, 6.065],...
%     [0, 9, 6.4, 6.65; 0, 9, 5.52, 5.57];...
%     [0, 5, 1.5, 1.6; 0, 5, 1.05, 1.07],...
%     [0, 5, 4.9, 5.1; 0, 5, 3.98, 4.02]};
% % each row of cells is for a different gel. Each cell in that row is a
% % different grating. The first line of the array in the cell is the floor and the second
% % is the top for vectors describing the window in (x1, x2, y1, y2) for the
% % corresponding profilometry file.
% 
% %% run main for loop
% for i = 1: num_gels
%     for j = 1:num_widths
%         cd mat_files\
%         load(file_names{i,j})
%         cd ..
%         gel.profile(gel.profile<0) = 0; % by histogram inspection
%         gel.profile = gel.profile.*gel_constant;
%         gel_vecs = area_vecs{i,j};
%         gel_floor = gel_vecs(1,:);
%         gel_top = gel_vecs(2,:);
%         take_input = 0;
%         gels(i,j) = peakAndValley(take_input, gel, gel_floor, gel_top);
%     end
% end
% 
% %% plot gel data
% figure;
% hold on
% color_str = ['b', 'r', 'k'];
% width_2mm = 2; %mm
% width_1mm = 1;
% x = [width_1mm, width_2mm]';
% m = zeros(num_gels+1); %slopes
% n = m; %intercepts
% ax = {}; %trendlines
% for i = 1:num_gels
%     y = gels(i, :)'; 
%     scatter(x,y, color_str(i));
%     fit_ob = fit(x,y, 'poly1');
%     m(i) = fit_ob.p1;
%     n(i) = fit_ob.p2;
%     ax{i} = plot(fit_ob, color_str(i));
% end
% 
% %% add gibson and craig data
% x = [0.35, 0.5, 0.75, 1, 1.2, 1.5]'; % in mm, from gibson and craig '06, fig 5
% y = [0.06, 0.1, 0.17, 0.28, 0.366, 0.4]';
% i = i+1;
% scatter(x,y, color_str(i));
% fit_ob = fit(x,y, 'poly1');
% m(i) = fit_ob.p1;
% n(i) = fit_ob.p2;
% ax{i} = plot(fit_ob, color_str(i));
% 
% legend('1:32 gel', '1:32 gel trendline', '1:36 gel', '1:36 gel trendline', 'Gibson & Craig, 2006 (human)', 'G&B trendline');
% xlabel("Grating Width (mm)")
% ylabel("Conformance Depth (mm)")
% title("Width-Conformance relationship at 200g force")
% 
% 
