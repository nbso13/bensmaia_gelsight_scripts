%% Measuring conformance
% do to this, you need to visualize the profile and write down the ranges
% for min and max by eye
clear
close all
load('colorscheme.mat')
gel_constant = 1.56;
force = 200; %g
num_gels = 5; %1:32, 1:35_201116, 1:36 , 1:35_201118 .... 1:36old
num_widths = 7; %0.35,0.5,0.75,1,1.2,1.5,2 
gels = zeros(num_gels, num_widths);
plot_gels = [0, 1, 1, 1, 1]; %determines which gels are visualized (see order below in filenames)

%% set empirically chosen parameters %CRAIG STIMS 
file_names = {'201111_craig_stim_32_gel_processed', ...
    '201116_craig_stim_35_gel_processed', '201111_craig_stim_36_gel_processed', '201118_craig_35_gel_processed', '201204_craig_36_gel_processed'};
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
%% run main for loop
plot_flag = 1;
for i = 1: num_gels
    cd mat_files\
    load(file_names{i})
    cd ..
    gel.profile(gel.profile<0) = 0; % by histogram inspection
    gel.profile = gel.profile.*gel_constant;
    gels(i, :) = find_grating_differences(gel, area_vecs{i,3},...
                                               area_vecs{i,1}, ...
                                               area_vecs{i,2}, ...
                                               plot_flag)';
end

%% plot gel data
figure;
hold on
dot_size = 15;
color_str = ['b', 'g', 'r', 'c', 'k'];
x_norm = [0.35, 0.5, 0.75, 1, 1.2, 1.5, 2]';
x_36 = [0.5, 0.75, 1, 1.2, 1.5]';
m = zeros(num_gels+1); %slopes
n = m; %intercepts
ax = {}; %trendlines
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
    end
    scatter(x,y, dot_size, colorscheme(i, :));
    fit_ob = fit(x,y, 'poly1');
    m(i) = fit_ob.p1;
    n(i) = fit_ob.p2;
    ax{i} = plot(fit_ob);
    ax_handles = ax{i};
    set(ax_handles,'color',colorscheme(i, :));
end

%% add gibson and craig data
x = [0.35, 0.5, 0.75, 1, 1.2, 1.5]'; % in mm, from gibson and craig '06, fig 5
y = [0.06, 0.1, 0.17, 0.28, 0.366, 0.4]';
i = i+1;
scatter(x,y, dot_size, colorscheme(i+1, :));
fit_ob = fit(x,y, 'poly1');
m(i) = fit_ob.p1;
n(i) = fit_ob.p2;
ax{i} = plot(fit_ob);

%'1:32 gel', '1:32 gel trendline', ...
legend('1:35 gel', '1:35 gel trendline', '1:36 gel', '1:36 gel trendline', ...
    '1:35_201118 gel', '1:35_201118 trendline', 'Gibson & Craig, 2006 (human)', ...
    'G&B trendline', '1:36 gel dec', '1:36 dec trendline');
xlabel("Grating Width (mm)")
ylabel("Conformance Depth (mm)")
title("Width-Conformance relationship at 200g force")


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
