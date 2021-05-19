%% Process Gain
clear
close all
% Try convolution!!!
% Try going pixel by pixel
cd ../../mat_files
load('201204_gain_35_gel_processed');
new_gel = gel;
load('201119_gain_no_gel_processed');
load('201119_gain_gel_processed');
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts


%MAKE COLOR BARS ALL THE SAME SCALE

% 
% gel.profile = gel.profile-0.035;
% no_gel.profile = no_gel.profile-0.04;
% gel.profile(gel.profile<0) = 0;
% no_gel.profile(no_gel.profile<0) = 0;


gel = cropProfile(gel, 'bottom', 3, 'px');
no_gel = cropProfile(no_gel, 'left', 3, 'px');
no_gel = cropProfile(no_gel, 'right', 4, 'px');
no_gel = cropProfile(no_gel, 'left', 0.026, 'mm');
gel = cropProfile(gel, 'right', 0.025, 'mm');


% 
% compare_lines(gel, no_gel, 10.7, 10.7, 'h')
% 
% compare_lines(gel, no_gel, 13.7, 13.7, 'h')


ground_swath = [4,5,12,13];
no_gel_bottom = sectionbyMm(no_gel, ground_swath);
gel_bottom = sectionbyMm(gel, ground_swath);

no_gel_bottom_mean = mean(no_gel_bottom, 'all');
gel_bottom_mean = mean(gel_bottom, 'all');
no_gel.profile = no_gel.profile -(no_gel_bottom_mean-gel_bottom_mean);

new_gel.profile = imrotate(new_gel.profile, 90);
temp = new_gel.x_axis;
new_gel.x_axis = new_gel.y_axis;
new_gel.y_axis = temp;
new_gel = cropProfile(new_gel, 'bottom', 0.32, 'mm');
new_gel = cropProfile(new_gel, 'top', 2.198, 'mm');
new_gel = cropProfile(new_gel, 'left', 0.407, 'mm');
new_gel = cropProfile(new_gel, 'right', 0.217, 'mm');
new_no_gel = cropProfile(no_gel, 'left', 1.514, 'mm');
new_no_gel = cropProfile(new_no_gel, 'top', 3.488, 'mm');



new_no_gel_bottom = sectionbyMm(new_no_gel, [5,7,6,8]);
new_gel_bottom = sectionbyMm(new_gel, [5,7,6,8]);
new_no_gel_bottom_mean = mean(new_no_gel_bottom, 'all');
new_gel_bottom_mean = mean(new_gel_bottom, 'all');
new_no_gel.profile = new_no_gel.profile -(new_no_gel_bottom_mean-new_gel_bottom_mean);



if ~checkSizeMatch(new_gel, new_no_gel)
    [new_gel, new_no_gel] = resampleToMin(new_gel, new_no_gel); %resamples to the min resolution
    [new_gel, new_no_gel] = bruteCropFit(new_gel, new_no_gel); %crops to same size
end

% compare_lines(new_gel, new_no_gel, 7.5, 7.45, 'h')
% compare_lines(new_gel, new_no_gel, 1.3, 1.3, 'h')
% compare_lines(gel, no_gel, 7.7, 7.7, 'v')

location_mat = vertcat( [6, 6.2, 8.5, 9], [6.1, 6.2, 1.3, 1.4], ...
    [0.6, 0.8, 7.3, 7.6], [3.3, 3.5, 4.2, 4.4], [0.6, 0.7, 4.3, 4.5]);
bottom_mat = vertcat( [3.5, 4, 9, 9.5], [7, 7.4, 1, 2], ...
    [2.5, 3, 7, 7.5], [5, 5.5, 4, 4.5], [1.5, 2, 4, 4.5]);
gel_vals = zeros(1,size(location_mat, 1));
no_gel_vals = zeros(1, size(location_mat, 1));

figure;
subplot(2,2,1)
title("Gel A")
visualizeProfile(new_gel);
max_no_gel = max(new_no_gel.profile(:));
color_vec = ["r", "b", "k", "c", "g"];
subplot(2,2,2)
title("Raw")
visualizeProfile(new_no_gel);
for i = 1:size(location_mat, 1)
    bottom_gel_vals(i) = mean(sectionbyMm(new_gel, bottom_mat(i,:)), 'all');
    bottom_no_gel_vals(i) = mean(sectionbyMm(new_no_gel, bottom_mat(i,:)), 'all');
    gel_vals(i) = mean(sectionbyMm(new_gel, location_mat(i,:)), 'all') - bottom_gel_vals(i);
    no_gel_vals(i) = mean(sectionbyMm(new_no_gel, location_mat(i,:)), 'all') - bottom_no_gel_vals(i);
    
    loc_vec = location_mat(i,:);
    new_loc = [loc_vec(1), loc_vec(3), loc_vec(2)-loc_vec(1), loc_vec(4)-loc_vec(3)];
    loc_vec = bottom_mat(i,:);
    new_bottom_loc = [loc_vec(1), loc_vec(3), loc_vec(2)-loc_vec(1), loc_vec(4)-loc_vec(3)];
    
    
    subplot(2,2,1); hold on; rectangle('Position', new_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    rectangle('Position', new_bottom_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    subplot(2,2,2); hold on; rectangle('Position', new_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    rectangle('Position', new_bottom_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    
    subplot(2,2,3)
    hold on
    scatter(no_gel_vals(i), gel_vals(i), color_vec(i), 'filled');
end

x = no_gel_vals;
y = gel_vals;

subplot(2,2,3)
hold on

fit_ob = fit(x',y', 'poly1');
m = fit_ob.p1;
n = fit_ob.p2;
ax = plot(fit_ob);
legend(ax, sprintf('y=%f*x + %f',m, n))
xlabel("Dot Height No Gel (mm)")
ylabel("Dot Height Gel (mm)")
title("Measuring Gel Gain")

subplot(2,2,4);
diff = new_no_gel.profile - (1/m).*new_gel.profile;
new_gel.profile = diff;
min_diff = min(diff(:));
visualizeProfile(new_gel);
title("Scaled Difference")
sgtitle("Gel A Gain Analysis")

subplot(2,2,1); caxis([0 max_no_gel]);
subplot(2,2,2); caxis([0 max_no_gel]);
subplot(2,2,4); caxis([min_diff max_no_gel+min_diff]);
% for j = [1] %%trying to see what would minimize variance
%     x = 0.5:0.02:2;
%     y = zeros(size(x));
%     z = y;
%     count = 1;
%     for i = x
%         if j == 1
%         diff = gel;
%         diff.profile = no_gel.profile - gel.profile./i;
%         else
%             diff = new_gel;
%             diff.profile = new_no_gel.profile - new_gel.profile./i;
%         end
%         y(count) = mean(diff.profile, 'all');
%         z(count) = std(diff.profile(:));
%         count = count+1;
%         
%     end
%     figure
%     hold on
%     scatter(x,y, 5, 'b')
%     scatter(x, z, 5, 'r')
% end

%imshowpair(no_gel.profile, gel.profile,'Scaling', 'joint')

% visualizeProfile(gel)
% visualizeProfile(no_gel)
% predicted_ratio = 0.66666 ;
% 
% resolution = 0.008;
% gainAnalysis(gel, no_gel, resolution, predicted_ratio)

%% comparing peaks
location_mat = vertcat( [7.6, 7.9, 12.2, 12.5], [2.1, 2.4, 10.7, 10.9], ...
    [4.9, 5.2, 7.6, 7.9], [7.7, 7.8, 4.7, 4.8], [2.2, 2.3, 13.7, 13.8]);
bottom_mat = vertcat( [8.4, 8.8, 8.5, 9], [4, 4.5, 11, 11.5], ...
    [6.5, 7, 6, 6.5], [8.5, 9, 5, 5.5], [3, 3.5, 12.5, 13]);
gel_vals = zeros(1,size(location_mat, 1));
no_gel_vals = zeros(1, size(location_mat, 1));

figure;
subplot(2,2,1)
title("Gel B")
visualizeProfile(gel);
max_no_gel = max(no_gel.profile(:));
color_vec = ["r", "b", "k", "c", "g"];
subplot(2,2,2)
title("Raw")
visualizeProfile(no_gel);
for i = 1:size(location_mat, 1)
    bottom_gel_vals(i) = mean(sectionbyMm(gel, bottom_mat(i,:)), 'all');
    bottom_no_gel_vals(i) = mean(sectionbyMm(no_gel, bottom_mat(i,:)), 'all');
    gel_vals(i) = mean(sectionbyMm(gel, location_mat(i,:)), 'all') - bottom_gel_vals(i);
    no_gel_vals(i) = mean(sectionbyMm(no_gel, location_mat(i,:)), 'all') - bottom_no_gel_vals(i);
    
    loc_vec = location_mat(i,:);
    new_loc = [loc_vec(1), loc_vec(3), loc_vec(2)-loc_vec(1), loc_vec(4)-loc_vec(3)];
    loc_vec = bottom_mat(i,:);
    new_bottom_loc = [loc_vec(1), loc_vec(3), loc_vec(2)-loc_vec(1), loc_vec(4)-loc_vec(3)];
    
    
    subplot(2,2,1); hold on; rectangle('Position', new_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    rectangle('Position', new_bottom_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    subplot(2,2,2); hold on; rectangle('Position', new_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    rectangle('Position', new_bottom_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    
    subplot(2,2,3)
    hold on
    scatter(no_gel_vals(i), gel_vals(i), color_vec(i), 'filled');
end

x = no_gel_vals;
y = gel_vals;

subplot(2,2,3)
hold on

fit_ob = fit(x',y', 'poly1');
m = fit_ob.p1;
n = fit_ob.p2;
ax = plot(fit_ob);
legend(ax, sprintf('y=%f*x + %f',m, n))
xlabel("Dot Height No Gel (mm)")
ylabel("Dot Height Gel (mm)")
title("Measuring Gel Gain")

subplot(2,2,4);
diff = no_gel.profile - (1/m).*gel.profile;
gel.profile = diff;
min_diff = min(diff(:));
visualizeProfile(gel);
title("Scaled Difference")
sgtitle("Gel B Gain Analysis")

subplot(2,2,1); caxis([0 max_no_gel]);
subplot(2,2,2); caxis([0 max_no_gel]);
subplot(2,2,4); caxis([min_diff max_no_gel+min_diff]);

