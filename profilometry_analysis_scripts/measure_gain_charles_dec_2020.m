%% Process Gain
clear
close all
% Try convolution!!!
% Try going pixel by pixel
cd ..
cd mat_files
load('201204_gain_35_gel_processed');
new_gel = gel;
load('201119_gain_no_gel_processed');
load('201119_gain_gel_processed');
cd ..
cd bensmaia_gelsight_scripts

% 
% gel.profile = gel.profile-0.035;
% no_gel.profile = no_gel.profile-0.04;
% gel.profile(gel.profile<0) = 0;
% no_gel.profile(no_gel.profile<0) = 0;


gel = cropProfile(gel, 'bottom', 3, 'px');
no_gel = cropProfile(no_gel, 'left', 3, 'px');
no_gel = cropProfile(no_gel, 'right', 4, 'px');

compare_lines(gel, no_gel, 10.7, 10.7, 'h')

compare_lines(gel, no_gel, 13.7, 13.7, 'h')


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
new_gel = cropProfile(new_gel, 'bottom', 0.25, 'mm');
new_gel = cropProfile(new_gel, 'top', 2.198, 'mm');
new_gel = cropProfile(new_gel, 'left', 0.407, 'mm');
new_gel = cropProfile(new_gel, 'right', 0.217, 'mm');
new_no_gel = cropProfile(no_gel, 'left', 1.514, 'mm');
new_no_gel = cropProfile(new_no_gel, 'top', 3.488, 'mm');

% new_no_gel_bottom = sectionbyMm(new_no_gel, [5,7,6,8]);
% new_gel_bottom = sectionbyMm(new_gel, [5,7,6,8]);
% new_no_gel_bottom_mean = mean(new_no_gel_bottom, 'all');
% new_gel_bottom_mean = mean(new_gel_bottom, 'all');
% new_no_gel.profile = new_no_gel.profile -(new_no_gel_bottom_mean-new_gel_bottom_mean);



if ~checkSizeMatch(new_gel, new_no_gel)
    [new_gel, new_no_gel] = resampleToMin(new_gel, new_no_gel); %resamples to the min resolution
    [new_gel, new_no_gel] = bruteCropFit(new_gel, new_no_gel); %crops to same size
end

compare_lines(new_gel, new_no_gel, 7.5, 7.45, 'h')
compare_lines(new_gel, new_no_gel, 1.3, 1.3, 'h')
compare_lines(gel, no_gel, 7.7, 7.7, 'v')
y = [0.75, 0.25, 0.25, 0.5, 0.5 1];
x = [0.52, 0.175, 0.17, 0.35, 0.345, 0.67];

% y = [y, 0, 0, 0, 0, 0, 0];
% x = [x, 0, 0, 0, 0, 0, 0];

figure;
hold on
scatter(x, y);

fit_ob = fit(x',y', 'poly1');
m = fit_ob.p1;
n = fit_ob.p2;
ax = plot(fit_ob);
legend(ax, sprintf('y=%f*x + %f',m, n))
xlabel("Dot Height Gel (mm)")
ylabel("Dot Height No Gel (mm)")
title("Measuring Gel Gain")


for j = [1] %%trying to see what would minimize variance
    x = 0.5:0.02:2;
    y = zeros(size(x));
    z = y;
    count = 1;
    for i = x
        if j == 1
        diff = gel;
        diff.profile = no_gel.profile - gel.profile./i;
        else
            diff = new_gel;
            diff.profile = new_no_gel.profile - new_gel.profile./i;
        end
        y(count) = mean(diff.profile, 'all');
        z(count) = std(diff.profile(:));
        count = count+1;
        
    end
    figure
    hold on
    scatter(x,y, 5, 'b')
    scatter(x, z, 5, 'r')
end

%imshowpair(no_gel.profile, gel.profile,'Scaling', 'joint')

% visualizeProfile(gel)
% visualizeProfile(no_gel)
predicted_ratio = 0.66666 ;

resolution = 0.008;
gainAnalysis(gel, no_gel, resolution, predicted_ratio)

%% comparing peaks

gel_cross_height = mean(sectionbyMm(gel, [7.6, 7.9, 12.2, 12.5]), 'all');
no_gel_cross_height = mean(sectionbyMm(no_gel, [7.6, 7.9, 12.2, 12.5]), 'all');
gel_tall_height = mean(sectionbyMm(gel, [2.1, 2.4, 10.7, 10.9]), 'all');
no_gel_tall_height = mean(sectionbyMm(no_gel, [2.1, 2.4, 10.7, 10.9]), 'all');
gel_med_height = mean(sectionbyMm(gel, [4.9, 5.2, 7.6, 7.9]), 'all');
no_gel_med_height = mean(sectionbyMm(no_gel, [4.8, 5.1, 1.5, 1.9]), 'all');
gel_small_height = mean(sectionbyMm(gel, [7.7, 7.8, 4.7, 4.8]), 'all');
no_gel_small_height = mean(sectionbyMm(no_gel, [7.7, 7.8, 4.7, 4.8]), 'all');
gel_small_height_2 = mean(sectionbyMm(gel, [2.2, 2.3, 13.7, 13.8]), 'all');
no_gel_small_height_2 = mean(sectionbyMm(no_gel, [2.2, 2.3, 13.7, 13.8]), 'all');

gel_vals = [gel_cross_height, gel_tall_height, gel_med_height, gel_small_height, gel_small_height_2];
no_gel_vals = [no_gel_cross_height, no_gel_tall_height, no_gel_med_height, no_gel_small_height, no_gel_small_height_2];
y = gel_vals - gel_bottom_mean;
x = no_gel_vals - no_gel_bottom_mean;

y = [y, 0, 0, 0, 0, 0, 0];
x = [x, 0, 0, 0, 0, 0, 0];

figure;
hold on
scatter(x, y);

fit_ob = fit(x',y', 'poly1');
m = fit_ob.p1;
n = fit_ob.p2;
ax = plot(fit_ob);
legend(ax, sprintf('y=%f*x + %f',m, n))
xlabel("Dot Height No Gel (mm)")
ylabel("Dot Height Gel (mm)")
title("Measuring Gel Gain")
%% Plot relationship
% 
% x = no_gello(:);
% y = gello(:);
% 
% figure;
% hold on
% scatter(x,y);
% 
% fit_ob = fit(x,y, 'poly1');
% m = fit_ob.p1;
% n = fit_ob.p2;
% ax = plot(fit_ob);
% legend(ax, sprintf('y=%f*x',m))
% xlabel("Dot Height No Gel (mm)")
% ylabel("Dot Height Gel (mm)")
% title("Measuring Gel Gain")
% 


%% Now for new gel


visualizeProfile(new_no_gel)
visualizeProfile(new_gel)

%imshowpair(new_no_gel.profile, new_gel.profile,'Scaling', 'joint')


resolution = 0.008;
gainAnalysis(new_gel, new_no_gel, resolution, predicted_ratio)

plotLine(new_gel, 9.5, 'h')
plotLine(new_no_gel, 9.5, 'h')