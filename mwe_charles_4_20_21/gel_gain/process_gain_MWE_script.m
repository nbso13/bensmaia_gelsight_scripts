%% Process Gain

clear
close all
% 
path_to_data_str = '../../../mwe_data/gain_data';
path_back = '../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/gel_gain';

cd(path_to_data_str);
% load("210423_gel_7_gain_stim_200_grams_processed");
load("201119_gain_no_gel_processed")
cd(path_back);

addpath("helper_functions")
color_vec = ["r", "b", "k", "c", "g"];


% 
% % gel 7: align profiles
% 
% 
% 
% no_gel = cropProfile(no_gel, 'top', 1.2, 'mm');
% no_gel = cropProfile(no_gel, 'right', 0.6, 'mm');
% no_gel = cropProfile(no_gel, 'left', 0.7, 'mm');
% % no_gel = cropProfile(no_gel, 'left', 0.026, 'mm');
% 
% gel = cropProfile(gel, 'bottom', 0.025, 'mm');
% 
% %subtract out estimated floor
% ground_swath = [4,5,12,13];
% no_gel_bottom = sectionbyMm(no_gel, ground_swath);
% gel_bottom = sectionbyMm(gel, ground_swath);
% 
% no_gel_bottom_mean = mean(no_gel_bottom, 'all');
% gel_bottom_mean = mean(gel_bottom, 'all');
% no_gel.profile = no_gel.profile -(no_gel_bottom_mean-gel_bottom_mean);
% 
% if ~checkSizeMatch(gel, no_gel)
%     [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
%     [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
% end
% 
% 
% 
% location_mat = vertcat( [7, 7.2, 10.5, 11], [1.5, 1.8, 9.4, 9.8], ...
%     [1.5, 1.6, 12.6, 12.7], [4.1, 4.5, 6.4, 6.8], [4.05, 4.2, 0.3, 0.8]);
% bottom_mat = vertcat( [4, 4.5, 10.5, 11], [3.5, 4, 9, 9.8], ...
%    [3, 3.5, 12.6, 12.7], [6, 6.5, 6, 6.5], [2, 2.5, 1, 2]);
% [scale_factor_7, r_gel7] = gain_process(gel, no_gel, location_mat, bottom_mat, color_vec);

%%  gel 11: align profiles
% 
% cd data
% load("210420_gain_stim_gel_11_200_grams_processed");
% load("201119_gain_no_gel_processed")
% cd ..
% 
% gel.name = "Gel 11 Gain Stim";
% no_gel = cropProfile(no_gel, 'top', 0.95, 'mm');
% no_gel = cropProfile(no_gel, 'left', 1.2, 'mm');
% no_gel = cropProfile(no_gel, 'bottom', 0.2, 'mm');
% 
% % gel = cropProfile(gel, 'bottom', 0.025, 'mm');
% 
% 
% %subtract out estimated floor
% ground_swath = [3,4,12,13];
% no_gel_bottom = sectionbyMm(no_gel, ground_swath);
% gel_bottom = sectionbyMm(gel, ground_swath);
% 
% no_gel_bottom_mean = mean(no_gel_bottom, 'all');
% gel_bottom_mean = mean(gel_bottom, 'all');
% no_gel.profile = no_gel.profile -(no_gel_bottom_mean-gel_bottom_mean);
% 
% if ~checkSizeMatch(gel, no_gel)
%     [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
%     [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
% end
% 
% 
% location_mat = vertcat( [7, 7.2, 10.5, 11], [1.5, 1.8, 9.4, 9.8], ...
%     [1.5, 1.6, 12.6, 12.7], [4.1, 4.5, 6.4, 6.8], [4.15, 4.25, 0.3, 0.8]);
% bottom_mat = vertcat( [4, 4.5, 10.5, 11], [3.5, 4, 9, 9.8], ...
%    [3, 3.5, 12.6, 12.7], [6, 6.5, 6, 6.5], [2, 2.5, 1, 2]);
% location_mat(:, 1:2) = location_mat(:, 1:2)-0.5;
% location_mat(:, 3:4) = location_mat(:, 3:4)+0.1;
% bottom_mat(:, 1:2) = bottom_mat(:, 1:2)-0.5;
% bottom_mat(:, 3:4) = bottom_mat(:, 3:4)+0.1;
% [scale_factor_11, r_gel11] = gain_process(gel, no_gel, location_mat, bottom_mat, color_vec);
% 
% 
% %% gel 18: align profiles
% 
% cd data
% load("210415_gainstim_gel_18_200_grams_processed");
% load("201119_gain_no_gel_processed")
% cd ..
% 
% no_gel = cropProfile(no_gel, 'bottom', 0.1, 'mm');
% no_gel = cropProfile(no_gel, 'left', 1.1, 'mm');
% no_gel = cropProfile(no_gel, 'top', 1, 'mm');
% no_gel = cropProfile(no_gel, 'right', 0.2, 'mm');
% 
% %subtract out estimated floor
% ground_swath = [4,5,12,13];
% no_gel_bottom = sectionbyMm(no_gel, ground_swath);
% gel_bottom = sectionbyMm(gel, ground_swath);
% 
% no_gel_bottom_mean = mean(no_gel_bottom, 'all');
% gel_bottom_mean = mean(gel_bottom, 'all');
% no_gel.profile = no_gel.profile -(no_gel_bottom_mean-gel_bottom_mean);
% 
% 
% if ~checkSizeMatch(gel, no_gel)
%     [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
%     [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
% end
% 
% 
% location_mat = vertcat( [7, 7.2, 10.5, 11], [1.5, 1.8, 9.4, 9.8], ...
%     [1.5, 1.6, 12.6, 12.7], [4.1, 4.5, 6.4, 6.8], [4.15, 4.25, 0.3, 0.8]);
% bottom_mat = vertcat( [4, 4.5, 10.5, 11], [3.5, 4, 9, 9.8], ...
%    [3, 3.5, 12.6, 12.7], [6, 6.5, 6, 6.5], [2, 2.5, 1, 2]);
% location_mat(:, 1:2) = location_mat(:, 1:2)-0.5;
% location_mat(:, 3:4) = location_mat(:, 3:4)+0.1;
% bottom_mat(:, 1:2) = bottom_mat(:, 1:2)-0.5;
% bottom_mat(:, 3:4) = bottom_mat(:, 3:4)+0.1;
% [scale_factor_18, r_gel18] = gain_process(gel, no_gel, location_mat, bottom_mat, color_vec);
% 

% 
% %% gel 18: align profiles
% 
% cd data
% load("210428_gain_stim_gel_19_200_grams_processed");
% load("201119_gain_no_gel_processed")
% cd ..
% 
% no_gel = cropProfile(no_gel, 'bottom', 0.1, 'mm');
% no_gel = cropProfile(no_gel, 'left', 1.1, 'mm');
% no_gel = cropProfile(no_gel, 'top', 1, 'mm');
% no_gel = cropProfile(no_gel, 'right', 0.2, 'mm');
% gel = cropProfile(gel, 'top', 0.9, 'mm');
% 
% %subtract out estimated floor
% ground_swath = [4,5,12,13];
% no_gel_bottom = sectionbyMm(no_gel, ground_swath);
% gel_bottom = sectionbyMm(gel, ground_swath);
% 
% no_gel_bottom_mean = mean(no_gel_bottom, 'all');
% gel_bottom_mean = mean(gel_bottom, 'all');
% no_gel.profile = no_gel.profile -(no_gel_bottom_mean-gel_bottom_mean);
% 
% gel.name = "Gel 19 Gain Stim";
% 
% % if ~checkSizeMatch(gel, no_gel)
% %     [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
% %     [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
% % end
% 
% location_mat = vertcat( [7, 7.2, 10.5, 11], [1.5, 1.8, 9.4, 9.8], ...
%     [1.5, 1.6, 12.6, 12.7], [4.1, 4.5, 6.4, 6.8], [4.15, 4.25, 0.3, 0.8]);
% bottom_mat = vertcat( [4, 4.5, 10.5, 11], [3.5, 4, 9, 9.8], ...
%    [3, 3.5, 12.6, 12.7], [6, 6.5, 6, 6.5], [2, 2.5, 1, 2]);
% location_mat(:, 1:2) = location_mat(:, 1:2)-0.5;
% location_mat(:, 3:4) = location_mat(:, 3:4)+0.1;
% bottom_mat(:, 1:2) = bottom_mat(:, 1:2)-0.5;
% bottom_mat(:, 3:4) = bottom_mat(:, 3:4)+0.1;
% [scale_factor_19, r_gel19] = gain_process(gel, no_gel, location_mat, bottom_mat, color_vec);
% 


%% gel B5_3: align profiles

cd(path_to_data_str)
load("210727_gain_gel_B5_3_processed");
load("201119_gain_no_gel_processed")
cd(path_back)

no_gel = cropProfile(no_gel, 'bottom', 0.1, 'mm');
no_gel = cropProfile(no_gel, 'left', 1.1, 'mm');
no_gel = cropProfile(no_gel, 'top', 1, 'mm');
no_gel = cropProfile(no_gel, 'right', 0.2, 'mm');
gel = cropProfile(gel, 'top', 0.9, 'mm');

%subtract out estimated floor
ground_swath = [4,5,12,13];
no_gel_bottom = sectionbyMm(no_gel, ground_swath);
gel_bottom = sectionbyMm(gel, ground_swath);

no_gel_bottom_mean = mean(no_gel_bottom, 'all');
gel_bottom_mean = mean(gel_bottom, 'all');
no_gel.profile = no_gel.profile -(no_gel_bottom_mean-gel_bottom_mean);

gel.name = "Gel B5_2 Gain Stim";

if ~checkSizeMatch(gel, no_gel)
    [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
    [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
end

location_mat = vertcat( [7, 7.2, 10.5, 11], [1.5, 1.8, 9.4, 9.8], ...
    [1.5, 1.6, 12.6, 12.7], [4.1, 4.5, 6.4, 6.8], [4.15, 4.25, 0.3, 0.8]);
bottom_mat = vertcat( [4, 4.5, 10.5, 11], [3.5, 4, 9, 9.8], ...
   [3, 3.5, 12.6, 12.7], [6, 6.5, 6, 6.5], [2, 2.5, 1, 2]);
location_mat(:, 1:2) = location_mat(:, 1:2)-0.6;
location_mat(:, 3:4) = location_mat(:, 3:4)+0.1;
bottom_mat(:, 1:2) = bottom_mat(:, 1:2)-0.5;
bottom_mat(:, 3:4) = bottom_mat(:, 3:4)+0.1;
[scale_factor_b5_2, r_gelb5_2] = gain_process(gel, no_gel, location_mat, bottom_mat, color_vec);


%% Older gels

load('201204_gain_35_gel_processed');
new_gel = gel;
new_gel.name = "Gel 35 Gain Stim";
load('201119_gain_no_gel_processed');
load('201119_gain_gel_processed');
gel.name = "Gel 36 Gain Stim";
cd ..

%% Align profiles
gel = cropProfile(gel, 'bottom', 3, 'px');
no_gel = cropProfile(no_gel, 'left', 3, 'px');
no_gel = cropProfile(no_gel, 'right', 4, 'px');
no_gel = cropProfile(no_gel, 'left', 0.026, 'mm');
gel = cropProfile(gel, 'right', 0.025, 'mm');


%subtract out estimated floor
ground_swath = [4,5,12,13];
no_gel_bottom = sectionbyMm(no_gel, ground_swath);
gel_bottom = sectionbyMm(gel, ground_swath);

no_gel_bottom_mean = mean(no_gel_bottom, 'all');
gel_bottom_mean = mean(gel_bottom, 'all');
no_gel.profile = no_gel.profile -(no_gel_bottom_mean-gel_bottom_mean);


if ~checkSizeMatch(new_gel, new_no_gel)
    [new_gel, new_no_gel] = resampleToMin(new_gel, new_no_gel); %resamples to the min resolution
    [new_gel, new_no_gel] = bruteCropFit(new_gel, new_no_gel); %crops to same size
end

%% Finding height differences at peaks


location_mat = vertcat( [6, 6.2, 8.5, 9], [6.1, 6.2, 1.3, 1.4], ...
    [0.6, 0.8, 7.3, 7.6], [3.3, 3.5, 4.2, 4.4], [0.6, 0.7, 4.3, 4.5]);
bottom_mat = vertcat( [3.5, 4, 9, 9.5], [7, 7.4, 1, 2], ...
    [2.5, 3, 7, 7.5], [5, 5.5, 4, 4.5], [1.5, 2, 4, 4.5]);
[scale_factor] = gain_process(new_gel, new_no_gel, location_mat, bottom_mat, color_vec);



%% For new gel
location_mat = vertcat( [7.6, 7.9, 12.2, 12.5], [2.1, 2.4, 10.7, 10.9], ...
    [4.9, 5.2, 7.6, 7.9], [7.7, 7.8, 4.7, 4.8], [2.2, 2.3, 13.7, 13.8]);
bottom_mat = vertcat( [8.4, 8.8, 8.5, 9], [4, 4.5, 11, 11.5], ...
    [6.5, 7, 6, 6.5], [8.5, 9, 5, 5.5], [3, 3.5, 12.5, 13]);

[scale_factor] = gain_process(gel, no_gel, location_mat, bottom_mat, color_vec);


