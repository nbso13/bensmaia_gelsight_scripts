%% Main Script for Analysis
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry
% Date: June 26 2020

clear 
close all

%% set vars
gel_constant = 1.5;
vline = 600;
hline = 600;
plotflag = 0;
log = 1; %yes/no we want the color map to be log scale
cax = "max"; %we want the max val to be the range
max_freq = 5; %1 dot per mm is upper freq limit
one_dim = 0; % yes, this is one dimensional and grating goes horizontal.
% filename_gel = "201118_corduroy_35_gel_trimmed";
% filename_nogel = "201118_corduroy_no_gel_trimmed";
% filename_gel = "201116_2mm_grating_35_gel_processed";
% filename_nogel = "201019_no_gel_2mm_grating";
filename_gel = "201116_1mm_grating_35_gel_processed";
filename_nogel = "201019_no_gel_1mm_grating";

%% Generate sim texture
% height = 450;
% diam = 600;
% res = 10; %microns
% freq = 0.5; %dots per mm
% windowsize = 21000; %microns
% dot_pattern = generate_texture("dots", freq, height, diam, res, windowsize);
% visualizeProfile(dot_pattern);
% plot_flag = 1;
% [ind, sa] = surfaceArea4Ind(dot_pattern, plot_flag);
%

%% Load data process data
cd ../../mat_files/
load(filename_gel);
load(filename_nogel);
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts

gel.profile = gel.profile.*gel_constant; %scale up

if ~checkSizeMatch(gel, no_gel)
    [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
    [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
end

figure
visualizeProfile(gel);
figure
visualizeProfile(no_gel);

% if filename_gel == "10x_Gel_Dot_200524_processed_aligned"
%     gel = gel_dot_200524;
%     gel.profile = gel.profile./1000;
%     gel.x_axis = gel.x_axis/1000;
%     gel.y_axis = gel.y_axis/1000;
%     no_gel = nogel_dot_200524;
%     no_gel.profile = no_gel.profile./1000;
%     no_gel.x_axis = no_gel.x_axis/1000;
%     no_gel.y_axis = no_gel.y_axis/1000;
% end
% 
% if filename_gel == "10x_Gel_Dot_200923_processed_aligned"
%     gel.profile = gel.profile - 0.25;
%     gel.profile(gel.profile<0)=0;
% end
%plot_flag = 1;
%surfaceArea4Ind(no_gel, plot_flag);


%% Prepare for TouchSim
% pins_per_mm = 8;
% [shape, offset] = profilometry2shape(no_gel, pins_per_mm);
% no_gel_ts = struct;
% no_gel_ts.shape = shape;
% no_gel_ts.offset = offset;
% %save("no_gel_ts", "no_gel_ts");

%% Analyze Empirical Data and Characterize Filter
[gel_amp_ratio_mat, gel_amp_ratio_fx, gel_amp_ratio_fy] = characterizeFilterFull(no_gel, ...
    gel, vline, hline, cax, log, plotflag, max_freq, one_dim);

%% compare to calculated filter
% sigma = 200;
% %filtered_dots_profile = gelsightFilterFull(no_gel, gel_amp_ratio_mat, gel_amp_ratio_fx, gel_amp_ratio_fy);
% %testing filter_characterization
% filtered_dots_profile = imgaussfilt(no_gel.profile, sigma);
% filt_dots = no_gel;
% filt_dots.profile = real(filtered_dots_profile);
% plot_flag = 1;
% characterizeFilterFull(no_gel, filt_dots, vline, hline, cax, log, plot_flag, max_freq, one_dim);


%% generate touchsim models
ppm = 7;
plot_flag = 1;
[new_gel, touchsim_gel, new_no_gel] = TouchSimSkin(gel, no_gel, ppm, plot_flag);

%show the profiles
figure
visualizeProfile(touchsim_gel);
figure
visualizeProfile(new_gel);
figure
visualizeProfile(new_no_gel);


[shape, offset] = profilometry2shape(new_no_gel, ppm);
no_gel_ts = struct;
no_gel_ts.shape = shape;
no_gel_ts.offset = offset;
no_gel_ts.name = "no_gel_ts";

[shape, offset] = profilometry2shape(new_gel, ppm);
gel_ts = struct;
gel_ts.shape = shape;
gel_ts.offset = offset;
gel_ts.name = "gel_ts";

[shape, offset] = profilometry2shape(touchsim_gel, ppm);
touchsim_gel_ts = struct;
touchsim_gel_ts.shape = shape;
touchsim_gel_ts.offset = offset;
touchsim_gel_ts.name = "touchsim_gel_ts";

a = affpop_hand('D2d', 0.4, 'SA1');

cd ../touchsim
setup_path;

ts_structs = [touchsim_gel_ts, gel_ts, no_gel_ts];
speed = 80; %mm/s.
len = 1; % s
loc = [0 0];
samp_freq = 200; % hz
ramp_len = 0.2;

for i = 1:length(ts_structs)
    amp = max(ts_structs(i).offset);
    s = stim_scan_shape(ts_structs(i).shape, ts_structs(i).offset, ppm, len, samp_freq, amp, speed, gel_flag);
    figure
    plot(s)
    r = a.response(s);
    %take out neurons that fire less than 2 spikes per second
    r_new = excludeNeurons(r, 2);
    figure
    plot(r_new)
    title(ts_structs(i).name);
end

%% characterize filters

plot_flag = 1;
vline = 10;
hline = 10;

characterizeFilterFull(new_no_gel, touchsim_gel, vline, hline, cax, log, plot_flag, max_freq, one_dim);
characterizeFilterFull(new_no_gel, new_gel, vline, hline, cax, log, plot_flag, max_freq, one_dim);

%% Make Simulated dots
ratio_mat_size = size(gel_amp_ratio_mat);
window_size = 1000;
dot_height = [25, 20, 15, 10];
dot_freq = 1:30;
%sigma = [1, 5, 10, 15, 20, 50];
sigma = 1;
amp_avs = {};
plotflag = 0;
freq_limit = 30;

%get size of amp ratio mats and power spec mats
dot_pattern = generate_texture("dots", 4, 100, 100, window_size);
filtered_dots_profile = imgaussfilt(dot_pattern.profile, 1);
filt_dots = dot_pattern;
filt_dots.profile = real(filtered_dots_profile);
[amp_ratio_mat, ~, ~] = characterizeFilterFull(dot_pattern, filt_dots, 1, 1, cax, log, plotflag);
rat_mat_size = size(amp_ratio_mat); %what size is the ratio mat?
[power_mat, ~, ~] = powerSpectrumFull(dot_pattern);
power_mat_size = size(power_mat);
power_spec_sum = zeros(power_mat_size); 
add_count_power = 0;

%averageing amp ratio mats across all textures for each different sigma.
figure;
for j = 1:length(sigma)
    amp_mat_sum = zeros(rat_mat_size); 
    add_count =0; %count how many times, so we can avearge later

  
    for k = 1:length(dot_height) %for every height
        for i =  1:length(dot_freq) %for every frequency
            if dot_height(k)*2*dot_freq(i) > 1000 %if the size and freq would mean overcrowding
                continue %skip
            end
            % generate texture
            dot_pattern = generate_texture("dots", dot_freq(i), dot_height(k), window_size);
            %filtered_dots_profile = gelsightFilterFull(dot_pattern, gel_amp_ratio_mat, gel_amp_ratio_fx, gel_amp_ratio_fy);
            filtered_dots_profile = imgaussfilt(dot_pattern.profile, sigma(j));
            filt_dots = dot_pattern;
            filt_dots.profile = real(filtered_dots_profile);
            [amp_ratio_mat, fx, fy] = characterizeFilterFull(dot_pattern, filt_dots, 1, 1, cax, log, plotflag);
            amp_mat_sum = amp_mat_sum + amp_ratio_mat; %summing to total to average later
            add_count = add_count+1;
            
            %go thru textures and get average power spectra
            if j == 1
                %power Spectrum
                [power_spec_mat, fpx, fpy] = powerSpectrumFull(dot_pattern);
                power_spec_sum = power_spec_sum + power_spec_mat;
                add_count_power = add_count_power+1;
            end
            
        end
    end
    %average
    amp_avs{j} = amp_mat_sum./add_count;
    
    %plot
    subplot(2,3,j)
    fx = fx(fx<freq_limit); %we only care about some frequencies
    fy = fy(fy<freq_limit);
    amp_avs{j} = amp_avs{j}(1:length(fy), 1:length(fx)); 
    
    imagesc(fx, fy, amp_avs{j})
    colorbar;
    caxis([0 2]);
    ax = gca;
    xlabel("Frequency (1/mm)")
    ylabel("Frequency (1/mm)")
    ax.YDir = 'normal';
    title(strcat("Average Amp Ratio for Sigma = ", num2str(sigma(j))))
end
set(gcf, 'position', [100 100 1300 800]);

%average
power_average = power_spec_sum./add_count_power;
%plot power spectrum average
figure
fpx = fpx(fpx<freq_limit); %we only care about some frequencies
fpy = fpy(fpy<freq_limit);
fpx = fpx(2:end);
fpy = fpy(2:end);
power_average = power_average(2:length(fpy), 2:length(fpx)); 

imagesc(fpx, fpy, power_average)
c=colorbar;
ylabel(c, 'Power (mm^2/Hz)');
%calculate max for display
maxo = max(power_average, [], 'all');
caxis([0 maxo]);
ax = gca;
%set(ax,'ColorScale','log')
xlabel("Frequency (1/mm)")
ylabel("Frequency (1/mm)")
ax.YDir = 'normal';
title("Average Power Spectrum for all textures");
set(gcf, 'position', [100 100 1300 800]);

%% testing different textures
i = 3;
plotflag = 1;
switch i
    case 1
        dot_pattern = generate_texture("flat", dot_freq, dot_height, window_size);
    case 2
        dot_pattern = generate_texture("unit_impulse", dot_freq, dot_height, window_size);
    case 3
        dot_pattern = generate_texture("unit_step", dot_freq, dot_height, window_size);
end
% filter texture based on empirical gelsight filter
filtered_dots_profile = gelsightFilterFull(dot_pattern, gel_amp_ratio_mat, gel_amp_ratio_fx, gel_amp_ratio_fy);
filt_dots = dot_pattern;
filt_dots.profile = real(filtered_dots_profile);
hline = findBestLine(filt_dots.profile, 'horiz');
vline = findBestLine(filt_dots.profile, 'vert');
characterizeFilterFull(dot_pattern, filt_dots, vline, hline, cax, log, plotflag);

