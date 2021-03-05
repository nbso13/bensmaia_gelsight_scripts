%% Main Script for Analysis
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry
% Date: June 26 2020
% cd ~/Documents/bensmaia_lab/bensmaia_gelsight_scripts/profilometry_analysis_scripts
clear 
close all

%% set vars
ppm = 10;
save_figures = 0;
texture_name = "sueded_cuddle";
figure_dir = strcat('../../pngs/feb_23_charles_checkin/', texture_name);
gel_constant = 1.48;
log = 1; %yes/no we want the color map to be log scale
cax = "max"; %we want the max val to be the range
max_freq = 5; %1 dot per mm is upper freq limit
one_dim = 0; % yes, this is one dimensional and grating goes horizontal.

%% Texture options

% COMPLIANT
% % sueded cuddle
% filename_gel = "210219_sueded_cuddle_gel_7_processed";
% filename_nogel = "210222_sueded_cuddle_no_gel_processed";

% wool_blend
% filename_gel = "210217_wool_blend_gel_7_processed";
% filename_nogel = "210216_wool_blend_no_gel_processed";

% hucktowel
% filename_gel = "210209_hucktowel_gel_11_processed";
% filename_nogel = "210204_hucktowel_nogel_processed";

%velvet REDO
% filename_gel = "210122_velvet_gel_3_processed";
% filename_gel = "210122_velvet_gel_4_processed";
% filename_nogel = "210121_velvet_no_gel_processed";

% CORDUROY REDO
% filename_gel = "201118_corduroy_35_gel_trimmed";
% filename_nogel = "201118_corduroy_no_gel_trimmed";

% %upholstery 1 on gel 1
% filename_gel = "210113_upholstry_36_gel_1_processed";
% filename_nogel = "210113_upholstry_no_gel_processed";

%upholstery gel 2
% filename_gel = "210112_upholstry_36_gel_2_processed";
% filename_nogel = "210111_upholstery_no_gel_processed";

%upholstery gel 1
% filename_gel = "210112_upholstry_36_gel_1_processed";
% filename_nogel = "210111_upholstery_no_gel_processed";

% TO DO
% Empire Velveteen
% Taffeta
% Wool Gabardine
% Flag/Banner
% Premier Velvet
% Wool Crepe
% Chiffron
% Swimwear Lining
% Blizzard Fleece
% Drapery Tape(Foam Side)


% NON COMPLIANT

% TO DO
% 12 mm grating
% 8mm grating
% 5mm grating
% 2mm embossed dots
% 3mm embossed dots
% 4mm embossed
% 5
% 6

%cross
% filename_gel = "201119_cross_gel_processed";
% filename_gel = "210219_cross_gel_7_processed";
% filename_nogel = "201119_cross_no_gel_processed";

% gain_stim
% filename_gel = "201119_gain_gel_processed";
% filename_nogel = "201119_gain_no_gel_processed";

% 3/05 DOTS
% filename_gel = "200305_dots_gel_processed";
% filename_nogel = "200305_dots_no_gel_processed";

% 5/24 DOTS
% filename_gel = "200524_dots_gel_processed_aligned";
% filename_nogel = "200524_dots_no_gel_processed";

% 9/23 DOTS
% filename_gel = "200923_dots_gel_processed_aligned";
% filename_nogel = "200923_dots_no_gel";

% 9/25 DOTS
% filename_gel = "200925_dots_gel_processed_aligned";
% filename_nogel = "200925_dots_no_gel";

% 01/19/21 DOTS
% filename_gel = "210119_dots_gel_3_processed";
% filename_nogel = "210120_dots_no_gel_processed";

% 2MM GRATING THIN
% filename_gel = "201116_2mm_grating_35_gel_processed";
% filename_nogel = "201019_no_gel_2mm_grating";

% 1MM GRATING THIN
% filename_gel = "201116_1mm_grating_35_gel_processed";
% filename_nogel = "201021_1mm_grating_no_gel";


% 2/16/21 3mm grating,
% filename_gel = "210216_3mm_grating_gel_7_processed"; %gel 7
% filename_gel = "210216_3mm_grating_gel_11_processed"; %gel 11
% filename_nogel = "210212_3_mm_grating_no_gel_processed";


% 1MM GRATING 2/23
% filename_gel = "210223_1mm_grating_gel_7_processed";
% filename_gel = "210223_1mm_grating_gel_11_processed";
% filename_nogel = "210115_1mm_grating_no_gel_processed";


% 1MM GRATING 2/26
% filename_gel = "210226_1mm_grating_gel_7_100_grams_processed";
% filename_gel = "210226_1mm_grating_gel_7_200_grams_processed";
% filename_nogel = "210115_1mm_grating_no_gel_processed";

% 3MM GRATING 2/23
% filename_gel = "210226_3mm_grating_gel_7_100_grams_processed";
filename_gel = "210226_3mm_grating_gel_7_200_grams_processed";
filename_nogel = "210212_3_mm_grating_no_gel_processed";



%% Load data process data
cd ../../mat_files/
load(filename_gel);
load(filename_nogel);
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts

if ~isfield(gel, 'scaled')
    gel.profile = gel.profile.*gel_constant; %scale up
end
% 
% gel = rotateProfilometry(gel, 270);
% no_gel = rotateProfilometry(no_gel, 270);

if ~checkSizeMatch(gel, no_gel)
    [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
    [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
end

figure
visualizeProfile(gel);
figure
visualizeProfile(no_gel);

%% generate touchsim models
% gel_area = gel.x_axis(end)*gel.y_axis(end);

plot_flag = 1;
pin_radius = 0.025;
[new_gel_ts, new_no_gel_ts, skin_surface_ts, ...
    gel_fig, no_gel_fig, skin_fig] = TouchSimSkin(gel, no_gel, ppm, pin_radius, plot_flag);

%save surfs
if save_figures
    cd(figure_dir)
    cd('gel')
    temp = char(filename_gel);
    date_gel = temp(1:6);
    saveas(gel_fig, strcat("surf_gel_", texture_name, "_", string(date), '.png'));
    cd ..
    cd('real')
    temp = char(filename_nogel);
    date_no_gel = temp(1:6);
    saveas(no_gel_fig, strcat("surf_real_", texture_name, "_", string(date), '.png'));
    cd ..
    cd('ts')
    saveas(skin_fig, strcat("surf_ts_", texture_name, "_", string(date), '.png'));
    cd ../../../../bensmaia_gelsight_scripts/profilometry_analysis_scripts %out of ts, hucktowel, _checkin, pngs,
end


%get profiles
touchsim_gel = shape2profilometry(skin_surface_ts.shape, ...
    skin_surface_ts.offset, skin_surface_ts.pins_per_mm);
new_gel = shape2profilometry(new_gel_ts.shape, ...
    new_gel_ts.offset, new_gel_ts.pins_per_mm);
new_no_gel = shape2profilometry(new_no_gel_ts.shape, ...
    new_no_gel_ts.offset, new_no_gel_ts.pins_per_mm);


% [~, ~, fig_handles_gel] = characterizeFilterNew(new_no_gel, new_gel, vline, hline, plot_flag, [1, 2]);
% [~, ~, fig_handles_ts] = characterizeFilterNew(new_no_gel, touchsim_gel, vline, hline, plot_flag, [1, 2]);
%save freq figures
% if save_figures
%     cd(figure_dir)
%     cd('gel')
%     temp = char(filename_gel);
%     date_gel = temp(1:6);
%     saveas(fig_handles_gel(1), strcat("gel_filt_", texture_name, "_", string(date), '.png'));
%     saveas(fig_handles_gel(2), strcat("gel_freq_", texture_name, "_", string(date), '.png'));
%     cd ..
%     cd('ts')
%     temp = char(filename_nogel);
%     date_no_gel = temp(1:6);
%     saveas(fig_handles_ts(1), strcat("ts_filt_", texture_name, "_", string(date), '.png'));
%     saveas(fig_handles_ts(2), strcat("ts_freq_", texture_name, "_", string(date), '.png'));
%     cd ../../../../bensmaia_gelsight_scripts/profilometry_analysis_scripts %out of ts, hucktowel, _checkin, pngs,
% end

% COMPARE GRATING CONFORMANCE
[diff1, trace, axis, fig_handle_gel] = displayDirecAverage(new_gel, 'v');
title("3mm Grating, Gel 7, 200 grams")
[diff2, ~, ~, fig_handle] = displayDirecAverage(touchsim_gel, 'v');
disp(strcat("difference in maxima: ", num2str(diff2-diff1)));

%show the profiles
% figure
% visualizeProfile(touchsim_gel);
% figure
% visualizeProfile(new_gel);
% figure
% visualizeProfile(new_no_gel);

cd ../touchsim_gelsight
setup_path;
cd ../profilometry_analysis_scripts/

%% pressure matching procedure
% amplitudes = 0:0.2:1;
% gel_mass = 200; %200 grams used
% gel_area = 27*12; %mm sq
% plot_flag = 1;
% ts_struct = skin_surface_ts;
% pressures = ampCurve(ts_struct, pin_radius, 12*20, gel_mass, amplitudes, plot_flag);

%% Neural Data
% amps = 0.05:0.05:0.4;
% gel_frs = zeros(length(amps), 6);
% ts_frs = zeros(length(amps), 6);
plot_flag = 1;
% for k = 1:length(amps)
aff_pop = affpop_hand('D2d', 0.8);

skin_surface_ts.amp = max(skin_surface_ts.offset); %max(skin_surface_ts.offset);
% skin_surface_ts.amp = 1.95; %1/19 dots
new_gel_ts.amp = max(new_gel_ts.offset); %max(new_gel_ts.offset);
new_no_gel_ts.amp = 0;
ts_structs = [skin_surface_ts, new_gel_ts];
speed = 80; %mm/s.
len = 0.4; % s
loc = [0 0];
samp_freq = 2000; % hz
ramp_len = 0.01;
%% calc_responses
[FRs, r, a, figure_handles] = calcResponses(ts_structs, aff_pop, ppm, speed, len, loc, samp_freq, ramp_len, plot_flag);
% ts_frs(k,:) = [FRs{1,4}', FRs{1,5}'];
% gel_frs(k,:) = [FRs{2,4}', FRs{2,5}'];
% end


if save_figures
    cd(figure_dir)
    direcs = ["ts", "gel"];
    dates = [string(date_no_gel), string(date_gel)];
    for i = 1:size(figure_handles, 1)
        cd(direcs(i))
        %saveas(figure_handles{i,1}, strcat("stim_", direcs(i), "_", texture_name, "_", dates(i), '.png'));
        saveas(figure_handles{i,2}, strcat("response_", direcs(i), "_", texture_name, "_", dates(i), '.png'));
        saveas(figure_handles{i,3}, strcat("mean_responses_", direcs(i), "_", texture_name, "_", dates(i), '.png'));
        saveas(figure_handles{i,4}, strcat("force_profile_", direcs(i), "_", texture_name, "_", dates(i), '.png'));
        saveas(figure_handles{i,5}, strcat("trace_profile_", direcs(i), "_", texture_name, "_", dates(i), '.png'));
        cd ..
    end
    cd ../../../bensmaia_gelsight_scripts/profilometry_analysis_scripts %out of ts, hucktowel, _checkin, pngs,
end

%% characterize filters
% plot force profiles


plot_flag = 1;
vline = 10;
hline = 10;
cd ../profilometry_analysis_scripts/
characterizeFilterNew(new_no_gel, new_gel, vline, hline, plot_flag, [1, 2]);
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
%% compare to calculated filter
% sigma = 200;
% %filtered_dots_profile = gelsightFilterFull(no_gel, gel_amp_ratio_mat, gel_amp_ratio_fx, gel_amp_ratio_fy);
% %testing filter_characterization
% filtered_dots_profile = imgaussfilt(no_gel.profile, sigma);
% filt_dots = no_gel;
% filt_dots.profile = real(filtered_dots_profile);
% plot_flag = 1;
% characterizeFilterFull(no_gel, filt_dots, vline, hline, cax, log, plot_flag, max_freq, one_dim);

%% Analyze Empirical Data and Characterize Filter
vline = 1;
hline = 1;
log = 1;
cax = "max";
max_freq = 2;
plot_flag = 1;
% [gel_amp_ratio_mat, gel_amp_ratio_fx, gel_amp_ratio_fy] = characterizeFilterFull(no_gel, ...
%     gel, vline, hline, cax, log, plot_flag, max_freq, one_dim);


% Thesis thinking;

% Characterizing the action of the skin.
% Need to reread skin mechanics papers. And touchsim. And Gelsight. and tribology.
% 
% 1. Is it well approximated by a linear filter? Compare touchsim and gelsight frequency responses. freqz2
% 2. Is the filter linear or non-linear?
% 3. space variant or space invariant?
% 4. Causal or noncausal.
% 5. IIR or FIR?
% 6. Determine transfer function. OR determine convolution of the space-domain input with the filter's impulse response.
% 2. How about a gaussian filter?
% 2. Examine different textures and fourier transforms.
% 5. Give some impulse.
% 6. Give some step. input

