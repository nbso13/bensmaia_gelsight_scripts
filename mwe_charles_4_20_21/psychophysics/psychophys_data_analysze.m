%% psychophysics script
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry
% Date: July  2021

clear
close all


%% main

% load, and set up path
[local_data_path_str, local_path_back, colorscheme, roughData] = clear_close_path_setup();

% set vars
[scatter_size, pins_per_mm, stopBand, pin_radius, plot_flag] = setup_variables();

% choose which texture set to analyze
name_str = 'compliant'; % 'all', 'compliant', 'noncompliant', '100grams', 'test'
[texture_nums, filename_gel, filename_nogel] = assign_names(name_str);

%pull texture names
[texture_names, my_texture_names, rough_means, rough_sds] = pull_texture_names(texture_nums, roughData);

% timer on
tic

%main loop
num_textures = length(filename_gel);
[gel, no_gel, gel_ts, no_gel_ts, skin_surface_ts, ...
    x_axis] = calc_stats(num_textures, filename_gel, ...
    filename_nogel, local_data_path_str, local_path_back, stopBand, ...
    pins_per_mm, pin_radius, plot_flag, my_texture_names);

% timer off
time = toc;
disp(strcat("main loop time per texture: ", num2str(time/num_textures)))

%% plotting
k_len = 5;
plot_psychophys(k_len, x_axis, roughData, texture_nums, num_textures, ...
    scatter_size, colorscheme, my_texture_names)
%% functions

function [texture_nums, filename_gel, filename_nogel, texture_names, ...
    my_texture_names] = assign_names(name_str, htxt_name)
% assigns filenames based on name_str:
switch name_str
    case 'all'
        filename_gel = ["210223_velvet_gel_7_processed", ...
            "201118_corduroy_35_gel_trimmed", ...
            "210226_blizzard_fleece_gel_7_200_grams_processed", ...
            "210217_wool_blend_gel_7_processed", ...
            "210209_hucktowel_gel_11_processed",  ...
            "210423_gel_18_careerwear_flannel_200_grams_processed",...
            "210423_gel_18_thick_corduroy_200_grams_processed", ...
            "210423_gel_18_premier_velvet_200_grams_processed", ...
            "210423_gel_18_velour_200_grams_processed", ...
            "210423_gel_19_snowflake_knitside_200_grams_processed",...
            "210728_gel_B5_3_empire_velveteen_processed", ...
            "210728_gel_B5_3_wool_gabardine_processed",...
            "210729_denim_gel_B5_3_processed", ...
            "210729_snowflake_fuzzyside_gel_B5_3_processed",...
            "210729_microsuede_gel_B5_3_processed", ...
            "210730_satin_gel_B5_3_processed",...
            "210730_chiffon_gel_B5_3_processed",...
            "210216_3mm_grating_gel_7_processed", ...
            "210223_1mm_grating_gel_11_processed", ...
            "210414_2mm_dots_gel_17_200_grams_processed"]; 
        
        filename_nogel = ["210121_velvet_no_gel_processed", ...
            "201118_corduroy_no_gel_trimmed", ...
            "210226_blizzard_fleece_no_gel_processed",...
            "210216_wool_blend_no_gel_processed", ...
            "210204_hucktowel_nogel_processed",  ...
            "210423_no_gel_careerwear_flannel_processed", ...
            "210423_no_gel_thick_corduroy_processed", ...
            "210423_no_gel_premier_velvet_processed", ...
            "210423_velour_no_gel_processed", ...
            "210423_no_gel_snowflake_knitside_processed",...
            "210715_empire_velveteen_no_gel_processed", ...
            "210729_wool_garb_no_gel_processed", ...
            "210729_denim_no_gel_processed", ...
            "210729_snowflake_fuzzyside_no_gel_processed", ...
            "210729_microsuede_no_gel_processed", ...
            "210730_satin_no_gel_processed.mat",...
            "210730_chiffon_no_gel_processed",...
            "210212_3_mm_grating_no_gel_processed", ...
            "201021_1mm_grating_no_gel",...
            "210414_2mm_dots_no_gel_processed"];
        
        texture_nums = [31, 25, 45, 7, 4, 42, 33, 9, 21, 38, 44, 6, 15, 2, 43, 20, 1, 12, 49, 50, 55];
        
    case 'compliant'
        filename_gel = ["210223_velvet_gel_7_processed", ...
            "201118_corduroy_35_gel_trimmed", ...
            "210226_blizzard_fleece_gel_7_200_grams_processed", ...
            "210217_wool_blend_gel_7_processed", ...
            "210209_hucktowel_gel_11_processed",  ...
            "210423_gel_18_careerwear_flannel_200_grams_processed",...
            "210423_gel_18_thick_corduroy_200_grams_processed", ...
            "210423_gel_18_premier_velvet_200_grams_processed", ...
            "210423_gel_18_velour_200_grams_processed", ...
            "210423_gel_19_snowflake_knitside_200_grams_processed", ...
            "210728_gel_B5_3_empire_velveteen_processed", ...
            "210728_gel_B5_3_wool_gabardine_processed", ...
            "210729_denim_gel_B5_3_processed", ...
            "210729_snowflake_fuzzyside_gel_B5_3_processed",...
            "210729_microsuede_gel_B5_3_processed", ...
            "210730_satin_gel_B5_3_processed",...
            "210730_chiffon_gel_B5_3_processed"];
        
        filename_nogel = ["210121_velvet_no_gel_processed", ...
            "201118_corduroy_no_gel_trimmed", ...
            "210226_blizzard_fleece_no_gel_processed",...
            "210216_wool_blend_no_gel_processed", ...
            "210204_hucktowel_nogel_processed",  ...
            "210423_no_gel_careerwear_flannel_processed", ...
            "210423_no_gel_thick_corduroy_processed", ...
            "210423_no_gel_premier_velvet_processed", ...
            "210423_velour_no_gel_processed", ...
            "210423_no_gel_snowflake_knitside_processed", ...
            "210715_empire_velveteen_no_gel_processed", ...
            "210729_wool_garb_no_gel_processed", ...
            "210729_denim_no_gel_processed", ...
            "210729_snowflake_fuzzyside_no_gel_processed", ...
            "210729_microsuede_no_gel_processed", ...
            "210730_satin_no_gel_processed.mat",...
            "210730_chiffon_no_gel_processed"];
        texture_nums = [31, 25, 45, 7, 4, 33, 9, 21, 38, 44, 6, 15, 2, 43, 20, 1, 12]; % adding chiffon, suedeside sueded cuddle, satin, 
        
    case 'noncompliant'
        filename_gel = ["210216_3mm_grating_gel_7_processed", ...
            "210223_1mm_grating_gel_11_processed", ...
            "210414_2mm_dots_gel_17_200_grams_processed"];
        filename_nogel = ["210212_3_mm_grating_no_gel_processed", ...
            "201021_1mm_grating_no_gel",...
            "210414_2mm_dots_no_gel_processed"];
        texture_nums = [49, 50, 55];
        
        
    case '100grams'
        filename_gel = ["210304_blizzard_fleece_gel_11_100_grams_processed", ...
            "210304_hucktowel_gel_11_100_grams_processed", ...
            "210304_velvet_gel_11_100_grams_processed", ...
            "210310_wool_blend_gel_11_100_grams_processed", ...
            "210428_velour_gel_19_100_grams_processed",...
            "210428_thick_corduroy_gel_19_100_grams_processed", ...
            "210428_empire_velveteen_gel_19_100_grams_processed"];
        filename_nogel = ["210226_blizzard_fleece_no_gel_processed",...
            "210204_hucktowel_nogel_processed",  ...
            "210121_velvet_no_gel_processed", ...
            "210216_wool_blend_no_gel_processed", ...
            "210423_velour_no_gel_processed", ...
            "210423_no_gel_thick_corduroy_processed", ...
            "210715_empire_velveteen_no_gel_processed"];
        texture_nums = [45, 4, 31, 7, 38, 9, 6];
    case 'test'
        filename_gel = ["210728_gel_B5_3_empire_velveteen_processed", ...
            "210728_gel_B5_3_wool_gabardine_processed", ...
            "210729_denim_gel_B5_3_processed"];
        filename_nogel = ["210715_empire_velveteen_no_gel_processed", ...
            "210729_wool_garb_no_gel_processed", ...
            "210729_denim_no_gel_processed"];
        texture_nums = [6, 15, 2];
    otherwise
        error("No mode match found - choose test, 100grams, compliant, noncompliant, or all")
end

% filename_gel = ["201118_corduroy_35_gel_trimmed", ...
%     "210226_blizzard_fleece_gel_7_200_grams_processed", ...
%     "210223_1mm_grating_gel_11_processed"];


% filename_nogel = ["201118_corduroy_no_gel_trimmed", ...
%     "210226_blizzard_fleece_no_gel_processed",...
%     "201021_1mm_grating_no_gel"];

% texture_nums = [25, 45, 50];
% texture_nums = [49, 50, 55];
end


function [local_data_path_str, local_path_back, colorscheme, roughData] = clear_close_path_setup()
local_data_path_str = "../../../mwe_data/";
local_path_back = "../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/psychophysics";

cd(char(strcat(local_data_path_str, 'psycho_data')))
load("roughData.mat")
cd(char(local_path_back));

cd ../simulating_responses


local_data_path_str = "../../../mwe_data/";
local_path_back = "/../bensmaia_gelsight_scripts/mwe_charles_4_20_21/simulating_responses";
addpath('helper_functions')
cd helper_functions/touchsim_gelsight
setup_path;
cd ../..

load('colorscheme')
end

function [scatter_size, pins_per_mm, stopBand, pin_radius, plot_flag] = setup_variables()
scatter_size = 50;
pins_per_mm = 8;
stopBand = 0.3; %frequencies below 0.5 are noise
pin_radius = 0.025;
plot_flag = 0;
end

function [texture_names, my_texture_names, rough_means, ...
    rough_sds] = pull_texture_names(texture_nums, roughData)

texture_names = {};
for i = 1:length(roughData.textureNames)
    texture_names{i} = roughData.textureNames{i};
end

my_texture_names = {};
for i = 1:length(texture_nums)
    my_texture_names{i} = texture_names{texture_nums(i)};
end

rough_means = roughData.roughMean(texture_nums);
rough_sds = roughData.roughMean(texture_nums);

end

function [sd_ratio, ...
    comp_ratio, sd, sd_ratio_ts, comp_ratio_ts,...
    raw_gel_sd_ratio, ts_sd, ts_rms, gel_rms, gel_sd] = setup_loop_vars(num_textures)

compressions = zeros(num_textures, 1);
sd_ratio = compressions;
comp_ratio = compressions;
sd = compressions;
sd_ratio_ts = compressions;
comp_ratio_ts = compressions;
raw_gel_sd_ratio = compressions;
ts_sd = compressions;
ts_rms = compressions;
gel_sd = compressions;
gel_rms = compressions;
end

function [gel, no_gel, gel_ts, no_gel_ts, skin_surface_ts, ...
    x_axis] = calc_stats(num_textures, filename_gel, ...
    filename_nogel, local_data_path_str, local_path_back, stopBand, ...
    pins_per_mm, pin_radius, plot_flag, my_texture_names)

[sd_ratio, ...
    comp_ratio, sd, sd_ratio_ts, comp_ratio_ts,...
    raw_gel_sd_ratio, ts_sd, ts_rms, gel_rms,...
    gel_sd] = setup_loop_vars(num_textures);

for i = 1:num_textures
    %load
    disp(strcat("Loading data from ", filename_gel(i)));
    cd(char(strcat(local_data_path_str, "sim_data")));
    load(filename_gel(i), "gel");
    load(filename_nogel(i), "no_gel");
    cd(char(strcat("..", local_path_back)));
    
    %truncate low values 
%     bottom_flag = 1;
%     sd = 1.5;
%     gel = truncateProfile(gel, sd, bottom_flag);
%     no_gel = truncateProfile(no_gel, sd, bottom_flag);
    
    %high pass filter
    disp(strcat("Highpass filter at ", num2str(stopBand), " per mm."));
    gel = removeLowFreq(gel, stopBand, 'charles');
    no_gel = removeLowFreq(no_gel, stopBand, 'charles');
    
    % touchsim skin
    disp("Calculating TouchSim surfaces")
    flip_flag = 0;
    [gel_ts, no_gel_ts, skin_surface_ts] = TouchSimSkin(gel, ...
    no_gel, pins_per_mm, pin_radius, flip_flag, plot_flag);
    gel_ts = shape2profilometry(gel_ts.shape, gel_ts.offset, pins_per_mm);
    skin_surface_ts = shape2profilometry(skin_surface_ts.shape, skin_surface_ts.offset, pins_per_mm);
    
    % show example figures
    if i == 1
        fig = figure;
        subplot(2,2, 1); visualizeProfile(gel_ts);
        title("Gel")
        subplot(2,2,2); visualizeProfile(skin_surface_ts);
        title("TouchSim")
        subplot(2,2,3); visualizeProfile(no_gel);
        title("Raw Profile")
        sgtitle(my_texture_names{i})
    end
    
    %calc compression stats
    disp(strcat("Calculating compression difference."));
    gel_iqr = iqr(gel_ts.profile, 'all');
    ts_iqr = iqr(skin_surface_ts.profile, 'all');
    no_gel_iqr = iqr(no_gel.profile, 'all');
    comp_ratio(i) = gel_iqr/no_gel_iqr;
    comp_ratio_ts(i) = ts_iqr/no_gel_iqr;
    sd_ratio(i) = std(gel_ts.profile(:))/std(no_gel.profile(:));
    raw_gel_sd_ratio(i) = std(gel.profile(:))/std(no_gel.profile(:));
    sd_ratio_ts(i) = std(skin_surface_ts.profile(:))/std(no_gel.profile(:));
    sd(i) = std(no_gel.profile(:));
    ts_sd(i) = std(skin_surface_ts.profile(:));
    gel_sd(i) = std(gel_ts.profile(:));
    gel_rms(i) = rms(gel_ts.profile(:));
    ts_rms(i) = rms(skin_surface_ts.profile(:));
    
    
end

x_axis = {sd, ts_sd, gel_sd, ts_rms, gel_rms, comp_ratio, sd_ratio, comp_ratio_ts, ...
    sd_ratio_ts, raw_gel_sd_ratio};

end

function [] = plot_psychophys(k_len, x_axis, roughData, texture_nums, ...
    num_textures, scatter_size, colorscheme, my_texture_names)

titstrs = ["Raw SD", "TS SD", "Gel SD", "TS RMS", "Gel RMS", "Gel IQR ratio", ...
    "Gel SD ratio", "TS IQR ratio", "TS SD ratio", "Gel Full SD Ratio"];
num_participants = 8;

correlation_results = zeros(length(x_axis), 2); %first entry mean, then SD
all_corrs = zeros(length(x_axis), num_participants);
mean_ratings = mean(roughData.allData, 2);
sd_ratings = std(roughData.allData, [], 2);

for k = 1:k_len
    x_ax = x_axis{k};
    title_str = titstrs(k);
    fig = figure;
    corrs = zeros(size(roughData.allData, 2), 1);
    slopes = corrs;
    for i = 1:size(roughData.allData, 2)
        subplot(3,4, i);
        hold on
        rough_data = roughData.allData(:,i);
        rough_data = rough_data(texture_nums);
        for j = 1:num_textures
            scatter(x_ax(j), rough_data(j), scatter_size, colorscheme(j, :), 'filled');
        end
        xlabel(titstrs(k))
        ylabel("rated roughness")
        ylim([0 2.5])
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        
        if i == 1
            strs = my_texture_names';
            colors = colorscheme(1:size(strs,1), :);
            leg = legend([color_legend(strs, colors)]);
            leg.Box = 'off';
        end
        
        coefficients = polyfit(x_ax, rough_data, 1);
        slopes(i) = coefficients(1);
        xFit = linspace(min(x_ax), max(x_ax), 1000);
        yFit = polyval(coefficients , xFit);
        plot(xFit, yFit, 'k-', 'LineWidth', 2); % Plot fitted line.
        corr = corrcoef(x_ax, rough_data);
        corr_coef = corr(1,2);
        corrs(i) = corr_coef;
        title(strcat("Ppnt. ", num2str(i), " r=", num2str(corr_coef)));
    end
    
    %plot mean ratings
    subplot(3, 4, i+1)
    hold on
    rough_data = mean_ratings;
    rough_data = rough_data(texture_nums);
    error_data = sd_ratings(texture_nums);
    for j = 1:num_textures
        scatter(x_ax(j), rough_data(j), scatter_size, colorscheme(j, :), 'filled');
    end
    er = errorbar(x_ax, rough_data, error_data, 'vertical', 'linestyle', 'none', 'HandleVisibility','off'); %real
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    xlabel(titstrs(k))
    ylabel("mean rated roughness")
    ylim([0 2.5])
    ax = gca;
    ax.FontSize = 12;
    ax.FontWeight = 'bold';
    
    coefficients = polyfit(x_ax, rough_data, 1);
    slopes(i) = coefficients(1);
    xFit = linspace(min(x_ax), max(x_ax), 1000);
    yFit = polyval(coefficients , xFit);
    plot(xFit, yFit, 'k-', 'LineWidth', 2); % Plot fitted line.
    corr = corrcoef(x_ax, rough_data);
    corr_coef = corr(1,2);
    corrs(i) = corr_coef;
    title(strcat("Mean, r=", num2str(corr_coef)));
    
    disp(mean(corrs))
    
    correlation_results(k,1) = nanmean(corrs, 'all');
    correlation_results(k,2) = std(corrs);
    all_corrs(k, :) = corrs;
    disp(strcat(num2str(mean(slopes)), " +/- ", num2str(std(slopes)))); 
    fig.Position = [-113 320 1987 676];
    sgtitle(titstrs(k));
end

sd_strs = titstrs(1:3);
fig = figure;
hold on
x = 1:3;
b = bar(x, correlation_results(1:3, 1));
errorbar(x, correlation_results(1:3, 1), correlation_results(1:3, 2), 'k', 'linestyle', 'none');
title("Mean Corrcoeff For Profile z-axis Standard Deviation")
xticks([1, 2, 3])
xticklabels(sd_strs)
ylabel("Mean+/-SD Correlation Coefficient Across Participants")

for i = 1:3
    bee_swarm(all_corrs(i, :), i, 'r', length(all_corrs(i, :)));
end

fig = figure;
hold on
ts_rms = x_axis{4};
gel_rms = x_axis{5};
title_str = "TS vs Gel RMS";
for j = 1:num_textures
    scatter(ts_rms(j), gel_rms(j), scatter_size, colorscheme(j, :), 'filled');
end
xlabel("TS RMS")
ylabel("Gel RMS")
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';

strs = my_texture_names';
colors = colorscheme(1:size(strs,1), :);
leg = legend([color_legend(strs, colors)]);
leg.Box = 'off';

coefficients = polyfit(ts_rms, gel_rms, 1);
xFit = linspace(min(ts_rms), max(ts_rms), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 2); % Plot fitted line.
corr = corrcoef(x_ax, rough_data);
corr_coef = corr(1,2);
title(strcat("TS vs Gel RMS, r=", num2str(corr_coef)));
fig.Position = [-113 320 1987 676];


    
end

%% scrap code

% fig = figure;
% hold on
% for j = 1:num_textures
%     scatter(x_ax(j), rough_means(j), scatter_size, colorscheme(j, :), 'filled'); 
%     er = errorbar(x_ax(j), rough_means(j), rough_sds(j), 'vertical', 'linestyle', 'none', 'HandleVisibility','off'); %real
%     er.Color = [0 0 0];
%     er.LineStyle = 'none';
% end
% xlabel(tit)
% ylabel("Mean rated roughness")
% title(strcat("Roughness ratings vs ", tit));
% ax = gca;
% ax.FontSize = 12;
% ax.FontWeight = 'bold';
% strs = my_texture_names';
% colors = colorscheme(1:size(strs,1), :);
% leg = legend([color_legend(strs, colors)]);
% leg.Box = 0;
% mdl = fitlm(x_ax, rough_means);
% disp(mdl)
