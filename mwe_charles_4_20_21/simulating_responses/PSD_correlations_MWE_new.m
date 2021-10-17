%% Main Script for Analysis
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry

clear
close all
local_data_path_str = "../../../mwe_data/";
local_path_back = "../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/simulating_responses";
addpath("helper_functions")

%% set vars

cd(strcat(local_data_path_str, "neural_data"))
load("RawPAFData")
load("TextureNames")
cd(local_path_back)


load('colorscheme');
scatter_size = 50;


good_neurons = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 22 25 28 33 34];
texture_nums = [31, 25, 45, 7, 4, 33, 9, 21, 38, 44, 6, 15, 2, 43, 20, 1, 12];

% texture_nums = [45, 4, 31, 42, 7]; %100 gram textures
% texture_nums = [25, 45, 7];
% texture_nums = [25, 45, 50];
% texture_nums = [49, 50, 55];

%pull names
for i = 1:length(htxt_name)
    texture_names(i) = string(htxt_name{i});
end

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
        
num_textures = length(filename_gel);
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

%HYPERPARAMS
stopBand = 0.3; %frequencies below 0.5 are noise
time_samp_period = 0.001; %millisecond time resolution
speed = 80; %mm/s

%% pull activities and real rates for each texture 

my_texture_names = texture_names(texture_nums);

disp(texture_names(texture_nums));
if ~((length(texture_nums) == num_textures) && (num_textures == length(filename_nogel)))
    error("Filenames and/or Texture numbers don't match up.")
end

neuron_identities = {iPC, iRA, iSA};
excludeNeurons = 1; %don't average neurons that don't fire
[activities, spike_trains, space_vec] = pullRealActivities(rates, spikes, ...
    my_texture_names, good_neurons, neuron_identities, texture_nums, ...
    speed, excludeNeurons, time_samp_period);


%% Run Loop
correlations = {};
gel_no_gel_freq_ratios = {};
gel_no_gel_freq_ratio_norms = {};
f_no_gels = {};
% textures x afferent classes array([gel no_gel] x trial)

%Correlating
disp(strcat("Calculating spike trains PSD"));
[spike_psds, f_rate] = rateSpectralAnalysis(spike_trains, space_vec);

tic
for i = 1:num_textures
    %load
    cd(strcat(local_data_path_str, "sim_data"))
    load(filename_gel(i), "gel");
    load(filename_nogel(i), "no_gel");
    cd(local_path_back)
    
    disp(strcat("Highpass filter at ", num2str(stopBand), " per mm."));
    gel = removeLowFreq(gel, stopBand, 'charles');
    no_gel = removeLowFreq(no_gel, stopBand, 'charles');
    
    %profile and gel power spectra
    disp(strcat("Calculating power spectra of profiles"));
    [gel_psd, f_gel, no_gel_psd, f_no_gel, gel_to_nogel_ratio, ...
        gel_to_nogel_ratio_norm] = profSpectralAnalysis(gel, no_gel); %ratio freq axis is no gel.
    gel_no_gel_freq_ratios{i} = gel_to_nogel_ratio;
    gel_no_gel_freq_ratio_norms{i} = gel_to_nogel_ratio_norm;
    f_no_gels{i} = f_no_gel;
    
    %Correlating
    disp(strcat("Calculating correlations of rates to profile spectra"));
    
    aff_psds = {spike_psds{i, 1}, spike_psds{i, 2}, spike_psds{i, 3}};
    
    [cor] = spike_train_profile_corr(no_gel.name, gel_psd, f_gel, no_gel_psd, f_no_gel, aff_psds, f_rate); 
    % cor returned is aff_class x [gel no_gel] x max_trial num. Non used
    % trials are nans.
    correlations = vertcat(correlations, cor); %afferent classes array([gel no_gel] x trial)
    close all
end
total_time = toc;
disp(strcat("average time per texture: ", num2str(total_time/num_textures)))


%% find and plot mean ratio
mean_ratios = gel_no_gel_freq_ratios{1}';
mean_ratios_norm = gel_no_gel_freq_ratio_norms{1}';
for i = 2:num_textures
    interp_mean = interp1(f_no_gels{i}', gel_no_gel_freq_ratios{i}', f_no_gels{1}');
    mean_ratios = [mean_ratios; interp_mean];
    
    interp_mean_norm = interp1(f_no_gels{i}', gel_no_gel_freq_ratio_norms{i}', f_no_gels{1}');
    mean_ratios_norm = [mean_ratios_norm; interp_mean_norm];
end

mean_ratio = nanmean(mean_ratios, 1);

mean_ratio_norm = nanmean(mean_ratios_norm, 1);

f = figure;
f.Position = [100, 100, 300, 300];
plot(f_no_gels{1}, mean_ratio, 'k');
ylabel('Mean Raw Ratio'); 
yyaxis right
plot(f_no_gels{1}, mean_ratio_norm);
ylabel('Mean Norm Ratio'); 
xlabel('Spatial Frequency (1/mm)');
xlim([0 7]);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';


%% plotting correlations

% get mean corr and sd corr for each afferent for each texture for each gel
mean_corr_gel = zeros(size(correlations));
sd_corr_gel = mean_corr_gel;
mean_corr_no_gel = mean_corr_gel;
sd_corr_no_gel = mean_corr_gel;
p_vals = mean_corr_gel;
h_vals = mean_corr_gel;

for i = 1:size(correlations, 1)
    for j = 1:size(correlations, 2)
        corrs = correlations{i,j};
        mean_corr_gel(i,j) = mean(corrs(:,1)); % should instead look at best...
        mean_corr_no_gel(i,j) = mean(corrs(:,2));
        sd_corr_gel(i,j) = std(corrs(:,1));
        sd_corr_no_gel(i,j) = std(corrs(:,2));
        [h, p] = ttest(corrs(:,1), corrs(:,2));
        p_vals(i,j) = p;
        h_vals(i,j) = h;
    end
end

%% plotting corrs across textures


% more appropriate scatter plot

f = figure;
f.Position = [100, 100, 900, 400];
aff_names = ["PCs", "RAs", "SAs"];
sig_mat_textures = zeros(3, num_textures);
for i = 1:3 %plotting either scatter or bar
    subplot(1,3,i);
    hold on
    for j = 1:num_textures
        scatter(mean_corr_no_gel(j, i), mean_corr_gel(j, i), scatter_size, colorscheme(j, :), 'filled'); % touch_sim x, gel y
        er = errorbar(mean_corr_no_gel(j, i), mean_corr_gel(j, i), sd_corr_no_gel(j, i), 'horizontal', 'HandleVisibility','off'); %real
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er = errorbar(mean_corr_no_gel(j, i), mean_corr_gel(j, i), sd_corr_gel(j, i), 'HandleVisibility','off');
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    xlim([-1 1]);
    ylim([-1 1]);
    plot([-1, 1], [-1, 1], 'k');
    xlabel("No Gel PSD Correlations")
    ylabel("Gel PSD Correlations")
    title(aff_names(i));
    ax = gca;
    daspect([1 1 1]);
    ax.YAxis.Visible = 'off';
    if i == 1
        strs = my_texture_names';
        colors = colorscheme(1:size(strs,1), :);
        leg = legend([color_legend(strs, colors), "unity line"]);
        leg.Box = 0;
        ax.YAxis.Visible = 'on';
    end
    ax.FontSize = 12;
    ax.FontWeight = 'bold';
    
    
end


% attempts at a bar plot


PC_means = [mean_corr_gel(:,1), mean_corr_no_gel(:,1)];
RA_means = [mean_corr_gel(:,2), mean_corr_no_gel(:,2)];
SA_means = [mean_corr_gel(:,3), mean_corr_no_gel(:,3)];
aff_means = {PC_means, RA_means, SA_means};
PC_sds = [sd_corr_gel(:,1), sd_corr_no_gel(:,1)];
RA_sds = [sd_corr_gel(:,2), sd_corr_no_gel(:,2)];
SA_sds = [sd_corr_gel(:,3), sd_corr_no_gel(:,3)];
aff_sds = {PC_sds, RA_sds, SA_sds};
aff_names = ["PC PSDs", "RA PSDs", "SA PSDs"];

figure;
for i = 1:3
    subplot(3,1,i);
    aff_mean_mat = aff_means{i};
    % bars
    b = bar(aff_mean_mat);
    title(aff_names(i))
    xticklabels(my_texture_names);
    ylabel("Correlation Coefficient")
    ylim([-1 1]);
    %errorbars
    hold on
    
    errorbars_group(aff_mean_mat, aff_sds{i}, p_vals(:, i));
end


function [] = textBar(Bar, text_array)
    for k1 = 1:length(Bar)
        ctr(k1,:) = bsxfun(@plus, Bar(k1).XData, Bar(k1).XOffset');   
        ydt(k1,:) = Bar(k1).YData;                                        % Individual Bar Heights
    end
   
    for k1 = 1:size(ctr,2)
        text(ctr(:,k1), ydt(:,k1), text_array(k1), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
    end
end

function [] = errorbars_group(bar_data, error_data, significance)
    % Find the number of groups and the number of bars in each group
    [ngroups, nbars] = size(bar_data);
    % Calculate the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    
    % for text
    x = (1:ngroups) - groupwidth/2 + (2*1.5-1) * groupwidth / (2*nbars); %only works for 2 bar groups
    y = max(bar_data, [],  2);

    for i = 1:length(x)
        if significance(i) < 0.05
            if y(i) < 0
                text(x(i), y(i)/6, "*", 'FontSize',14, 'FontWeight', 'bold');        
            else
                text(x(i), y(i)+ y(i)/6, "*", 'FontSize',14, 'FontWeight', 'bold');        
            end
        end
    end
    
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, bar_data(:,i), error_data(:,i), 'k', 'linestyle', 'none');
    end
end
