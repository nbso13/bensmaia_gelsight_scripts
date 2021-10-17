%% Main Script for Analysis
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry
% Date: June 2021
clear
close all

%set loop vars

ppms = [9];
flip_flag = 1;
rads = 0.5*1./ppms;
root_means = {};
pearsons = {};
rati = {};
tpt = [];
m = 1;% for m = 1:length(ppms)

%% main
%setup vars
[local_data_path_str, local_path_back, colorscheme, figure_dir, ...
    aff_colors, neuron_id, rates, spikes, iPC, iRA, iSA, htxt_name] = setup_vars();

% choose which texture set to analyze
name_str = 'compliant'; % 'all', 'compliant', 'noncompliant', '100grams', 'test'
[texture_nums, filename_gel, filename_nogel, texture_names, ...
    my_texture_names] = assign_names(name_str, htxt_name);

% set same neuron flag - if true, use the same neurons for all all textures
% not choosing new neuron for each texture
same_neuron_flag = 1;
min_area_flag = 1; %if true, uses 'and' scanning area for textures. if false, uses 'or' area, where any scanned area is valid for neurons.
cross_fold = 1;
plot_flag = 0;
% choose good neurons
% good_neurons = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 22 25 28 33 34];
good_neurons = 1:39;

num_textures = length(filename_gel);

% set hyper params
[scatter_size, t_if_mean_f_if_min, ppm, top_neuron_number, ...,
    amplitude, aff_density, speed, pin_radius, gel_weight, gel_num, texture_type,...
    len, stopBand, time_samp_period] = set_hyperparams(ppms(m), rads(m));

% pull activities and real rates for each texture (for best neuron comparison)
neuron_identities = {iPC, iRA, iSA};
excludeNeurons = 0; %don't average neurons that don't fire
[activities, av_spike_trains, space_vec, raw_rates] = pullRealActivities(rates, spikes, ...
    my_texture_names, good_neurons, neuron_identities, texture_nums, ...
    speed, excludeNeurons, time_samp_period);

% set selection modes and initialize aff pop
neuron_selection_modes = set_modes();

distances_real_gel = {}; %texture_num x aff identity x [sim neurons x real neurons]
distances_real_ts = {};

[aff_pop] = make_aff_pop(aff_density);

gel_FRs = cell(3,1);
ts_FRs = cell(3,1);
locs = [];
res_cols = {};
% run loop

tic
for i = 1:num_textures
    %load
    disp(strcat("Loading data from ", filename_gel(i)));
    cd(strcat(local_data_path_str, "sim_data"));
    load(filename_gel(i), "gel");
    load(filename_nogel(i), "no_gel");
    cd(strcat("..", local_path_back));
    
    %truncate low values
    %     bottom_flag = 0;
    %     sd = 2;
    %     gel = truncateProfile(gel, sd, bottom_flag);
    %     no_gel = truncateProfile(no_gel, sd, bottom_flag);
    
    disp(strcat("Highpass filter at ", num2str(stopBand), " per mm."));
    gel = removeLowFreq(gel, stopBand, 'charles');
    no_gel = removeLowFreq(no_gel, stopBand, 'charles');
    
    disp(strcat("Calculating neural response."));
    texture_rates = {raw_rates{i,1}, raw_rates{i,2}, raw_rates{i,3}}; %cell array with three entires - vectors of rates for
    %PCs, RAs, and SAs, for this texture
    [FRs_ts, FRs_gel, loc, r, len_scan] = pullResponses(aff_pop, gel, ...
        no_gel, ppms(m), top_neuron_number, ...
        amplitude, len, speed, rads(m), ...
        texture_rates, neuron_selection_modes, flip_flag, same_neuron_flag, ...
        figure_dir);
    res_colls{i} = r;
    
    locs = [locs; loc];
    
    for p = 1:3
        aff_resps = FRs_gel{p};
        gel_FRs{p} = [gel_FRs{p}; aff_resps'];
        aff_resps = FRs_ts{p};
        ts_FRs{p} = [ts_FRs{p}; aff_resps'];
    end
    
    
    mean_ts = FRs_ts{4}';
    sem_ts = FRs_ts{5}';
    mean_gel = FRs_gel{4}';
    sem_gel = FRs_gel{5}';
    activities.ts(i,:) = [mean_ts, sem_ts];
    activities.gel(i,:) = [mean_gel, sem_gel];
    close all
    
    if plot_flag
        h3 = findobj('Type','figure');
        fig = figure(h3(3));
        fig.Position = [100 100 900 700];
        subplot(2,3,3); % plot real spikes using TouchSim plot_spikes function wrapper
        tsPlotSpikes(spikes, len_scan, good_neurons, htxt_name, neuron_identities, texture_nums(i), speed)
        title(strcat(no_gel.name, " Recorded Data"));
        ylabel("");
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        
        subplot(2,3, 6); plotFiringRates(rates, good_neurons, htxt_name, neuron_identities, texture_nums(i), speed, 1);
        %     title(strcat(no_gel.name, "Mean Rate Recorded Data"));
        title("")
        ylabel("")
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        close(h3(1), h3(2), h3(3), h3(4)); %close extraneous figures
        %
    end
    
    disp("Calculating distance metric between real and sim spike trains...")
    % spike distance
    distances_real_gel{i} = spike_dist_touchsim(spikes, r{2}.responses,... % double check r{2} is gel?
        len_scan, good_neurons, neuron_identities, texture_nums(i), speed);
    % each entry is, for each aff,num sim affs by num real affs
    distances_real_ts{i} = spike_dist_touchsim(spikes, r{1}.responses,...
        len_scan, good_neurons, neuron_identities, texture_nums(i), speed);
end

total_time = toc;
tpt(m) = total_time/num_textures;
disp(strcat("average time per texture: ", num2str(tpt(m))));


%% find affs that correlate well
aff_names = ["PCs", "RAs", "SAs"];
if same_neuron_flag %if we allow the same neurons to be chosen multiple times
    if min_area_flag % 
        %restrict affs to those in common scanning window for RA and SA.
        loc_mins = min(abs(locs)); loc_mins(1) = -loc_mins(1); loc_mins(3) = -loc_mins(3);
    else
        %restrict affs to those largest communal scanning window for RA and SA.
        loc_mins = max(abs(locs)); loc_mins(1) = -loc_mins(1); loc_mins(3) = -loc_mins(3);
    end
    FRs = {gel_FRs, ts_FRs};
    in_loc_indices = cell(2,3); % stores indices in lists for RA and SA affs that are within scan locations (and gel and ts)
    in_loc_indices{1,1} = 1:length(gel_FRs{1,1}); %all PC indices are "in location" bc PCs have large receptive fields
    in_loc_indices{2,1} = 1:length(gel_FRs{1,1});
    for i = 1:2 %for gel and touchsim
        for j = 2:3 % for RAs and SAs, smaller RF afferents for which scanning
            %area matters
            % 1's where those affs are in the scan area
            neuron_indices_mask = neuron_in_loc(aff_pop, loc_mins, j);
            in_loc_indices{i,j} = find(neuron_indices_mask);
            rates_for_aff_for_texture = FRs{i}{j};
            rates_for_aff_for_texture = rates_for_aff_for_texture(:, neuron_indices_mask);
            FRs{i}{j} = rates_for_aff_for_texture;
        end
    end
    % fix rates
    real_rates = cell(1,3);
    for i = 1:3
        for j = 1:num_textures %compiling real recorded rates together so that
            %for each aff class, you have a texture by aff number matrix
            %(entries are rates) to correlate with sim afferents
            
            real_rates{i} = [real_rates{i}; raw_rates{j, i}];
        end
    end
    % to do no cross fold, set cross fold var to 1
    [distances_real_gel, distances_real_ts, ...
        activities] = cross_fold_wrapper(real_rates, ...
        FRs, in_loc_indices, aff_pop, aff_names, num_textures, ...
        res_colls, spikes, len_scan, good_neurons,  htxt_name, ...
        neuron_identities, texture_nums, speed, no_gel, rates, activities, ...
        cross_fold);
end



%% Mother of all plots
%
if isstring(amplitude)
    amplitude = 0;
end
c= date;

%activity params: gel weight, gel_num, top_neuron_number,
%touchsim_amplitude, ppm, speed, aff_density, modes (PC,RA,SA)
%NOTE: touchsim amplitude is written as "0" if maximum. gel num is written
%as 0 if gel num varies.

% title_str = strcat(texture_type, "_activities_", c, "_",...
%     num2str(gel_weight), "_", ...
%     num2str(gel_num), "_", ...
%     num2str(top_neuron_number), "_", ...
%     num2str(amplitude), "_",...
%     num2str(ppms(m)), "_", ...
%     num2str(speed), "_",...
%     num2str(aff_density), ...
%     neuron_selection_modes(1), ...
%     neuron_selection_modes(2), ...
%     neuron_selection_modes(3), ".mat");

% [rmses, rs, mean_ratios] = motherOfAllPlotsFunc(activities);
[rmses, rs, mean_ratios] = MOAPnew(activities);
root_means{m} = rmses;
pearsons{m} = rs;
rati{m} = mean_ratios;
% sgtitle(title_str, 'Interpreter', 'none') ;

%% Bar plot average spike distance

spike_distances = zeros(length(filename_gel), 6);
min_distances = {};
sd_spike_distances = spike_distances;
% first two cols, PC mean gel and PC mean no gel. second two, RA mean gel,
% RA mean no gel. last two, SA mean gel, SA mean no gel.

for i = 1:num_textures % for each texture
    gel_aff_distances = distances_real_gel{i}; % get distance matrix across afferents for this texture
    ts_aff_distances = distances_real_ts{i};
    for j = 1:3
        % find means for each texture for each aff class for gel and no gel
        gel_dis = gel_aff_distances{j}; %sim by real afferents
        ts_dis = ts_aff_distances{j};
        if t_if_mean_f_if_min
            spike_distances(i, j*2-1) = mean(gel_dis(:));
            spike_distances(i, j*2) = mean(ts_dis(:));
            sd_spike_distances(i, j*2-1) = std(gel_dis(:))/sqrt(length(gel_dis(:)));
            sd_spike_distances(i, j*2) = std(ts_dis(:))/sqrt(length(ts_dis(:)));
        else
            min_distances{i, j*2-1} = min(gel_dis, [], 2); % we want MIN for every REAL afferent
            min_distances{i, j*2} = min(ts_dis, [], 2);
        end
        
    end
end

if ~t_if_mean_f_if_min
    aff_class_means = zeros(1, 6);
    aff_class_sd = aff_class_means;
    SAs = [];
    RAs = [];
    PCs = [];
    affs = {SAs, RAs, PCs};
    for i = 1:num_textures
        for j = 1:3
            affs{j} = vertcat(affs{j}, [min_distances{i, j*2-1}, min_distances{i, j*2}]);
        end
    end
    for j = 1:3
        aff_class_means((j*2-1):(j*2)) = mean(affs{j}, 1);
        aff_class_sds((j*2-1):(j*2)) = std(affs{j}, 1);
    end
else
    aff_class_means = mean(spike_distances, 1);
    aff_class_sds = std(spike_distances, 1);
end


figure;

sig_mat_textures = zeros(3, num_textures);
for i = 1:3 %plotting either scatter or bar
    
    if t_if_mean_f_if_min
        subplot(3,1,i);
    else
        subplot(1,3,i);
    end
    
    if t_if_mean_f_if_min % if taking the mean
        significance = zeros(1,num_textures);
        for j = 1:num_textures
            gel_aff_distances = distances_real_gel{j}; distance_mat_gel = gel_aff_distances{i};
            ts_aff_distances = distances_real_ts{j}; distance_mat_ts = ts_aff_distances{i};
            dist_mat_gel = distance_mat_gel(:); dist_mat_ts = distance_mat_ts(:);
            if ne(length(dist_mat_gel), length(dist_mat_ts))
                min_length = min(length(dist_mat_gel), length(dist_mat_ts));
                dist_mat_gel = dist_mat_gel(1:min_length);
                dist_mat_ts = dist_mat_ts(1:min_length);
            end
            [~, sig] = ttest(dist_mat_gel, dist_mat_ts);
            significance(j) = sig;
        end
        sig_mat_textures(i, :) = significance;
        
        distance_means = spike_distances(:, (i*2-1):(i*2));
        distance_sds = sd_spike_distances(:, (i*2-1):(i*2));
        b = bar(distance_means); % gel and ts for each afferent
        b(1).FaceColor = 'c'; b(2).FaceColor = 'r';
        title(strcat(aff_names(i), " Spike Distances"));
        xticklabels(my_texture_names);
        ylabel("Mean Spike Distance")
        hold on
        errorbars_group(distance_means, distance_sds, significance);
        if i == 1
            legend(["Gel", "No Gel"])
        end
    else %plotting scatter of best distances
        hold on
        for j = 1:num_textures
            scatter(min_distances{j, i*2}, min_distances{j, i*2-1}, scatter_size, colorscheme(j, :), 'filled'); % touch_sim x, gel y
        end
        if i == 1
            strs = my_texture_names';
            colors = colorscheme(1:size(strs,2), :);
            leg = legend([color_legend(strs', colors)]);
            leg.Box = 0;
        end
        xlabel("No Gel Spike Distance")
        ylabel("Gel Spike Distance")
        xlim([0 1]);
        ylim([0 1]);
        plot([0, 1], [0, 1], 'k');
        title(strcat(aff_names(i), " Spike Distances to Real Data"));
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
    end
end

% plotting across all textures and just for each class.
figure;
sig_vec_aff = zeros(1,3);
aff_mean_distances = zeros(3,2);
aff_sd_distances = zeros(3,2);
for i = 1:3
    aff_mean_distances(i, :) = aff_class_means((i*2-1):(i*2));
    aff_sd_distances(i, :) = aff_class_sds((i*2-1):(i*2));
    distance_means = spike_distances(:, (i*2-1):(i*2));
    [~, sig] = ttest(distance_means(:, 1), distance_means(:, 2));
    sig_vec_aff(i) = sig;
end
b = bar(aff_mean_distances);
b(1).FaceColor = colorscheme(3, :); b(2).FaceColor = colorscheme(4,:);
title("Mean Spike Distances for each afferent class")
xticklabels(aff_names);
ylabel("Spike Distance (+/- SEM)")
hold on
errorbars_group(aff_mean_distances, aff_sd_distances, sig_vec_aff);
strs = {'Gel', 'Raw'}';
colors = colorscheme(3:4, :);
leg = color_legend(strs, colors);
leg = legend(leg);
leg.Box = 0;
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
%% functions

function [local_data_path_str, local_path_back, colorscheme, figure_dir, ...
    aff_colors, neuron_id, rates, spikes, iPC, iRA, iSA, htxt_name] = setup_vars()
local_data_path_str = "../../../mwe_data/";
local_path_back = "/../bensmaia_gelsight_scripts/mwe_charles_4_20_21/simulating_responses";
addpath('helper_functions')
cd helper_functions/touchsim_gelsight
setup_path;
cd ../..

load('colorscheme')

figure_dir = 0; % do not save figures

PC_COLOR =  [255 127 0]/255;
RA_COLOR =  [30 120 180]/255;
SA_COLOR = [50 160 40]/255;

aff_colors = {PC_COLOR, RA_COLOR, SA_COLOR};

cd(strcat(local_data_path_str, "neural_data"))
load("RawPAFData")
load("TextureNames")
cd(strcat("..", local_path_back));

end



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
        filename_gel = ["210729_snowflake_fuzzyside_gel_B5_3_processed", ...
            "210729_microsuede_gel_B5_3_processed",...
            "210730_chiffon_gel_B5_3_processed"];
        filename_nogel = [ "210729_snowflake_fuzzyside_no_gel_processed", ...
            "210729_microsuede_no_gel_processed",...
            "210730_chiffon_no_gel_processed"];
        texture_nums = [43, 20, 12];
    otherwise
        error("No mode match found - choose test, 100grams, compliant, noncompliant, or all")
end

texture_names = string(htxt_name);

my_texture_names = texture_names(texture_nums);
num_textures = length(filename_gel);
disp(texture_names(texture_nums));
if ~((length(texture_nums) == num_textures) && (num_textures == length(filename_nogel)))
    error("Filenames and/or Texture numbers don't match up.")
end

end

% filename_gel = ["201118_corduroy_35_gel_trimmed", ...
%     "210226_blizzard_fleece_gel_7_200_grams_processed", ...
%     "210223_1mm_grating_gel_11_processed"];


% filename_nogel = ["201118_corduroy_no_gel_trimmed", ...
%     "210226_blizzard_fleece_no_gel_processed",...
%     "201021_1mm_grating_no_gel"];

% texture_nums = [25, 45, 50];
% texture_nums = [49, 50, 55];



function [scatter_size, t_if_mean_f_if_min, ppm, top_neuron_number, ...,
    amplitude, aff_density, speed, pin_radius, gel_weight, gel_num, texture_type,...
    len, stopBand, time_samp_period] = set_hyperparams(ppm, pin_radius)
scatter_size = 50;
t_if_mean_f_if_min = 0; %if 1, mean. if 0, min. FOR SPIKE TIMEs
% ppm = 10;
top_neuron_number = [60 60 60]; %PC, RA, SA
% top_neuron_number = [30 30 30]; %PC, RA, SA
amplitude = "max"; % "max" or value - if value, add in difference between median texture value and this value!
aff_density = 1; %afferent population density
speed = 80; %mm/s
% pin_radius = 0.05;% mm
gel_weight = 200;
gel_num = 0;
texture_type = "compliant"; %compliant, noncompliant, or combined
len = "full"; %seconds. 12 mm length / 80 mm/s so no edge scan.
stopBand = 0.3; %frequencies below 0.5 are noise
time_samp_period = 0.001; %millisecond time resolution
% ramp len
end

function [aff_pop] = make_aff_pop(aff_density)
% generate population of simulated neurons
aff_pop = affpop_hand('D2d',aff_density);

plot(aff_pop)
title("aff pop")
end

function [modes] = set_modes()
% three entries for three afferent types: PC, RA, SA.
% modes: "area" - all afferents in area
%        "top"   - average of n=top_neuron_number neurons that respond

% neuron_selection_modes =
%        "best" -  average of closest n=top_neuron_number
%        "best_area" - average of closest n=top neuron number in texture
%        "area_rand" - random n affs in area (that fire) NOTE these are
%        good for complete aff correlation where you choose matching affs
%        not texture by texture but across all responses
%        "all_rand" - random n affs (that fire)
%        "best_area_match" - match on an afferent by afferent basis, same
%        as best area
%        "best_match" - match on an aff by aff basis, same as best
% neuron_selection_modes = ["best", "best_area", "best_area"]; "best_match", "best_area_match", "best_area_match"
modes = ["all_rand", "area_rand", "area_rand"];
end

function [max_inds_no_repeats] = diff_max_for_each_col(temp)
temp(isnan(temp)) = -1; % don't chose NaN correlations
[~, sort_keys] = sort(temp, 1, 'descend');
max_inds_no_repeats = zeros(1,size(temp, 2));
for i = 1:size(temp, 2)
    sorting_inds = sort_keys(:, i);
    for j = 1:length(sorting_inds)
        if ~ismember(sorting_inds(j), max_inds_no_repeats)
            max_inds_no_repeats(i) = sorting_inds(j);
            break
        else
            continue
        end
    end
    if j == length(sorting_inds)
        error("all sim neurons somehow already selected to match real neurons")
    end
end
end

function [max_grand_inds] = define_indices(real_rates, FRs, in_loc_indices, ...
    aff_pop, aff_names, train_inds, test_inds)
be_plottin = 1;
corrs = cell(2,3);
neur_indices = cell(2,3);
for m = 1:2 % for gel and touchsim
    for i = 1:3 % for every aff class
        temp = zeros(size(real_rates{i}, 2), size(FRs{m}{i}, 2)); % real by sim affs for this aff type for gel/ts
        temp_tester = zeros(size(real_rates{i}, 2), size(FRs{m}{i}, 2));
        for j = 1:size(real_rates{i}, 2) % for every real afferent
            for k = 1:size(FRs{m}{i}, 2) % for every simmed afferent
                real_responses = real_rates{i}(:,j); % rates for every texture
                res_train = real_responses(train_inds);
                sim_responses = FRs{m}{i}(:,k);
                sim_train = sim_responses(train_inds);
                    
                tempest = corrcoef(res_train, sim_train);
                temp(j,k) = tempest(1,2);
                tempest_test = corrcoef(real_responses(test_inds), sim_responses(test_inds)); 
                temp_tester(j,k) = tempest_test(1,2);
                
            end
        end
        
        if isempty(train_inds) % i.e., no cross val
            [inds] = diff_max_for_each_col(temp_tester');
            corrs{m,i} = temp_tester;
        else
            [inds] = diff_max_for_each_col(temp'); %choose a new best simulated aff for every real aff
            corrs{m,i} = temp;
        end
            
        % on average, for test indices/textures, how good are the
        % correlations (not directly maxed?)
        test_corrs  = temp_tester(1:size(temp_tester,1), inds); % chosen sim afferents based on corrs for train indices
        
        disp(strcat(" average correlation for test textures: ", ...
            num2str(nanmean(test_corrs, 'all')), ", aff_index: ", num2str(i)));
        if be_plottin
            figure;
            title(strcat("first ", aff_names(i), " (blue) firing over all textures and most corr sim aff (red)"));
            hold on
            plot(1:length(FRs{1}{i}(:,inds(1))), FRs{1}{i}(:,inds(1)), 'r')
            plot(1:length(real_rates{i}(:,1)), real_rates{i}(:,1), 'b')
            close
        end
        neur_indices{m,i}  = inds;
        
    end
end

aff_inds = {[aff_pop.afferents.iPC], [aff_pop.afferents.iRA], [aff_pop.afferents.iSA1]};
max_grand_inds = cell(1,2);

% for gel and touchsim
for m = 1:2
    inds_to_accumulate = [];
    for i = 1:3
        grand_inds = find(aff_inds{i}); % find inds for this aff among all
        grand_loc_inds = grand_inds(in_loc_indices{m,i});
        max_corr_inds = neur_indices{m,i}; % get the indices for best
        % matches relative to list of IN LOC afferents
        inds_to_accumulate = [inds_to_accumulate, grand_loc_inds(max_corr_inds)]; % now chosen top neurons
        % for this aff class relative to all afferents
        % and add them to a list.
    end
    max_grand_inds{m} = inds_to_accumulate; % these are the indices for best correlated matching neurons.
end

end



function [distances_real_gel, distances_real_ts, activities] = per_texture_loop(max_grand_inds, ...
    num_textures, res_colls, aff_pop, spikes, len_scan, good_neurons, ...
    htxt_name, neuron_identities, texture_nums, speed, no_gel, rates, activities)

reselected_rs = cell(1,num_textures);
for i = 1:num_textures
    gel_and_ts_old_rs = res_colls{i};
    gel_and_ts_old_rs = {gel_and_ts_old_rs{2}, gel_and_ts_old_rs{1}};
    % switch so gel then ts
    
    new_rs = cell(1,2);
    for m = 1:2 %reselecting all fields
        copy_aff_pop = AfferentPopulation(aff_pop.afferents);
        old_res_col = gel_and_ts_old_rs{m};
        new_res_col = ResponseCollection(copy_aff_pop, old_res_col.responses, old_res_col.stimulus);
        new_res_col.affpop.afferents = copy_aff_pop.afferents(:, max_grand_inds{m});
        new_res_col.responses = new_res_col.responses(:, max_grand_inds{m});
        new_rs{m} = new_res_col;
    end
    
    reselected_rs{i} = new_rs;
    plot_flag = 0;
    [new_FRs_ts, new_FRs_gel] = new_calc_responses(new_rs, plot_flag, ["names", "names"]);
    mean_ts = new_FRs_ts{4}';
    sem_ts = new_FRs_ts{5}';
    mean_gel = new_FRs_gel{4}';
    sem_gel = new_FRs_gel{5}';
    activities.ts(i,:) = [mean_ts, sem_ts];
    activities.gel(i,:) = [mean_gel, sem_gel];
    
    if plot_flag
        h3 = findobj('Type','figure');
        fig = figure(h3(1));
        fig.Position = [100 100 900 700];
        subplot(2,3,3); % plot real spikes using TouchSim plot_spikes function wrapper
        tsPlotSpikes(spikes, len_scan, good_neurons, htxt_name, neuron_identities, texture_nums(i), speed)
        title(strcat(no_gel.name, " Recorded Data"));
        ylabel("");
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';

        subplot(2,3, 6); plotFiringRates(rates, good_neurons, htxt_name, neuron_identities, texture_nums(i), speed, 1);
        %     title(strcat(no_gel.name, "Mean Rate Recorded Data"));
        title("")
        ylabel("")
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
    end
    %         close(h3(1)); %close extraneous figures
    %
    
    disp("Calculating distance metric between real and sim spike trains...")
    % spike distance
    distances_real_gel{i} = spike_dist_touchsim(spikes, new_rs{2}.responses,... % double check r{2} is gel?
        len_scan, good_neurons, neuron_identities, texture_nums(i), speed);
    % each entry is, for each aff,num sim affs by num real affs
    distances_real_ts{i} = spike_dist_touchsim(spikes, new_rs{1}.responses,...
        len_scan, good_neurons, neuron_identities, texture_nums(i), speed);
end
end


function [distances_real_gel, distances_real_ts, ...
    activities_out] = cross_fold_wrapper(real_rates, ...
    FRs, in_loc_indices, aff_pop, aff_names, num_textures, ...
    res_colls, spikes, len_scan, good_neurons,  htxt_name, ...
    neuron_identities, texture_nums, speed, no_gel, rates, activities, ...
    cross_fold)


cross_val_div = round(num_textures/cross_fold);
cross_val_devs = 1:cross_val_div:num_textures; 
cross_val_devs(end+1) = num_textures+1; % for each run thru the loop we
% are correlating for everything BUT the textures within these indices so
% we can use them as a test set.

activities_out = activities;

for x = 1:cross_fold % for every cross fold batch
    
    if x == 1 % if the first
        train_inds = cross_val_devs(2):(cross_val_devs(end)-1); % take every batch but the first for training
    elseif x == cross_fold % if the last
        train_inds = 1:(cross_val_devs(end-1)-1); %every batch but the last for training
    else
        before_inds = 1:(cross_val_devs(x)-1); %inds before the test batch
        after_inds = cross_val_devs(x+1):(cross_val_devs(end)-1); %inds after test_batch
        train_inds = [before_inds, after_inds];
    end
    
    test_inds = cross_val_devs(x):cross_val_devs(x+1)-1;
    
    [max_grand_inds] = define_indices(real_rates, FRs, in_loc_indices, aff_pop, aff_names, train_inds, test_inds); 
    % these are inds of best correlating neurons but ONLY based on train
    % textures - not test set - cross_val_dev(x):cross_val_dev(x+1)-1
    
    new_num_textures = length(test_inds);
    new_texture_nums = texture_nums(test_inds);
    new_activities = activities;
    new_activities.names = activities.names(test_inds);
    new_activities.real = activities.real(test_inds, :);
    new_activities.ts = activities.ts(test_inds, :);
    new_activities.gel = activities.gel(test_inds, :);

    [partial_distances_real_gel, partial_distances_real_ts, partial_activities] = per_texture_loop(max_grand_inds, ...
        new_num_textures, res_colls, aff_pop, spikes, len_scan, good_neurons, ...
        htxt_name, neuron_identities, new_texture_nums, speed, no_gel, rates, new_activities);
    close all;
    activities_out.real(test_inds, :) = partial_activities.real;
    activities_out.gel(test_inds, :) = partial_activities.gel;
    activities_out.ts(test_inds, :) = partial_activities.ts;
    distances_real_gel(test_inds) = partial_distances_real_gel; 
    distances_real_ts(test_inds) = partial_distances_real_ts;
end

end





% textures to do
% - satin
% organza
% - chiffon
% taffeta
% - flag banner
% crinkled silk
% designer wool
% - drapery tape (foam side)
% corrugated paper
% wool felt
% metallic silk
% silk jacquard
% wool crepe