%% Main Script for Analysis
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry
% Date: June 2021
clear
close all
local_data_path_str = "../../../mwe_data/sim_data";
local_path_back = "../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/simulating_responses";
addpath("helper_functions")

%% set vars

filename_gel = ["201118_corduroy_35_gel_trimmed", ...
    "210226_blizzard_fleece_gel_7_200_grams_processed", ...
    "210223_1mm_grating_gel_11_processed"];



filename_nogel = ["201118_corduroy_no_gel_trimmed", ...
    "210226_blizzard_fleece_no_gel_processed",...
    "201021_1mm_grating_no_gel"];
    
%HYPERPARAMS
stopBand = 0.3; %frequencies below 0.5 are noise
scale_bar_loc = [1 1];
total_texture_number = length(filename_gel);

%% Preprocess data
% for i = 1:length(filename_gel)
%     gel_flag = 1;
%     processAndUpdate(filename_gel(i), gel_flag);
%     gel_flag = 0;
%     processAndUpdate(filename_nogel(i), gel_flag);
% end

%% Run Loop

fig = figure;
fig.Position = [100, 100, 900, 600];
% tic
% for i = 1:length(filename_gel)
%     %load
%     disp(strcat("Loading data from ", filename_gel(i)));
%     cd data
%     load(filename_gel(i), "gel");
%     load(filename_nogel(i), "no_gel");
%     cd ..
%     
%     disp(strcat("Highpass filter at ", num2str(stopBand), " per mm."));
%     
%     gel = removeLowFreq(gel, stopBand, 'charles');
%     no_gel = removeLowFreq(no_gel, stopBand, 'charles');
%     
%     gel.profile = gel.profile + (max(no_gel.profile(:)) - max(gel.profile(:))); % indentation map
%     
%     plotExProf(total_texture_number, i, gel, no_gel, scale_bar_loc);
%     
% end
% total_time = toc;
% disp(strcat("average time per texture: ", num2str(total_time/length(filename_gel))))
% 


filename_nogel = "210222_sueded_cuddle_no_gel_processed";
filename_gel = "210217_wool_blend_gel_7_processed";
%load
disp(strcat("Loading data from ", filename_gel));
cd(local_data_path_str)
load(filename_nogel, "no_gel");
cd(local_path_back)

disp(strcat("Highpass filter at ", num2str(stopBand), " per mm."));

no_gel = removeLowFreq(no_gel, stopBand, 'charles');

visualizeProfile(no_gel)
yticks([]); yticklabels({})
xticks([]); xticklabels({})
ylabel("");
xlabel("");
daspect([1 1 1])
