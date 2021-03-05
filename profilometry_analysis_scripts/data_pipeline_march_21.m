%% Main Script for Analysis
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry
% Date: March 2 2020
% cd ~/Documents/bensmaia_lab/bensmaia_gelsight_scripts/profilometry_analysis_scripts
clear 
close all

load("activities_max_indentation");

activities.gel = [];
activities.ts = [];
%% set vars
% filename_gel = ["210219_sueded_cuddle_gel_7_processed", "210217_wool_blend_gel_7_processed", ...
%     "210209_hucktowel_gel_11_processed", "210122_velvet_gel_4_processed", ...
%     "201118_corduroy_35_gel_trimmed", "210219_cross_gel_7_processed", ...
%     "201116_2mm_grating_35_gel_processed", "201116_1mm_grating_35_gel_processed", ...
%     "210216_3mm_grating_gel_7_processed", "210216_3mm_grating_gel_11_processed" ];
% filename_nogel = ["210222_sueded_cuddle_no_gel_processed", "210216_wool_blend_no_gel_processed", ...
%     "210204_hucktowel_nogel_processed", "210121_velvet_no_gel_processed", ...
%     "201118_corduroy_no_gel_trimmed", "201119_cross_no_gel_processed", ...
%     "201019_no_gel_2mm_grating", "201021_1mm_grating_no_gel", ...
%     "210212_3_mm_grating_no_gel_processed", ];

filename_gel = [ "210217_wool_blend_gel_7_processed", ...
    "210223_velvet_gel_7_processed", ...
    "210209_hucktowel_gel_11_processed",  ...
    "210219_sueded_cuddle_gel_7_processed",...
    "201118_corduroy_35_gel_trimmed", ...
    "210226_blizzard_fleece_gel_7_200_grams_processed",...
    "210223_1mm_grating_gel_11_processed",...
    "210216_3mm_grating_gel_7_processed"];
   
filename_nogel = [ "210216_wool_blend_no_gel_processed", ...
    "210121_velvet_no_gel_processed", ...
    "210204_hucktowel_nogel_processed",  ...
    "210222_sueded_cuddle_no_gel_processed",...
    "201118_corduroy_no_gel_trimmed", ...
    "210226_blizzard_fleece_no_gel_processed",...
    "201021_1mm_grating_no_gel",...
    "210212_3_mm_grating_no_gel_processed"];

% figure_dir = "../../pngs/feb_23_charles_checkin/sueded_cuddle";
figure_dir = "E:\Users\somlab\Documents\bensmaia_lab\pngs\march_3_check_in\blizzard_fleece"; % do not save figures
ppm = 7;
ts_amplitude = "max";

% for i = 1:length(filename_gel)
%     processAndUpdate(filename_gel(i), 1);
%     close all
%     processAndUpdate(filename_nogel(i), 0);
%     close all
% end

[FRs_ts, FRs_gel, r, a] = pullResponses("210226_blizzard_fleece_gel_7_200_grams_processed", "210226_blizzard_fleece_no_gel_processed", ppm, ts_amplitude, 0);


for i = 1:length(filename_gel)
    [FRs_ts, FRs_gel, r, a] = pullResponses(filename_gel(i), filename_nogel(i), ppm, ts_amplitude, figure_dir);
    mean_ts = FRs_ts{4}';
    sem_ts = FRs_ts{5}';
    mean_gel = FRs_gel{4}';
    sem_gel = FRs_gel{5}';
    activities.ts(i,:) = [mean_ts, sem_ts];
    activities.gel(i,:) = [mean_gel, sem_gel];
    close all
end



%% Texture options

% COMPLIANT
% sueded cuddle
filename_gel = "210219_sueded_cuddle_gel_7_processed";
filename_nogel = "210222_sueded_cuddle_no_gel_processed";

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

% 2/22/21 1mm grating
% filename_gel = 
% filename_nogel = "210222_1mm_grating_no_gel_processed";
