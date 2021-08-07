%% Main Script for Analysis
% Author: Nick Ornstein
% Group: Bensmaia Lab
% Project: Gelsight Profilometry
% Date: June 2021
clear
close all

local_data_path_str = "../../../mwe_data/";
local_path_back = "../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/simulating_responses";
addpath("helper_functions")

%% set vars
% 
filename_gel = ["201118_corduroy_35_gel_trimmed", ...
    "210226_blizzard_fleece_gel_7_200_grams_processed", ...
    "210223_1mm_grating_gel_11_processed"];
% 
filename_nogel = ["201118_corduroy_no_gel_trimmed", ...
    "210226_blizzard_fleece_no_gel_processed",...
    "201021_1mm_grating_no_gel"];

% filename_gel = ["210310_sueded_cuddle_gel_11_100_grams_processed"];
% 
% filename_nogel = ["210222_sueded_cuddle_no_gel_processed"];
    
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


subplot_dim_1 = 2;
subplot_dim_2 = 2;

ind_gel = 2;
ind_no_gel = 1;
tic
for i = 1:length(filename_gel)
    fig = figure;
    fig.Position = [100, 100, 900, 600];
    %load
    disp(strcat("Loading data from ", filename_gel(i)));
    cd(strcat(local_data_path_str, "sim_data"))
    load(filename_gel(i), "gel");
    load(filename_nogel(i), "no_gel");
    cd(local_path_back)
    
    
    fig.Name = no_gel.name;
    disp(strcat("Highpass filter at ", num2str(stopBand), " per mm."));
    
    gel = removeLowFreq(gel, stopBand, 'charles');
    no_gel = removeLowFreq(no_gel, stopBand, 'charles');
    
%     gel.profile = gel.profile + (max(no_gel.profile(:)) - max(gel.profile(:))); % indentation map
    
    plotExProf(total_texture_number, gel, no_gel, scale_bar_loc, subplot_dim_1, subplot_dim_2, ind_gel, ind_no_gel)
    ax_raw = gca;
    cbar_raw = fig.Children(1);
    cbar_gel = fig.Children(3);
    
    box = annotation("textbox", ...
        [0.246486486486487 0.0761285722121672 0.498378378378379 0.0889967637540443],...
        'String', ...
        strcat('The surface of a swath of ', lower(no_gel.name), ...
                ' comprises small compliant ', ...
                'elements (top, height range, 0-', ...
                num2str(round(100*cbar_gel.Limits(2))/100), ...
                'mm). The gel compresses those ', ...
                'elements as would the finger, revealing the corresponding pattern of ', ...
                ' deformations in the skin (bottom, height range, 0-', ...
                num2str(round(100*cbar_raw.Limits(2))/100), 'mm).'),...
        'FontWeight','bold',...
        'EdgeColor','none');
    cbar_raw.Visible = 'off'; cbar_gel.Visible = 'off';
    
end
total_time = toc;
disp(strcat("average time per texture: ", num2str(total_time/length(filename_gel))))
%%
FolderName = 'D:\Users\Somlab\Documents\Nick_Ornstein\figures';
cd(FolderName)
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  saveas(FigHandle, strcat(FigName, '.png'));
end


% cd("D:\Users\Somlab\Documents\Nick_Ornstein\figures")
% im = imread("blizzard_image.jpg");
% im = imrotate(im, 90);
% gcf;subplot(subplot_dim_1, subplot_dim_2, 3); imshow(im)
% title("Blizzard Fleece Color Image")
% 

% 
% filename_nogel = "210222_sueded_cuddle_no_gel_processed";
% filename_gel = "210217_wool_blend_gel_7_processed";
% load
% disp(strcat("Loading data from ", filename_gel));
% cd(local_data_path_str)
% load(filename_nogel, "no_gel");
% cd(local_path_back)
% 
% disp(strcat("Highpass filter at ", num2str(stopBand), " per mm."));
% 
% no_gel = removeLowFreq(no_gel, stopBand, 'charles');
% 
% visualizeProfile(no_gel)
% yticks([]); yticklabels({})
% xticks([]); xticklabels({})
% ylabel("");
% xlabel("");
% daspect([1 1 1])
