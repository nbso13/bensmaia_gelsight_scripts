close all
clear

gel_constant = 1.48; %empirically derived factor to scale profilometry through gel
    


%% Load data and process data

filename_gel = "210217_wool_blend_gel_7_processed";
disp(strcat("Loading data from ", filename_gel));

cd ../../mat_files/
load(filename_gel, "gel");
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts


% freq = 0.12; %2/mm
% amp = 1500;
% res = 0.8;
% window_size = 5000;
% gel = generate_texture("grating", freq, amp, 10, res, window_size);


[pxx_gel, f_gel] = welchProfile(gel);
figure; 
subplot(2,2,1)
plotWelch(pxx_gel, f_gel);
title("Before filtering")

subplot(2,2,2)
visualizeProfile(gel);
title("Gel")


stopBand =0.5; %anything under 0.5/mm ie happens once every two mm. freq is too low.
gel_new = removeLowFreq(gel, stopBand);

[pxx_gel_new, f_gel_new] = welchProfile(gel_new);
subplot(2,2,3)
plotWelch(pxx_gel_new, f_gel_new);
title("After filtering")

subplot(2,2,4)
visualizeProfile(gel_new);
title("Gel")
