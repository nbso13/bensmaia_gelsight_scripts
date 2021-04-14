close all
clear




%% Load data and process data
% "210217_wool_blend_gel_7_processed"
% "210121_velvet_no_gel_processed"
filename_gel = "210223_velvet_gel_7_processed";
disp(strcat("Loading data from ", filename_gel));

cd ../../mat_files/
load(filename_gel, "gel");
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
prof = gel;

% freq = 0.12; %2/mm
% amp = 1500;
% res = 0.8;
% window_size = 5000;
% prof = generate_texture("grating", freq, amp, 10, res, window_size);


[pxx_prof, f_prof] = welchProfile(prof);
figure; 
subplot(2,2,1)
plotWelch(pxx_prof, f_prof);
title("Before filtering")

subplot(2,2,2)
visualizeProfile(prof);
title("Gel")


stopBand =1; %anything under 0.5/mm ie happens once every two mm. freq is too low.
prof_new = removeLowFreq(prof, stopBand, 'charles');

[pxx_prof_new, f_prof_new] = welchProfile(prof_new);
figure(1);
subplot(2,2,3)
plotWelch(pxx_prof_new, f_prof_new);
title("After filtering")

subplot(2,2,4)
visualizeProfile(prof_new);
title("Gel")
