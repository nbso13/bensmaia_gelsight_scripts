function [] = profSpectralAnalysis(filename_gel, ...
    filename_nogel, stopBand, figure_dir)
%pullResponses: given struct filenames and other hyperparams, calc firing
%rates. filenames indicate mat file name. ppm is pins per millimeter for
%touchsim model. ts amplitude indicates how much of the texture to input to
%skin mechanics for touchsim (must use entire gelsight profile). if figure
%dir is a string, save figures there.
gel_constant = 1.48; %empirically derived factor to scale profilometry through gel


%% Load data and process data

save_figures = 0;
if isstring(figure_dir)
    dir_list = split(figure_dir, "/");
    texture_name = dir_list(end);
    save_figures = 1;
    disp(strcat("Saving figures to ", figure_dir));
end

disp(strcat("Loading data from ", filename_gel));

cd ../../mat_files/
load(filename_gel, "gel");
load(filename_nogel, "no_gel");
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts

if ~isfield(gel, 'scaled')
    disp(strcat(filename_gel, " apparently not scaled - scaling up by gel constant."));
    gel.profile = gel.profile.*gel_constant; %scale up
    
end

% freq = 0.1; %2/mm
% amp = 1500;
% res = 10;
% window_size = 5000;
% gel = generate_texture("grating", freq, amp, 10, res, window_size);


% disp("Adjusting crop and sample rate...")
% if ~checkSizeMatch(gel, no_gel)
%         [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
%     [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
% end


figure;
subplot(4,3,1)
visualizeProfile(gel);
title("Gel Before filtering")
subplot(4,3,2)
visualizeProfile(no_gel);
title("No Gel before filtering")

[pxx_gel, f_gel] = welchProfile(gel);
[pxx_no_gel, f_no_gel] = welchProfile(no_gel);

subplot(4, 3, 4)
plotWelch(pxx_gel, f_gel)
title("Gel Before filtering")

subplot(4, 3, 5)
plotWelch(pxx_no_gel, f_no_gel)
title("No Gel Before filtering")


subplot(4,3,6)
interp_gel = interp1(f_gel, pxx_gel, f_no_gel);
plotWelch(interp_gel./pxx_no_gel, f_no_gel)
title("Gel : No Gel ratio")
ylabel("Ratio")



gel = removeLowFreq(gel, stopBand);

no_gel = removeLowFreq(no_gel, stopBand);

subplot(4,3,7)
visualizeProfile(gel);
title("Gel after filtering")
subplot(4,3,8)
visualizeProfile(no_gel);
title("No Gel after filtering")

[pxx_gel, f_gel] = welchProfile(gel);
[pxx_no_gel, f_no_gel] = welchProfile(no_gel);

subplot(4, 3, 10)
plotWelch(pxx_gel, f_gel)
title("Gel after filtering")

subplot(4, 3, 11)
plotWelch(pxx_no_gel, f_no_gel)
title("No Gel after filtering")


subplot(4,3,12)
interp_gel = interp1(f_gel, pxx_gel, f_no_gel);
plotWelch(interp_gel./pxx_no_gel, f_no_gel)
title("Gel : No Gel ratio")
ylabel("Ratio")

sgtitle(strcat(no_gel.name, " Spectral Analysis"));
end