function [FRs_ts, FRs_gel, response_collection, aff_pop_final, welch_results] = pullResponses(filename_gel, ...
    filename_nogel, ppm, stopBand, top_neuron_number, ts_amplitude, len, speed, pin_radius, aff_density, ...
     texture_rates, neuron_selection_modes, figure_dir)
%pullResponses: given struct filenames and other hyperparams, calc firing
%rates. filenames indicate mat file name. ppm is pins per millimeter for
%touchsim model. ts amplitude indicates how much of the texture to input to
%skin mechanics for touchsim (must use entire gelsight profile). if figure
%dir is a string, save figures there.

% presets
samp_freq = 2000; % hz
ramp_len = 0.01; %length of ramping on, in seconds
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


disp("Adjusting crop and sample rate...")
if ~checkSizeMatch(gel, no_gel)
%     [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
    %[gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
end


gel = removeLowFreq(gel, stopBand);
no_gel = removeLowFreq(no_gel, stopBand);

subplot(1,2,1)
visualizeProfile(gel);
title("Gel after filtering")
subplot(1,2,2)
visualizeProfile(no_gel);
title("No Gel after filtering")

[pxx_gel, f_gel] = welchProfile(gel);
[pxx_no_gel, f_no_gel] = welchProfile(no_gel);

% calculate length of time of scan

if isstring(len) %then we do two full scans through
    len = 2*gel.x_axis(end)/speed;
end

%% build models
% str = input("View touchsim surfaces? (y/n)", 's');
str = 'y';
disp("Building surface models...")
if str == "y"
    plot_flag = 1;
else
    plot_flag = 0;
end

% plot_flag = 0;
if plot_flag
[new_gel_ts, new_no_gel_ts, skin_surface_ts, ...
    surf_figures] = TouchSimSkin(gel, no_gel, ppm, pin_radius, plot_flag);
    gcf;
    sgtitle(filename_gel);
else
    [new_gel_ts, new_no_gel_ts, skin_surface_ts] = TouchSimSkin(gel, no_gel, ppm, pin_radius, plot_flag);
end

if save_figures
    sgtitle(texture_name);
    cd(figure_dir)
    temp = char(filename_gel);
    date_gel = temp(1:6);
    saveas(surf_figures, strcat("surfs_", texture_name, "_", string(date), '.png'));
    cd ../../../bensmaia_gelsight_scripts/profilometry_analysis_scripts %out of ts, hucktowel, _checkin, pngs,
end

%% set up touchsim
% density = input("Afferent population density on distal digit 2? (enter # 0.1 - 1.0)");
cd ../touchsim_gelsight
setup_path;
cd ../profilometry_analysis_scripts/

if ts_amplitude == "max"
    skin_surface_ts.amp = max(skin_surface_ts.offset); %max(skin_surface_ts.offset);
else
    skin_surface_ts.amp = ts_amplitude; %max(skin_surface_ts.offset);
end

% skin_surface_ts.amp = 1.95; %1/19 dots
new_gel_ts.amp = max(new_gel_ts.offset); %max(new_gel_ts.offset);
new_no_gel_ts.amp = 0;

% calculate welch's method
skin_surface_profile = shape2profilometry(skin_surface_ts.shape, skin_surface_ts.offset, ppm);
skin_surface_profile = rotateProfilometry(skin_surface_profile, 90);
[pxx_ts, f_ts] = welchProfile(skin_surface_profile);

%% calc_responses
% str = input("Calculating neural responses. Display figures? (y/n)", 's');
str = 'y';
if str == "y"
    plot_flag = 1;
else
    plot_flag = 0;
end

if plot_flag
    welch_fig = figure;
    interp_gel = interp1(f_gel, pxx_gel, f_no_gel);
    interp_ts = interp1(f_ts, pxx_ts, f_no_gel);
    subplot(2,3,1)
    plotWelch(pxx_no_gel, f_no_gel);
    title("No Gel")
    subplot(2,3,2)
    plotWelch(pxx_gel, f_gel);
    title("Gel")
    subplot(2,3,3)
    plotWelch(pxx_ts, f_ts);
    title("TouchSim")
    subplot(2,3,5)
    plotWelch(interp_gel./pxx_no_gel, f_no_gel)
    title("Gel : No Gel ratio")
    ylabel("Ratio")
    yticklabels('auto')
    subplot(2,3,6)
    plotWelch(interp_ts./pxx_no_gel, f_no_gel)
    ylabel("Ratio")
    yticklabels('auto')
    title("TS : No Gel ratio")
    sgtitle(strcat("Normalized Profile Power Spectra"));
    
%     fftfig = figure;
%     subplot(1,3,1)
%     calcPlot2dFFT(gel);
%     title("Gel FFT")
%     subplot(1,3,2)
%     calcPlot2dFFT(no_gel);
%     title("No Gel 2D FFT")
%     subplot(1,3, 3)
%     calcPlot2dFFT(skin_surface_profile);
%     title("TouchSim 2D FFT")
%     sgtitle(strcat("2D FFTs"));
end
    
[FRs_ts, FRs_gel, response_collection, aff_pop_final, figure_handles] = calcResponses(skin_surface_ts,...
    new_gel_ts, aff_density, ppm, speed, len, samp_freq, top_neuron_number, ...
    texture_rates, neuron_selection_modes, plot_flag);

if save_figures
    cd(figure_dir)
    direcs = ["ts", "gel"];
    dates = [string(date_no_gel), string(date_gel)];
    for i = 1:size(figure_handles, 1)
        cd(direcs(i))
        saveas(figure_handles{i,1}, strcat("stim_", direcs(i), "_", texture_name, "_", dates(i), '.png'));
        saveas(figure_handles{i,2}, strcat("response_", direcs(i), "_", texture_name, "_", dates(i), '.png'));
        cd ..
    end
%     saveas(welch_fig, strcat("welch_result", direcs(i), "_", texture_name, "_", dates(i), '.png'));
    cd ../../../bensmaia_gelsight_scripts/profilometry_analysis_scripts %out of ts, hucktowel, _checkin, pngs,
end

welch_results = {};
welch_results{1,1} = pxx_gel; welch_results{1,2} = f_gel;
welch_results{2,1} = pxx_no_gel; welch_results{2,2} = f_no_gel;
welch_results{3,1} = pxx_ts; welch_results{3,2} = f_ts;

end

