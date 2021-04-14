function [FRs_ts, FRs_gel, response_collection, aff_pop_final] = pullResponses(gel, ...
    no_gel, ppm, top_neuron_number, amplitude, len, speed, pin_radius, aff_density, ...
     texture_rates, neuron_selection_modes, figure_dir)
%pullResponses: given struct filenames and other hyperparams, calc firing
%rates. filenames indicate mat file name. ppm is pins per millimeter for
%touchsim model. ts amplitude indicates how much of the texture to input to
%skin mechanics for touchsim (must use entire gelsight profile). if figure
%dir is a string, save figures there.

% presets
samp_freq = 2000; % hz
ramp_len = 0.01; %length of ramping on, in seconds

%% Load data and process data

save_figures = 0;
if isstring(figure_dir)
    dir_list = split(figure_dir, "/");
    texture_name = dir_list(end);
    save_figures = 1;
    disp(strcat("Saving figures to ", figure_dir));
end


% freq = 0.1; %2/mm
% amp = 1500;
% res = 10;
% window_size = 5000;
% gel = generate_texture("grating", freq, amp, 10, res, window_size);
% disp("Adjusting crop and sample rate...")
% if ~checkSizeMatch(gel, no_gel)
%     [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
%     [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
% end

% calculate length of time of scan

if isstring(len) %then we do two full scans through
    len = 2*gel.x_axis(end)/speed;
end

%% build models
% str = input("View touchsim surfaces? (y/n)", 's');
str = 'n';
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

skin_surface_ts.offset = skin_surface_ts.offset- min(skin_surface_ts.offset);
new_gel_ts.offset = new_gel_ts.offset- min(new_gel_ts.offset);
gel_amplitude = 0;
ts_amplitude = 0;
 %figure out amplitude
if amplitude ~= "max"
     %if it is max, just completely indent both profiles from min to max, but no added amp.
    if (amplitude > median(new_gel_ts.offset)) && (amplitude > median(skin_surface_ts.offset))
        disp("amplitudes both low, adjusting to given value")
        gel_amplitude = amplitude - median(new_gel_ts.offset);
        ts_amplitude = amplitude - median(skin_surface_ts.offset);
    end
end

% skin_surface_ts.amp = 1.95; %1/19 dots
skin_surface_ts.offset = skin_surface_ts.offset + ts_amplitude;
new_gel_ts.offset = new_gel_ts.offset + gel_amplitude;

% calculate welch's method
skin_surface_profile = shape2profilometry(skin_surface_ts.shape, skin_surface_ts.offset, ppm);
[pxx_ts, f_ts] = welchProfile(skin_surface_profile);

%% calc_responses
% str = input("Calculating neural responses. Display figures? (y/n)", 's');
str = 'n';
if str == "y"
    plot_flag = 1;
else
    plot_flag = 0;
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
end

