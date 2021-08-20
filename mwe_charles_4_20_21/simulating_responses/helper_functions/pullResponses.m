function [FRs_ts, FRs_gel, loc, response_collection, len_scan] = pullResponses(aff_pop_in, gel, ...
    no_gel, ppm, top_neuron_number, amplitudes, len, speed, pin_radius, ...
     texture_rates, neuron_selection_modes, flip_flag, figure_dir)
%pullResponses: given struct filenames and other hyperparams, calc firing
%rates. filenames indicate mat file name. ppm is pins per millimeter for
%touchsim model. ts amplitude indicates how much of the texture to input to
%skin mechanics for touchsim (must use entire gelsight profile). if figure
%dir is a string, save figures there.

% presets
samp_freq = 500; % hz

% calculate length of time of scan

if isstring(len) %then we do two full scans through
    len = 2*gel.x_axis(end)/speed;
end
len_scan = len;

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
    surf_figures] = TouchSimSkin(gel, no_gel, ppm, pin_radius, flip_flag, plot_flag);
    gcf;
    sgtitle("Gel, Texture, and TouchSim Profiles");
else
    [new_gel_ts, new_no_gel_ts, skin_surface_ts] = TouchSimSkin(gel, no_gel, ppm, pin_radius, flip_flag, plot_flag);
end

    
%% set up touchsim

skin_surface_ts.offset = skin_surface_ts.offset- min(skin_surface_ts.offset);
new_gel_ts.offset = new_gel_ts.offset- min(new_gel_ts.offset);
gel_amplitude = 0;
ts_amplitude = 0;
 %figure out amplitude
if amplitudes ~= "max"
     %if it is max, just completely indent both profiles from min to max, but no added amp.
    if (amplitudes > median(new_gel_ts.offset)) && (amplitudes > median(skin_surface_ts.offset))
        disp("amplitudes both low, adjusting to given value")
        gel_amplitude = amplitudes - median(new_gel_ts.offset);
        ts_amplitude = amplitudes - median(skin_surface_ts.offset);
    end
end

skin_surface_ts.offset = skin_surface_ts.offset + ts_amplitude;
new_gel_ts.offset = new_gel_ts.offset + gel_amplitude;

% calculate welch's method
% skin_surface_profile = shape2profilometry(skin_surface_ts.shape, skin_surface_ts.offset, ppm);
% [pxx_ts, f_ts] = welchProfile(skin_surface_profile);

%% calc_responses
% str = input("Calculating neural responses. Display figures? (y/n)", 's');
str = 'y';
if str == "y"
    plot_flag = 1;
else
    plot_flag = 0;
end
    
[FRs_ts, FRs_gel, loc, response_collection, figure_handles] = calcResponses(aff_pop_in, skin_surface_ts,...
    new_gel_ts, ppm, speed, len, samp_freq, top_neuron_number, ...
    texture_rates, neuron_selection_modes, plot_flag);
end


