%% Touchsim Force / Amplitude curve
% Jan 18 21

% Approach: to compare gelsight and touchsim neural simulations, we need to be 
% using the same force per unit area, or pressure. Touchsim takes in profiles
% and amplitudes of indentation. If we can build a pressure vs amplitude curve
% for a stimulus, we can then choose the amplitude that matches the pressure used
% with gelsight.
% 
% steps:
% 1) calculate the area of the texture
% 2) sum the forces at a range of amplitudes
% 3) plot

clear 
close all

%% set vars
ppm = 18;
gel_constant = 1.49;

%cross
filename_gel = "201119_cross_gel_processed";
filename_nogel = "201119_cross_no_gel_processed";

% CORDUROY
% filename_gel = "201118_corduroy_35_gel_trimmed";
% filename_nogel = "201118_corduroy_no_gel_trimmed";


%% Load data process data
cd ../../mat_files/
load(filename_gel);
load(filename_nogel);
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts

gel.profile = gel.profile.*gel_constant; %scale up

if ~checkSizeMatch(gel, no_gel)
    [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
    [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
end

% gel = rotateProfilometry(gel, 90);
% no_gel = rotateProfilometry(no_gel, 90);
gel_area = gel.x_axis(end)*gel.y_axis(end);
no_gel_area = no_gel.x_axis(end)*no_gel.y_axis(end);

%% generate touchsim models

plot_flag = 0;
pin_radius = 0.025;
[new_gel_ts, new_no_gel_ts, skin_surface_ts] = TouchSimSkin(gel, no_gel, ppm, pin_radius, plot_flag);

%get profiles
touchsim_gel = shape2profilometry(skin_surface_ts.shape, ...
    skin_surface_ts.offset, skin_surface_ts.pins_per_mm);
new_gel = shape2profilometry(new_gel_ts.shape, ...
    new_gel_ts.offset, new_gel_ts.pins_per_mm);
new_no_gel = shape2profilometry(new_no_gel_ts.shape, ...
    new_no_gel_ts.offset, new_no_gel_ts.pins_per_mm);

%show the profiles
figure
visualizeProfile(touchsim_gel);
figure
visualizeProfile(new_gel);
figure
visualizeProfile(new_no_gel);


cd ../touchsim_gelsight
setup_path;
cd ../profilometry_analysis_scripts
ts_structs = [skin_surface_ts, new_gel_ts, new_no_gel_ts];
amplitudes = 0.1:0.1:1;
gel_mass = 200; %200 grams used
plot_flag = 1;
ts_struct = skin_surface_ts;
[~] = ampCurve(ts_struct, pin_radius, gel_area, gel_mass, amplitudes, plot_flag);