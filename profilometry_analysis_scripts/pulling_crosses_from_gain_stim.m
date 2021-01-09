%% pulling crosses from gain stim 1/08/21 for embossed cross show
close all
clear
filename_gel = "201119_gain_gel_processed";
filename_nogel = "201119_gain_no_gel_processed";

cd ../../mat_files/
load(filename_gel);
load(filename_nogel);
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts

%gel.profile = gel.profile.*gel_constant; %scale up

if ~checkSizeMatch(gel, no_gel)
    [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
    [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
end

gel = rotateProfilometry(gel, 90);
no_gel = rotateProfilometry(no_gel, 90);

figure
visualizeProfile(gel);
figure
visualizeProfile(no_gel);
gel = cropProfile(gel, 'left', 10, 'mm');
no_gel = cropProfile(no_gel, 'left', 10, 'mm');
no_gel = cropProfile(no_gel, 'bottom', 5, 'mm');
gel = cropProfile(gel, 'bottom', 5, 'mm');
gel = cropProfile(cropProfile(gel, 'left', 0.65, 'mm'), 'bottom', 0.7, 'mm');
gel = cropProfile(gel, 'bottom', 0.35, 'mm');
no_gel = cropProfile(cropProfile(no_gel, 'left', 0.65, 'mm'), 'bottom', 0.9, 'mm');
figure
visualizeProfile(gel);
figure
visualizeProfile(no_gel);
cd ../../mat_files/
save('201119_cross_gel_processed', 'gel')
save('201119_cross_no_gel_processed', 'no_gel')