%% AlignGelImages
close all
clear
%% Load data
cd ../../mat_files
filename_gel = "210118_upholstry2_gel_3_processed";
filename_nogel = "210119_upholstry_2_no_gel_processed";

load(filename_gel);
load(filename_nogel);
% gel = gel_dot_200925;
% no_gel = no_gel_dot_200925;
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
before = gel.profile;

% gel.profile(gel.profile == 0) = nan;
% no_gel.profile(no_gel.profile == 0) = nan;
figure;
visualizeProfile(no_gel)
figure;
visualizeProfile(gel)

%% Attempt to Align images
%Nans throw an error when aligning images. So no nans?

[optimizer, metric] = imregconfig('monomodal');
% reducing optimizer precision for better run time
optimizer.MinimumStepLength = 2e-4;
optimizer.MaximumStepLength = 0.07;
movingRegistered = imregister(gel.profile, no_gel.profile, 'affine', ...
    optimizer, metric);
%we want translation and rotation, so 'rigid'.
figure;
title("Difference Between Gel and No Gel Before and After Alignment")
subplot(1,2,1)
imshowpair(no_gel.profile, gel.profile,'Scaling', 'joint')
subplot(1,2,2)
imshowpair(no_gel.profile, movingRegistered,'Scaling','joint')

figure;
title("Changes to gel.profile")
imshowpair(gel.profile, movingRegistered, 'Scaling', 'joint')



gel.profile = movingRegistered;
figure;
subplot(1,2,1);
imagesc(before);
title("Before")
subplot(1,2,2);
imagesc(gel.profile);
title("After");
filename = strcat(filename_gel, "_aligned");
cd mat_files
save(filename, "gel");
cd .. 
