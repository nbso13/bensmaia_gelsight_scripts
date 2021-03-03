function [prof] = processAndUpdate(filename_prof, gel_flag)
%detrends, scales, rotates, and saves profilometry inputs
gel_constant = 1.48;
cd ../../mat_files/
load(filename_prof);
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts

if gel_flag
    prof = gel;
else
    prof = no_gel;
end

if gel_flag && ~isfield(prof,'scaled') % if its a gel and not yet scaled
	prof.profile = prof.profile.*gel_constant; %scale up
    prof.scaled = gel_constant;
end

figure;
visualizeProfile(prof);
title("Profile")
rotate_gel = input("Scanning left to right. Rotate profile counterclockwise? (enter #, 0 - 360)");
if ~rotate_gel == 0
    prof = rotateProfilometry(prof, rotate_gel);
end

figure
visualizeProfile(prof);
str = input("detrend profile? (y/n)", 's');

%% remove trend
if str == "y"
    prof = removeTrend(prof);
    prof.detrended = 1;
    figure
    visualizeProfile(prof);
end

str = input("save profile? (y/n)", 's');
if str == "y"
    cd ../../mat_files
    if gel_flag
        gel = prof;
        save(filename_prof, 'gel')
    else
        no_gel = prof;
        save(filename_prof, 'no_gel')
    end
    cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
end

end
