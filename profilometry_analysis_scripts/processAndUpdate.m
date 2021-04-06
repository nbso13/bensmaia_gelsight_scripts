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



%% crop
crop = 1;

while crop
    figure;
    visualizeProfile(prof);
    title("Profile")
    str = input("crop profile? (y/n)", 's');
    
    if str == "y"
        direc = input("direction? ('top', 'bottom', 'left', 'right')", 's');
        unit = input("unit? ('mm', 'px')", 's');
        amount = input("amount? (number)", 's');
        
        prof = cropProfile(prof, direc, str2double(amount), unit);
        figure
        visualizeProfile(prof);
    elseif str == 'n'
        crop = 0;
    end
end


%% remove trend
figure
visualizeProfile(prof);
str = input("detrend profile? (y/n)", 's');
if str == "y"
    new_prof = removeTrend(prof);
    new_prof.detrended = 1;
    figure
    visualizeProfile(new_prof);
    str = input("undo detrend? (y/n)", 's');
        if str == "n"
            prof = new_prof;
        end
end

% take out minimum
prof.profile = prof.profile - min(min(prof.profile));
figure
visualizeProfile(prof);

%% save
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
