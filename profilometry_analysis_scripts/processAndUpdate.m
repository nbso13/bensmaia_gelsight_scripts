function [prof] = processAndUpdate(filename_prof, gel_flag)
%detrends, scales, rotates, truncates, and saves profilometry inputs

cd ../../mwe_data/sim_data
load(filename_prof);
cd ../../bensmaia_gelsight_scripts/profilometry_analysis_scripts
disp(strcat("PROCESSING ", filename_prof));
if gel_flag
    prof = gel;
else
    prof = no_gel;
end


figure;
visualizeProfile(prof);
title("Profile Before Scale")

if gel_flag % if its a gel
    if ~isfield(prof, "scaled") % if not scaled
        str = input("Gel not yet scaled. Scale factor? (multiplier or gel_identity)", 's');
        if str2double(str)>3
            scale_factor = gel_id_to_factor(str);
        else
            scale_factor = str2double(str);
        end
        prof.profile = prof.profile.*scale_factor; %scale up
        prof.scaled = scale_factor;
    else %it has been scaled
        str = input(strcat("Gel already scaled at ", num2str(prof.scaled), ". Rescale? (n or scale factor, or gel number)"), 's');
        if ~(str == "n") %if they do want to rescale
            if str2double(str)>3
                scale_factor = gel_id_to_factor(str);
            else
                scale_factor = str2double(str);
            end
            prof.profile = (prof.profile./prof.scaled).*scale_factor; %unscale earlier scaler and rescale
            prof.scaled = scale_factor;
        end
    end
end

figure;
visualizeProfile(prof);
title("Profile After Scale")
rotate_gel = input("Scanning left to right. Rotate profile counterclockwise? (enter #, 0 - 360)");
% rotate_gel = 0;
if ~rotate_gel == 0
    prof = rotateProfilometry(prof, rotate_gel);
end

%% name

if ~isfield(prof, 'name')
    str = input("set name? 'n' or 'name'", 's');
    if ~strcmp(str, "n")
        disp("NAMING")
        prof.name = str;
    end
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
% str = "n";
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
prof.profile = prof.profile - min(prof.profile(:));
figure
visualizeProfile(prof);
title(prof.name);


%% Truncate

figure;
visualizeProfile(prof);
title("Profile")
str = input("Truncate Profile? y/n", 's');
old_prof = prof;
% str = "n";
if str == "y"
    truncate_prof = 1;
    while truncate_prof
        if isfield(prof, 'truncated')
            disp("Already truncated. Don't truncate a profile multiple times.")
            truncate_prof = 0;
        else
            str = input("stds to cut? (number, default 3)", 's');
            stan_dev = str2num(str);
            prof = truncateProfile(prof, stan_dev);
            prof.truncated = stan_dev; %as in, three standard devs has been truncated
            gcf;
            visualizeProfile(prof);
            title(strcat("Truncated to ", str));
            str = input("Redo truncation? (y/n)", 's');
            if str == "y"
                prof = old_prof;
            elseif str == "n"
                truncate_prof = 0;
            end
        end
    end
end

%% fill missing
if sum(isnan(prof.profile(:))) > 0
    disp("INPAINTING NANS")
    prof.profile = inpaint_nans(prof.profile);
end

prof.profile = prof.profile-min(prof.profile(:));
%% save
% str = "y";
str = input("save profile? (y/n)", 's');
if str == "y"
    cd ../../mwe_data/sim_data
    if gel_flag
        gel = prof;
        save(filename_prof, 'gel')
    else
        no_gel = prof;
        save(filename_prof, 'no_gel')
    end
    cd ../../bensmaia_gelsight_scripts/profilometry_analysis_scripts
end

end

function [factor] = gel_id_to_factor(id_str)

gel_7_factor = 1.4479;
gel_11_factor = 1.4145;
gel_18_factor = 1.4418;
gel_19_factor  = 1.4388;
if strcmp(id_str, "7")
    factor = gel_7_factor;
elseif strcmp(id_str, "11")
    factor = gel_11_factor;
elseif strcmp(id_str, "18")
    factor = gel_18_factor;
elseif strcmp(id_str, "19")
    factor = gel_19_factor;
else
    error("Gel ID not recognized.")
end
end
