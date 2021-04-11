function [] = cropProfilePerm(filename, gel_id, profile_struct, crop_direction, pixel_value, unit)
%cropProfilePerm wrapper for cropProfile but saves result as new struct.
profile_struct = cropProfile(profile_struct, crop_direction, pixel_value, unit);
profile_struct.profile = profile_struct.profile - min(min(profile_struct.profile)); % zero out min
cd ../../mat_files
if gel_id == 1
    gel = profile_struct;
    save(filename, 'gel')
else
    no_gel = profile_struct;
    save(filename, 'no_gel')
    
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
disp(strcat("saved cropped profilometry struct as ", filename));
end

