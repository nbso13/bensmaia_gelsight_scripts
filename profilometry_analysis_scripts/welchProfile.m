function [pxx, f] = welchProfile(profile_struct)
%welchProfile computes the normalized power spectrum according to welch's method for the given profile. 
profile_struct = rotateProfilometry(profile_struct, 90);
fs = 1/profile_struct.y_res;
[pxx, f] = pwelch(profile_struct.profile, [], [], [], fs);
pxx = mean(pxx, 2);
end

