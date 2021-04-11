function [pxx, f] = welchProfile(profile_struct)
%welchProfile computes the normalized power spectrum according to welch's method for the given profile. 
% - thanks to charles greenspon for code - 

% demean
profile_struct.profile = profile_struct.profile - mean(profile_struct.profile(:));

profile_struct = rotateProfilometry(profile_struct, 90);
samp_period = profile_struct.y_res;
fs = 1/samp_period; % samples per mm
% Welch's implementation
n_windows = 2;
win_length = round(size(profile_struct.profile,1) / n_windows);
[temp_fft, f] = pwelch(profile_struct.profile, win_length,[],[], fs);
pxx = mean(temp_fft,2);
end

