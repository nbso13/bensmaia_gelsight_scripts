function prof_struct_out = truncateProfile(profile_struct, stan_dev, bottom_flag)
% Truncates input heightmap between designated range
% tHeightMap = h_truncate(obj, 'r'/'p', [min max]) or h_truncate(obj, 'r'/'p')
if isfield(profile_struct, "truncated")
    warning("profile already truncated!")
end
if nargin == 1
    mean_h = mean(profile_struct.profile(:));
    std_h = std(profile_struct.profile(:));
    range = [mean_h - std_h*3, mean_h + std_h*3];
    disp(['Auto truncating points below ', num2str(round(range(1),2)), ' and above ', num2str(round(range(2),2)), ' mm']);
else
    mean_h = mean(profile_struct.profile(:));
    std_h = std(profile_struct.profile(:));
    if bottom_flag
        range = [mean_h - std_h*stan_dev, max(profile_struct.profile(:))];
    else    
        range = [mean_h - std_h*stan_dev, mean_h + std_h*stan_dev];
    end
    disp(['Auto truncating points below ', num2str(round(range(1),2)), ' and above ', num2str(round(range(2),2)), ' mm']);
end

mask = profile_struct.profile < range(2) & profile_struct.profile > range(1);
total_interped = size(mask, 1)*size(mask, 2) - sum(mask(:));
disp(strcat("Total outliers re-interpolated during truncation: " + num2str(total_interped)));
pre_interp = profile_struct.profile;
pre_interp(mask == 0) = NaN;
if bottom_flag
    tHeightMap = fillmissing(pre_interp,'constant', min(pre_interp(:)));
else
    tHeightMap = fillmissing(pre_interp,'linear','EndValues','nearest');
    tHeightMap = inpaint_nans(tHeightMap);
end
prof_struct_out = profile_struct;
prof_struct_out.profile = tHeightMap- min(tHeightMap(:));
end
