function prof_struct_out = truncateProfile(profile_struct, stan_dev, range)
% Truncates input heightmap between designated range
% tHeightMap = h_truncate(obj, 'r'/'p', [min max]) or h_truncate(obj, 'r'/'p')
if isfield(profile_struct, "truncated")
    error("profile already truncated!")
end
if nargin == 1
    mean_h = mean(profile_struct.profile(:));
    std_h = std(profile_struct.profile(:));
    range = [mean_h - std_h*3, mean_h + std_h*3];
    disp(['Auto truncating points below ', num2str(round(range(1),2)), ' and above ', num2str(round(range(2),2)), ' mm']);
end

if nargin == 2
    mean_h = mean(profile_struct.profile(:));
    std_h = std(profile_struct.profile(:));
    range = [mean_h - std_h*stan_dev, mean_h + std_h*stan_dev];
    disp(['Auto truncating points below ', num2str(round(range(1),2)), ' and above ', num2str(round(range(2),2)), ' mm']);
end

mask = profile_struct.profile < range(2) & profile_struct.profile > range(1);
total_interped = size(mask, 1)*size(mask, 2) - sum(mask(:));
disp(strcat("Total outliers re-interpolated during truncation: " + num2str(total_interped)));
pre_interp = profile_struct.profile;
pre_interp(mask == 0) = NaN;
tHeightMap = fillmissing(pre_interp,'linear','EndValues','nearest');
tHeightMap = tHeightMap - min(tHeightMap(:));
prof_struct_out = profile_struct;
prof_struct_out.profile = tHeightMap;

end
