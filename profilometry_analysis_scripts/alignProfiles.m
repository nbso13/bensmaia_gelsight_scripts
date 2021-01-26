function [prof1_new, prof2_new] = alignProfiles(prof1,prof2, plot_flag)
%alignProfiles takes profile 1 and profile 2 and uses cross correlation to
%align them, cropping the first input to the second.
if ~(prof1.x_res ==prof2.x_res)
    [prof1, prof2] = resampleToMin(prof1, prof2); %resamples to the min resolution
end



c = normxcorr2(prof2.profile, prof1.profile);
if plot_flag
    figure, surf(c), shading flat
end

[~, imax] = max(abs(c(:)));
[ypeak, xpeak] = ind2sub(size(c),imax(1));
corr_offset = [xpeak, ypeak];

xoffset = corr_offset(1);
yoffset = corr_offset(2);

xbegin = round(xoffset+1);
xend   = round(xoffset+ size(prof1.profile, 2));
ybegin = round(yoffset+1);
yend   = round(yoffset+size(prof1_profile, 1));

if xoffset > 0 % if prof1 starts after prof2
    [prof1_new] = cropProfile(prof1, "left", -xbegin, "px"); %cut off
else
    [prof2_new] = cropProfile(prof2, "left", xbegin, "px"); %cut off
end
if yoffset > 0 % if prof1 starts after prof2
    [prof1_new] = cropProfile(prof1, "top", ybegin, "px"); %cut off
else
    [prof2_new] = cropProfile(prof2, "top", -ybegin, "px"); %cut off
end

end

