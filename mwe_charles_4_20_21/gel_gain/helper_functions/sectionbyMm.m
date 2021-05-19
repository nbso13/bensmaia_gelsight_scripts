function [output_arr] = sectionbyMm(gel, xy_vec)
%sectionbyMm takes in two x,y coordinate pairs, one lower left hand one
%upper right hand, in mm, and a gel or no_gel profilometry struct. Returns
%the profile between them.
if xy_vec(1)>xy_vec(2) || xy_vec(3)> xy_vec(4)
    error("check that first x and y coordinates are bigger than second. It goes [x1,x2,y1,y2]");
end

x1_px = floor(xy_vec(1)/gel.x_res)+1; %otherwise will round to 0 and you can't index with 0 in matlab
x2_px = floor(xy_vec(2)/gel.x_res);
y1_px = floor(xy_vec(3)/gel.y_res)+1;
y2_px = floor(xy_vec(4)/gel.y_res);
if x2_px > length(gel.x_axis)
    x2_px = length(gel.x_axis);
end
if y2_px > length(gel.y_axis)
    y2_px = length(gel.y_axis);
end
output_arr = gel.profile(y1_px:y2_px, x1_px:x2_px);
end

