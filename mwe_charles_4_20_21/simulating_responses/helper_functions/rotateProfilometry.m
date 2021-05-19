function [profilometry] = rotateProfilometry(profilometry_in, angle)
%rotateProfilometry rotates the input struct counterclockwise
new_prof = profilometry_in;
if angle ~= 90 && angle ~= 180 && angle ~= 270
    error("caution: angle must be divisible by 90 or axes are off")
end
new_prof.profile = imrotate(profilometry_in.profile, angle);
if angle ~= 180
    new_prof.x_axis = profilometry_in.y_axis;
    new_prof.y_axis = profilometry_in.x_axis;
    new_prof.x_res = profilometry_in.y_res;
    new_prof.y_res = profilometry_in.x_res;
end

profilometry = new_prof;
end

