function [P_comp] = PfromSkinProfile(shape,offset, pin_radius)
%PfromProfile takes in touchsim parameters that is a skin profile, either
%from GelSight or touchsim, and calculates and returns the local stresses.
%As calculated in CircLoad.

cd ../touchsim_gelsight/
setup_path;
cd ../profilometry_analysis_scripts/

y_axis = shape(1,2); %y coordinates are layed out right away
for i = 2:length(shape)
    if shape(i-1, 2) < shape(i, 2)
        y_axis = [y_axis, shape(i, 2)];
    else %if it's not increasing any more you've finished
        break
    end
end
xax = shape(:, 1); %x coordinates are first col of shape but skip size y axis each time
x_axis = xax(1:length(y_axis):length(xax));
[X, ~] = meshgrid(x_axis, y_axis);
gel_flag = 1; %this is already a skin profile
%calculate P backward using touchsim derived skin profile.
[P_comp, ~, ~] = CircIndent2LoadProfile(offset, shape, 1, pin_radius, gel_flag);
P_comp = reshape(P_comp ,size(X));
end

