function [gel_ts, no_gel_ts, skin_surface_ts, surf_figures, P] = TouchSimSkin(gel, ...
    no_gel, pins_per_mm, pin_radius, flip_flag, plot_flag)
%touchsimskin takes in two profilometry structs and returns the gel but
%downsampled as the rest of them are, the "new gel" or the touchsim version
%model of the skin mech and the new_no_gel, downsampled to touchsim level.
%P returned is local stresses as calculated for touchsim gel.
%amplitudes is amplitude gel then amplitude ts.

if flip_flag
    no_gel = rotateProfilometry(no_gel, 90);
end

gel_surf = 0;
no_gel_surf = 0;
skin_surf = 0;
mm_per_pin = 1/pins_per_mm;
[shape, offset] = profilometry2shape(no_gel, pins_per_mm);
no_gel_ts = struct;
no_gel_ts.shape = shape;
no_gel_ts.offset = offset;
no_gel_ts.pins_per_mm = pins_per_mm;
no_gel_ts.name = no_gel.name;
no_gel_ts.gel_flag = 0; %not a gel, but the actual textured surface

[shape, offset] = profilometry2shape(gel, pins_per_mm);
gel_ts = struct;
gel_ts.shape = shape;
gel_ts.offset = offset;
gel_ts.pins_per_mm = pins_per_mm;
gel_ts.name = gel.name;
gel_ts.gel_flag = 1; %gelsight gel profile

% Subtract min offset height

gel_ts.offset = gel_ts.offset - min(gel_ts.offset);
no_gel_ts.offset = no_gel_ts.offset - min(no_gel_ts.offset);


%% touchsim operation on shape, skin mechanics on NO gel

before_after_plot_flag = 0;
[new_offset, P] = skinModel(no_gel_ts.shape, no_gel_ts.offset, pin_radius, before_after_plot_flag);
skin_surface_ts = no_gel_ts;
if flip_flag
    [flip_shape, new_offset] = rotateTouchsim(no_gel_ts.shape, new_offset);
    skin_surface_ts.shape = flip_shape;
end
skin_surface_ts.offset = new_offset;
skin_surface_ts.name = no_gel.name;
skin_surface_ts.gel_flag = 1; %gelsight gel profile



if plot_flag
    surf_figures = figure;
    subplot(1,3,1)
    surfTouchSim(gel_ts.shape, gel_ts.offset);
    title(gel_ts.name)
    subplot(1,3,2)
    surfTouchSim(no_gel_ts.shape, no_gel_ts.offset);
    title(no_gel_ts.name)
    subplot(1,3,3)
    surfTouchSim(skin_surface_ts.shape, skin_surface_ts.offset);
    title(strcat(skin_surface_ts.name, " TouchSim"));
end

%zero out offsets

skin_surface_ts.offset = skin_surface_ts.offset - min(skin_surface_ts.offset);
skin_surface_ts.area = no_gel.x_axis(end)*no_gel.y_axis(end);
no_gel_ts.offset = no_gel_ts.offset - min(no_gel_ts.offset);
no_gel_ts.area = no_gel.x_axis(end)*no_gel.y_axis(end);
gel_ts.offset = gel_ts.offset - min(gel_ts.offset);
gel_ts.area = gel.x_axis(end)*gel.y_axis(end);
end

