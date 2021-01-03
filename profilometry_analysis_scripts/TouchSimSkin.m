function [gel_ts, no_gel_ts, skin_surface_ts, P] = TouchSimSkin(gel, no_gel, pins_per_mm, pin_radius, plot_flag)
%touchsimskin takes in two profilometry structs and returns the gel but
%downsampled as the rest of them are, the "new gel" or the touchsim version
%model of the skin mech and the new_no_gel, downsampled to touchsim level.
%P returned is local stresses as calculated for touchsim gel.
mm_per_pin = 1/pins_per_mm;
[shape, offset] = profilometry2shape(no_gel, pins_per_mm);
no_gel_ts = struct;
no_gel_ts.shape = shape;
no_gel_ts.offset = offset;
no_gel_ts.pins_per_mm = pins_per_mm;
no_gel_ts.name = "texture surface";
no_gel_ts.gel_flag = 0; %not a gel, but the actual textured surface

[shape, offset] = profilometry2shape(gel, pins_per_mm);
gel_ts = struct;
gel_ts.shape = shape;
gel_ts.offset = offset;
gel_ts.pins_per_mm = pins_per_mm;
gel_ts.name = "gel surface";
gel_ts.gel_flag = 1; %gelsight gel profile


%save(strcat(filename_nogel, "_ts"), "no_gel_ts");

%% touchsim operation on shape, skin mechanics on NO gel
cd ../touchsim_gelsight/
setup_path;
cd ../profilometry_analysis_scripts/
[new_offset, P] = skinModel(no_gel_ts.shape, no_gel_ts.offset, pin_radius, mm_per_pin, plot_flag);
skin_surface_ts = no_gel_ts;
skin_surface_ts.offset = new_offset;
skin_surface_ts.name = "touchsim skin surface";
skin_surface_ts.gel_flag = 1; %gelsight gel profile


if plot_flag
    surfTouchSim(gel_ts.shape, gel_ts.offset)
    title(gel_ts.name)
    surfTouchSim(no_gel_ts.shape, no_gel_ts.offset)
    title(no_gel_ts.name)
    surfTouchSim(skin_surface_ts.shape, skin_surface_ts.offset)
    title(skin_surface_ts.name)
end
end

