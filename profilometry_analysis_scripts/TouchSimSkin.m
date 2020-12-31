function [gel, new_gel, no_gel] = TouchSimSkin(gel, no_gel, pins_per_mm, plot_flag)
%touchsimskin takes in two profilometry structs and returns the gel but
%downsampled as the rest of them are, the "new gel" or the touchsim version
%model of the skin mech and the new_no_gel, downsampled to touchsim level.

mm_per_pin = 1/pins_per_mm;
pin_radius = mm_per_pin*0.5;
[shape, offset] = profilometry2shape(no_gel, pins_per_mm);
no_gel_ts = struct;
no_gel_ts.shape = shape;
no_gel_ts.offset = offset;
no_gel_ts.pins_per_mm = pins_per_mm;

[shape, offset] = profilometry2shape(gel, pins_per_mm);
gel_ts = struct;
gel_ts.shape = shape;
gel_ts.offset = offset;
gel_ts.pins_per_mm = pins_per_mm;

%save(strcat(filename_nogel, "_ts"), "no_gel_ts");

%% touchsim operation on shape, skin mechanics on NO gel
cd ../../touchsim/
setup_path;
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts/
new_offset = skinModel(no_gel_ts.shape, no_gel_ts.offset, pin_radius, mm_per_pin, plot_flag);
new_gel_ts = no_gel_ts;
new_gel_ts.offset = new_offset;

if plot_flag
    surfTouchSim(gel_ts.shape, gel_ts.offset)
    title("Gel")
    surfTouchSim(no_gel_ts.shape, no_gel_ts.offset)
    title("No Gel")
    surfTouchSim(new_gel_ts.shape, new_gel_ts.offset)
    title("TouchSim")
    disp("now skin modeling!")
end

%% back to profilometry
new_gel = shape2profilometry(new_gel_ts.shape, new_gel_ts.offset, new_gel_ts.pins_per_mm);
no_gel = shape2profilometry(no_gel_ts.shape, no_gel_ts.offset, no_gel_ts.pins_per_mm);
gel = shape2profilometry(gel_ts.shape, gel_ts.offset, gel_ts.pins_per_mm);

end

