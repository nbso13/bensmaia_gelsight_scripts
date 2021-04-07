%% converting between shape and profilometry
clear
close all

cd mat_files/
filename_gel = "201116_2mm_grating_35_gel_processed";
filename_nogel = "201019_no_gel_2mm_grating";
load(filename_gel);
load(filename_nogel);
cd ..

gel_constant = 1.5;
gel.profile = gel.profile.*gel_constant;
%% profilometry to shape

pins_per_mm = 18;
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
surfTouchSim(gel_ts.shape, gel_ts.offset)
title("Gel")

%save(strcat(filename_nogel, "_ts"), "no_gel_ts");
disp("now skin modeling!")
%% touchsim operation on shape, skin mechanics on NO gel
cd ../touchsim/
setup_path;

plot_flag = 1;
new_offset = skinModel(no_gel_ts.shape, no_gel_ts.offset, pin_radius, mm_per_pin, plot_flag);
new_no_gel_ts = no_gel_ts;
new_no_gel_ts.offset = new_offset;


%visTexture(new_no_gel_ts.shape, new_no_gel_ts.offset, new_no_gel_ts.pins_per_mm);
%visTexture(no_gel_ts.shape, no_gel_ts.offset, no_gel_ts.pins_per_mm);
%% back to profilometry

no_gel_mod = shape2profilometry(new_no_gel_ts.shape, new_no_gel_ts.offset, new_no_gel_ts.pins_per_mm);
no_gel = shape2profilometry(no_gel_ts.shape, no_gel_ts.offset, no_gel_ts.pins_per_mm);
gel = shape2profilometry(gel_ts.shape, gel_ts.offset, gel_ts.pins_per_mm);
cd ../current_scripts
