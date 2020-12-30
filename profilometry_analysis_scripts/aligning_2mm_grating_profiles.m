clear
close all
cd mat_files/



filename_gel = "201116_2mm_grating_35_gel_processed";
filename_nogel = "NoGel_grating_201019";
load(filename_gel);
load(filename_nogel);
cd ..


no_gel.profile = imrotate(no_gel.profile, 90);
temp = no_gel.x_axis;
no_gel.x_axis = no_gel.y_axis;
no_gel.y_axis = temp;
no_gel = cropProfile(no_gel, "left", 1.3, "mm");
no_gel = cropProfile(no_gel, "right", 3.4, "mm");

[gel, no_gel] = bruteCropFit(gel, no_gel);


visualizeProfile(gel)
visualizeProfile(no_gel)

cd mat_files/
save("201019_no_gel_2mm_grating", "no_gel")