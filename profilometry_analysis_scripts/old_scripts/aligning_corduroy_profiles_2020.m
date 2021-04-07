clear
close all
cd mat_files/



filename_gel = "201118_corduroy_35_gel_processed";
filename_nogel = "201118_corduroy_no_gel_processed";
load(filename_gel);
load(filename_nogel);
cd ..




gel = cropProfile(gel, "top", 0.4, "mm");
no_gel = cropProfile(no_gel, "bottom", 0.6, "mm");


visualizeProfile(gel)
visualizeProfile(no_gel)


cd mat_files/
save("201118_corduroy_no_gel_trimmed", "no_gel")
save("201118_corduroy_35_gel_trimmed", "gel")
cd ..