%% Testing touchsim skin mech control

% force profile P from treating no_gel->skin mech as skin and solving for P 
% i.e., S1*D^(-1) should be the same as
% using circ2load to get P out of no_gel profile
clear
close all
gel_constant = 1.49;


% using 1mm grating, get P from direct blocksolve.

% 1MM GRATING THIN
% filename_gel = "201116_1mm_grating_35_gel_processed";
% filename_nogel = "201021_no_gel_1mm_grating";

% 2MM GRATING THIN
filename_gel = "201116_2mm_grating_35_gel_processed";
filename_nogel = "201019_no_gel_2mm_grating";


%% Load data process data
cd ../../mat_files/
load(filename_gel);
load(filename_nogel);
cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts

gel.profile = gel.profile.*gel_constant; %scale up

if ~checkSizeMatch(gel, no_gel)
    [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
    [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
end

gel = rotateProfilometry(gel, 90);
no_gel = rotateProfilometry(no_gel, 90);

figure
visualizeProfile(gel);
figure
visualizeProfile(no_gel);

%% calculate P directly thru blocksolve
ppm = 11;
% calculate P directly through blocksolve
% calculation of P only is based on no gel profile.
plot_flag = 1;

mm_per_pin = 1/ppm;
pin_radius = mm_per_pin*0.5;

[gel_ts, no_gel_ts, skin_surface_ts, P] = TouchSimSkin(gel, no_gel, ppm, pin_radius, plot_flag);
% gel_ts is touchsim form of gel profile
% no gel ts is touchsim form of texture profile
% skin surface ts is touchsim skin profile

%% Calculate P indirectly from skin surface touchsim
% turn into touchsim format
P_comp = PfromSkinProfile(skin_surface_ts.shape, skin_surface_ts.offset, pin_radius);
P_diff = P - P_comp;

x = P(:);
y = P_comp(:);
figure;
hold on
scatter(x, y);

fit_ob = fit(x,y, 'poly1');
m = fit_ob.p1;
n = fit_ob.p2;
ax = plot(fit_ob);
legend(ax, sprintf('y=%f*x + %f',m, n))
xlabel("P blocksolve")
ylabel("P from D inverse")
title("TouchSim")

%% Calculate P indirectly from skin surface GELSIGHT
P_gel = PfromSkinProfile(gel_ts.shape, gel_ts.offset', pin_radius);
P_diff = P - P_gel;

x = P(:);
y = P_gel(:);
figure;
hold on
scatter(x, y);

fit_ob = fit(x,y, 'poly1');
m = fit_ob.p1;
n = fit_ob.p2;
ax = plot(fit_ob);
legend(ax, sprintf('y=%f*x + %f',m, n))
xlabel("P blocksolve")
ylabel("P from D inverse")
title("Gelsight")

figure;
imagesc(P_gel)
c = colorbar;
ylabel(c, 'mm');
caxis([0, max(P(:))]);
figure;

imagesc(P)
c = colorbar;
ylabel(c, 'mm');
caxis([0, max(P(:))]);