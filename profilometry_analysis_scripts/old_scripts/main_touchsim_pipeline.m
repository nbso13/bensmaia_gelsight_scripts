%% Scanning Textures
% Nick Ornstein

close all
clear
cd ../../touchsim
setup_path;

%load('force_profile.mat');
%% Plot firing rates for each indentation
% load('firing_rate.mat')
%
% figure
%
% x = firing_rate(1,:);
% errhigh = firing_rate(3,:);
% errlow = errhigh;
% mean_fr = firing_rate(2,:);
% bar(x, mean_fr)
%
% hold on
%
% er = errorbar(x, mean_fr, errlow, errhigh);
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% title('Mean FR For Simulated Population Over Indentations')
% xlabel('indentation (mm)')
% ylabel('FR (Hz)')

%% Afferent
a = affpop_hand('D2d', 0.4, 'SA1');
figure
title("Location of Receptors")
plot(a)


%% Generated Texture
% we want 2mm spatial period
% so dot freq is 1/2 per mm
% dot_freq = 0.5;
% dot_height = 740; % microns. From Weber et al 2013
% dot_diameter = 500; %
% window_size = 8000; %dots per mm and microns respectively
% 
% pins_per_mm = 13;
% res = 15; %10 micron resolution
% dot_pattern = generate_texture("dots", dot_freq, dot_height, dot_diameter, res, window_size);
% %dot_pattern = generate_texture("probe", dot_freq, dot_height, 5000, res, window_size);
% %visualizeProfile(dot_pattern)
% [shape, pin_offset] = profilometry2shape(dot_pattern, pins_per_mm);
% visTexture(shape,pin_offset, pins_per_mm)


%5mm period grating
% period  = 3; %mm
% height = 740;% microns
% resolution = 20; % microns per pixel
% window_size = 10000;%10 mms
% pins_per_mm = 7;
% grating = generate_texture("grating", period, height, 0, resolution, window_size);
% visualizeProfile(grating)
% [shape, pin_offset] = profilometry2shape(grating, pins_per_mm);
% visTexture(shape,pin_offset, pins_per_mm)


% shift = 0;
% period = 50; %10 pins, or 1mm
% depth = 0.74;
% noedge = 1;
% pins_per_mm = 5;
% [shape,pin_offset] = shape_square_grating(shift, pins_per_mm, period, [], [], depth, noedge);
% visTexture(shape,pin_offset, pins_per_mm)
%% Gelsight
is_gel = 0; %is a gel being used?
load("no_gel_grating_ts.mat");
prof = no_gel;
if ~isfield(prof, 'pins_per_mm')
    prof.pins_per_mm = 10;
end
shape = prof.shape;
pin_offset = prof.offset;

%% profilometry shape
pins_per_mm = prof.pins_per_mm;
if is_gel
    amp = max(pin_offset);
else
    amp = 0.8; %mm
end
speed = 80; %mm/s.
len = 0.5; % s
loc = [0 0];
samp_freq = 200; % hz
ramp_len = 0.2;
%s = stim_indent_shape_new(shape,stim_ramp(amp, len, loc, samp_freq, ramp_len), pin_offset);
s = stim_scan_shape(shape, pin_offset, pins_per_mm, len, samp_freq, amp, speed);
figure
plot(s)
    

%% Titrating amplitudes, getting responses for different amplitudes
rates = {};
rate_counter = 1;
amps = [ 0.4, 0.5, 0.6];
%forces = zeros(size(amps));
for amplitude = 1:length(amps)
    %disp(strcat("For ", num2str(amps(amplitude)), " mm indentation"));
    %determine amplitude -
    % in weber et al, 2013, 0.5 N force used on drum. Assume contact area 100mm sq,
    % pressure then would be 5kPa. Assume pressure increases linearly by 2546
    % Pa/mm indentation, you get
    amp = amps(amplitude);%5/2.546; %mm
    speed = 80; %mm/s.
    len = 1; % s
    loc = [0 0];
    samp_freq = 3000; % hz
    ramp_len = 0.2;
    %s = stim_indent_shape_new(shape,stim_ramp(amp, len, loc, samp_freq, ramp_len), pin_offset);
    s = stim_scan_shape(shape, pin_offset, pins_per_mm, len, samp_freq, amp, speed);
    figure
    plot(s)
    
    %total_force = sum(s.profile(10,:));
    %forces(amplitude) = total_force;
    %disp(strcat("Total force: ", num2str(total_force), " N"));
    %     calculate response
    r = a.response(s);
    pin_num = size(r.stimulus.location,1);

    
    %take out neurons that fire less than 2 spikes per second
    r_new = excludeNeurons(r, 2);
    %r_new = r;
    rates{rate_counter} = r_new.rate;
    figure
    plot(r_new)
    strang = strcat("SA activity at ", num2str  (amp)," mm indent.fig");
    title(strang)
    h=gcf;
    %savefig(h,strang)
    rate_counter = rate_counter+1;
end

%% get strain matrix
%time_val = 0.5;
%strain = reconstructStrain(r,time_val);
mean_fr = zeros(1,length(rates));
sd_fr = zeros(1,length(rates));
for i = 1:length(rates)
    mean_fr(i) = mean(rates{i});
    sd_fr(i) = std(rates{i});
end

%% plot firing rate
firing_rate = [amps; mean_fr; sd_fr;];
figure

x = firing_rate(1,:);
errhigh = firing_rate(3,:);
errlow = errhigh;
mean_fr = firing_rate(2,:);
bar(x, mean_fr)

hold on

er = errorbar(x, mean_fr, errlow, errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
title('Mean FR For Simulated Population Over Indentations')
xlabel('indentation (mm)')
ylabel('FR (Hz)')

