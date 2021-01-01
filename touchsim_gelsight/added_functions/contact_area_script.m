%% Pulling Contact Area
%% Scanning Textures
% Nick Ornstein
close all
clear
%% Square grating
% shift = 0;
% period = 20;
% width = 10;
% contact_radius = 4;
% depth = 1;
% noedge = 1;
% pins_per_mm = 10;
% [shape,pin_offset] = shape_square_grating(shift, pins_per_mm, period, width, contact_radius, depth, noedge);
% visTexture(shape,pin_offset, pins_per_mm)
%% Dot Pattern
%we want 2mm spatial period
% so dot freq is 1/2 per mm
dot_freq = 0.5;
dot_height = 740; % microns. From Weber et al 2013
dot_diameter = 500; %
window_size = 8000; %dots per mm and microns respectively
area = window_size*window_size/1000000;
disp(strcat("Area: ", num2str(area), " mm sq"))

pins_per_mm = 14;
res = 15; %10 micron resolution
dot_pattern = generate_texture("dots", dot_freq, dot_height, dot_diameter, res, window_size);
%dot_pattern = generate_texture("probe", dot_freq, dot_height, 5000, res, window_size);
%visualizeProfile(dot_pattern)
[shape, pin_offset] = profilometry2shape(dot_pattern, pins_per_mm);
visTexture(shape,pin_offset, pins_per_mm)


%% Gelsight
%load("no_gel_ts");
%shape = no_gel_ts.shape;
%pin_offset = no_gel_ts.offset;
%% profilometry shape
rates = {};
rate_counter = 1;
amps = [ 0.4, 0.5, 0.6];
areas = zeros(size(amps));

for amplitude = 1:length(amps)
    disp(strcat("For ", num2str(amps(amplitude)), " mm indentation"));
    %determine amplitude -
    % in weber et al, 2013, 0.5 N force used on drum. Assume contact area 100mm sq,
    % pressure then would be 5kPa. Assume pressure increases linearly by 2546
    % Pa/mm indentation, you get
    amp = amps(amplitude);%5/2.546; %mm
    len = 1; % s
    loc = [0 0];
    samp_freq = 3000; % hz
    ramp_len = 0.2;
    s = stim_indent_shape_new(shape,stim_ramp(amp, len, loc, samp_freq, ramp_len), pin_offset);
    figure
    plot(s)
    total_force = sum(s.profile(10,:));
    forces(amplitude) = total_force;
    disp(strcat("Total force: ", num2str(total_force), " N"));
end

%% plot force/ind curve and pressure/ind curve

% 
figure;
plot(amps, forces)
title('Force/Indentation Curve Dots 441 mm sq')
xlabel('Indentation (mm)')
ylabel('Force (N)')

plotflag=1;
[ind, sa] = surfaceArea4Ind(dot_pattern, plotflag);
sa_at_points = interp1(ind, sa, amps); %in mm sq
pressures = forces./(sa_at_points./1000); %in kPa
pressures(1) = 0;
pressures(isnan(pressures)) = forces(isnan(pressures))./(441/1000); % i.e, whole surface area. in kPa
figure
plot(amps, pressures)
title('Pressure/Indentation Curve Dots 441 mm sq')
xlabel('Indentation (mm)')
ylabel('Pressure (kPa)')
