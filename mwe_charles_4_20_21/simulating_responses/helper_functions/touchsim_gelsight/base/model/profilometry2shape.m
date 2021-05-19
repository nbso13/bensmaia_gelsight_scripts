function [shape, offset] = profilometry2shape(prof_struct, pins_per_mm)
%takes in profilometry struct in mm, formats into
%touchsim shape Nx2 matrix  and pin offset as a vector.
x_axis = prof_struct.x_axis;
y_axis = prof_struct.y_axis;
x_res = prof_struct.x_res;
y_res = prof_struct.y_res;
z_res = prof_struct.z_res;
heightMap = prof_struct.profile;

%% resample for pin density (cut off, no padding)
y_size = length(y_axis); x_size = length(x_axis);
%number of pins necessary based on pins per mm and size of map in mm
x_pins = floor(pins_per_mm*x_res*x_size);
y_pins = floor(pins_per_mm*y_res*y_size);
resamp_prof = zeros(y_pins, x_pins);
resamp_factor_x = floor(x_size/x_pins);
resamp_factor_y = floor(y_size/y_pins);

%resample at pin density
for i = 1:y_pins
    y_start = (i-1)*resamp_factor_y+1;
    y_end = i*resamp_factor_y;
    for j = 1:x_pins
        x_start = (j-1)*resamp_factor_x+1;
        x_end = j*resamp_factor_x;
        resamp_prof(i,j) = mean(heightMap(y_start:y_end, x_start:x_end), 'all');
    end
end

%now resamp _ prof is downsampled to pin resolution.
[y_size, x_size] = size(resamp_prof);

N = zeros(y_size*x_size, 3);

%go through height map and make a 3 column matrix with 3D points for each
%point
for y_val = 1:y_size
    for x_val = 1:x_size
        matrix_val = (y_val-1)*x_size+x_val; 
        N(matrix_val, 1) = y_val;
        N(matrix_val, 2) = x_val;
        N(matrix_val, 3) = resamp_prof(y_val, x_val);
    end
end

%% Formatting for TouchSim
shape = N(:,1:2);
shape(:,1) = shape(:,1)-y_size/2;
shape(:,2) = shape(:,2)-x_size/2;
shape = shape./pins_per_mm;
offset = N(:,3); % in mm
end

