function [prof_struct] = shape2profilometry(shape, pin_offset, pins_per_mm)
%shape2profilometry takes a touchsim shape and pin_offset and outputs a
%profilometry struct

%correct for negative coordinates in shape
minx = min(shape(:, 1), [], 'all');
shape(:, 1) = shape(:, 1) -minx;
miny = min(shape(:, 2), [], 'all');
shape(:, 2) = shape(:, 2) -miny;

%make x and y axes
y_axis = shape(1,2); %y coordinates are layed out right away
for i = 2:length(shape)
    if shape(i-1, 2) < shape(i, 2)
        y_axis = [y_axis, shape(i, 2)];
    else %if it's not increasing any more you've finished
        break
    end
end
xax = shape(:, 1); %x coordinates are first col of shape but skip size y axis each time
x_axis = xax(1:length(y_axis):length(xax));
texture = zeros([length(y_axis), length(x_axis)]);

for i = 1:length(x_axis)
    for j = 1:length(y_axis)
        texture(j, i) = pin_offset((i-1)*length(y_axis)+j);
    end
end


prof_struct = struct;
prof_struct.profile = texture;
prof_struct.x_axis = x_axis;
prof_struct.y_axis = y_axis;
prof_struct.x_res = 1/pins_per_mm; %mm per pin or sample frequency
prof_struct.y_res = prof_struct.x_res;


end

