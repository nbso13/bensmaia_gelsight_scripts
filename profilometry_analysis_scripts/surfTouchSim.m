function [surf_fig] = surfTouchSim(shape, pin_offset)
%surfTouchSim takes in a shape and pin offset pair and plots a surface
mini = min(pin_offset);
pin_offset = pin_offset - mini;
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


[X, Y] = meshgrid(x_axis, y_axis);
re_shaped  = reshape(pin_offset',size(X));
% plot result
surf_fig = figure;
surf(X, Y, re_shaped, 'edgecolor', 'none');
colormap parula;
end

