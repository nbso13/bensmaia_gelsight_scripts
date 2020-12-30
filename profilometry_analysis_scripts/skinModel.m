function [offset_out] = skinModel(shape, pin_offset, pin_radius, samp_freq, plot_flag)
%skinModel takes in a touchsim stimulus and applies the skin mechanics
%model, giving back the skin profile for that stimulus.

% initial indentation matrix

% final indentation matrix

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


[X, Y] = meshgrid(x_axis, y_axis);
re_shaped  = reshape(pin_offset',size(X));
% plot result
if plot_flag
    figure
    subplot(1,2,1)
    surf(X, Y, re_shaped);
    title("Before skin mechanics")
end

[~,~,S1]=CircIndent2LoadProfile(pin_offset', shape, 1, pin_radius); % set pin diam to 2mm

re_shaped  = reshape(S1,size(X));
% plot result
if plot_flag
    subplot(1,2,2)
    surf(X, Y, re_shaped);
    title("After skin mechanics")
end
offset_out = S1;
end

