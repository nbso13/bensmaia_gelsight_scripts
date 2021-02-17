function [abs_difference, trace, axis, fig_handle] = displayDirecAverage(gel_struct, direc)
%displayDirecAverage displays the average trace in a given direction for
%input gel struct. Good for 1 dimensional stimuli.
if direc == 'v'
    trace = mean(gel_struct.profile, 2);
    axis = gel_struct.y_axis;
elseif direc == 'h'
    trace = mean(gel_struct.profile, 1);
    axis = gel_struct.x_axis;
end
    
fig_handle = figure;
if ~(length(axis) == length(trace))
    error("Plot lengths not equal. Try changing direction.");
end
trace = trace - min(trace);
plot(axis, trace);
title("Average cross section")
xlabel("mm")
ylabel("mm")
abs_difference = max(trace);
end

