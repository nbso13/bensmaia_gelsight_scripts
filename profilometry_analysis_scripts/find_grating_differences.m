function [differences] = find_grating_differences(profilometry, direction, min_ranges, max_ranges, plot_flag)
%find_grating_differences is good for getting conformance values from a
%craig stimulus all in one profilometry struct. It takes in a profile,
%averages it down in the indicated direction to one average trace,
%calculates local min within min ranges and max for max ranges. ranges
%should be x by 2 in dimensions, where x is the number of points.
n = size(min_ranges, 1);
if ~(n == size(max_ranges, 1))
    error("must ask for the same number of min and maxes.");
end

differences = zeros(n,1);
factor = 0;
if direction == 'v'
    trace = mean(profilometry.profile, 1);
    factor = profilometry.x_res;
    if plot_flag
        figure
        plot(profilometry.x_axis, trace)
        title("Average vertical trace")
        xlabel("x loc")
        ylabel("height (mm)")
    end
elseif direction =='h'
    trace = mean(profilometry.profile, 2);
    factor = profilometry.y_res;
    if plot_flag
        figure
        plot(profilometry.y_axis, trace)
        title("Average vertical trace")
        xlabel("y loc")
        ylabel("height (mm)")
    end
else
    error("direction unrecognized. use h for horizontal (averaging horizontally) or v for vertical (averaging down)");
end

min_ranges = min_ranges./factor;
max_ranges = max_ranges./factor;
mins = zeros(n,1);
maxs = mins;
for i = 1:n
    mins(i) = min(trace(ceil(min_ranges(i,1)):ceil(min_ranges(i,2))));
    maxs(i) = max(trace(ceil(max_ranges(i,1)):ceil(max_ranges(i,2))));
end


min_ranges = min_ranges.*factor;
max_ranges = max_ranges.*factor;
mid_min = (min_ranges(:,1)+min_ranges(:,2))./2;
mid_max = (max_ranges(:,1)+max_ranges(:,2))./2;

if plot_flag
    hold on;
    scatter(mid_min, mins);
    scatter(mid_max, maxs);
end

differences = maxs-mins;
end

