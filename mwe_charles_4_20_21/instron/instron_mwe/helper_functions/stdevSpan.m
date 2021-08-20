function [x2, inBetween] = stdevSpan(mean_run, std_run)
%stdevSpan takes in a mean_run and std struct and gives back the x axis and the
%inBetween coordinates to shade around the mean, the std.
lower = mean_run.force-std_run.force;
upper = mean_run.force+std_run.force;
x2 = [mean_run.indentation, fliplr(mean_run.indentation)];
inBetween = [lower, fliplr(upper)];
inBetween(isnan(inBetween)) = 0;
end

