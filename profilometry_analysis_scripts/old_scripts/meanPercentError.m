function [mean_perc_error] = meanPercentError(gel, human)
%calculates mean percent error between gel and human points
errors = zeros(length(human), 1);
for i = 1:length(human)
    error_val = sqrt((human(i)-gel(i))*(human(i)-gel(i)));
    errors(i) = 100*error_val/human(i);
end
mean_perc_error = mean(errors);
end

