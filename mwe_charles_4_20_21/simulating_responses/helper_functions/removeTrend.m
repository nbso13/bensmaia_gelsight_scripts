function [prof_struct] = removeTrend(prof_struct_in)
%removeTrend fits a plane and subtracts from profile

x_axis = prof_struct_in.x_axis;
y_axis = prof_struct_in.y_axis;
y_size = length(y_axis);
x_size = length(x_axis);
new_window = prof_struct_in.profile;
N = zeros(y_size*x_size, 3);
count = 1;
reverseStr = '';
%go through height map and make a 3 column matrix with 3D points for each
%point
for x_val = 1:x_size
    x_amount = x_axis(x_val);
    msg = sprintf('    fitting line %d of %d', count, x_size);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    count = count + 1;
    for y_val = 1:y_size
        matrix_val = (x_val-1)*y_size+y_val;
        N(matrix_val, 1) = x_amount;
        N(matrix_val, 2) = y_axis(y_val);
        N(matrix_val, 3) = new_window(y_val, x_val);
    end
end
fprintf('\n');
%calculate plane and fit
[plane, ~] = fit([N(:,1), N(:,2)], N(:,3), 'poly11');
params = coeffvalues(plane);

%% Filling in New Plane and Subtracting
grid_form_plane = zeros(size(new_window));
%fill a new grid height map with the trend plane
count = 1;
reverseStr = '';
for x_val = 1:x_size
    msg = sprintf('    fitting line %d of %d', count, x_size);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    count = count + 1;
    for y_val = 1:y_size
        grid_form_plane(y_val, x_val) = x_axis(x_val)*params(2) + y_axis(y_val)*params(3);
    end
end

deplaned = new_window - grid_form_plane; % Subtract out
deplaned = deplaned - min(min(deplaned)); %min should still be zero
prof_struct = prof_struct_in;
prof_struct.profile = deplaned;


end

