function [shape_out, offset_out] = rotateTouchsim(shape, offset)
%rotateTouchsim flips touchsim.


old_x = shape(:,1);
% get sort indices to sort y 
[sorted_old_y, ind_vec] = sort(shape(:,2));
shape_out = [sorted_old_y, old_x(ind_vec)];
offset_out = offset(ind_vec);
end

