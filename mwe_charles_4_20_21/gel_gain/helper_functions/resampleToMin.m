function [gel_out, no_gel_out] = resampleToMin(gel, no_gel)
%resampleToMin resamples each profile down to the min resolution at each
%dimension

min_x_res = min([gel.x_res, no_gel.x_res]);
min_y_res = min([gel.y_res, no_gel.y_res]);

gel_x_axis = 0:min_x_res:gel.x_axis(end);
gel_y_axis = 0:min_y_res:gel.y_axis(end);
[gel_x_mat, gel_y_mat] = meshgrid(gel.x_axis, gel.y_axis);
gel_prof = interp2(gel_x_mat, gel_y_mat, gel.profile, gel_x_axis, gel_y_axis');

no_gel_x_axis = 0:min_x_res:no_gel.x_axis(end);
no_gel_y_axis = 0:min_y_res:no_gel.y_axis(end);
[no_gel_x_mat, no_gel_y_mat] = meshgrid(no_gel.x_axis, no_gel.y_axis);
no_gel_prof = interp2(no_gel_x_mat, no_gel_y_mat, no_gel.profile, no_gel_x_axis, no_gel_y_axis');

gel_out = gel;
gel_out.profile = gel_prof;
gel_out.x_axis = gel_x_axis;
gel_out.y_axis = gel_y_axis;
gel_out.x_res = min_x_res;
gel_out.y_res = min_y_res;

no_gel_out = no_gel;
no_gel_out.profile = no_gel_prof;
no_gel_out.x_axis = no_gel_x_axis;
no_gel_out.y_axis = no_gel_y_axis;
no_gel_out.x_res = min_x_res;
no_gel_out.y_res = min_y_res;
end

