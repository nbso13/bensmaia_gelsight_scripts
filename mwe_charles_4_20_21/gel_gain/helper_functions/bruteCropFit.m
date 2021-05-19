function [gel_out, no_gel_out] = bruteCropFit(gel, no_gel)
%bruteCropFit gets the min in each direction and crops everyone to it.
gel_y_len = length(gel.y_axis);
no_gel_y_len = length(no_gel.y_axis);
gel_x_len = length(gel.x_axis);
no_gel_x_len = length(no_gel.x_axis);
gel_out = gel;
no_gel_out = no_gel;

if gel_y_len > no_gel_y_len
    diff_y = gel_y_len-no_gel_y_len;
    gel_out = cropProfile(gel, 'bottom', diff_y, 'px');
elseif gel_y_len < no_gel_y_len
    diff_y = no_gel_y_len-gel_y_len;
    no_gel_out = cropProfile(no_gel, 'bottom', diff_y, 'px');
else %even
    no_gel_out = no_gel;
    gel_out = gel;
end
    
if gel_x_len > no_gel_x_len
    diff_x = gel_x_len-no_gel_x_len;
    gel_out = cropProfile(gel_out, 'right', diff_x, 'px');
    
elseif no_gel_x_len > gel_x_len
    diff_x = no_gel_x_len-gel_x_len;
    no_gel_out = cropProfile(no_gel_out, 'right', diff_x, 'px');
end

