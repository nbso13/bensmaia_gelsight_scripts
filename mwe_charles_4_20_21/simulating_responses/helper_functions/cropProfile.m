function [gel_struct_out] = cropProfile(gel_struct_in, crop_direction, pixel_value, unit)
%cropProfile crops the profilometry profile and modifies struct
%appropriately. takes in the struct, crop direction string (either left
%right top or bottom), and the number of pixels to cut off (crop
%value).
gel_struct_out = gel_struct_in;
    
if crop_direction == "top"
    if unit == "mm"
            pixel_value = floor(pixel_value/gel_struct_in.y_res); %how many pixels to chop off, assuming y_res is mm's per pixel
    end
    gel_struct_out.profile = gel_struct_in.profile(pixel_value:end, :);
    gel_struct_out.y_axis = gel_struct_in.y_axis(pixel_value:end);
elseif crop_direction == "bottom"
    if unit == "mm"
            pixel_value = floor(pixel_value/gel_struct_in.y_res); %how many pixels to chop off, assuming y_res is mm's per pixel
    end
    gel_struct_out.profile = gel_struct_in.profile(1:end-pixel_value, :);
    gel_struct_out.y_axis = gel_struct_in.y_axis(1:end-pixel_value);
elseif crop_direction == "right"
    if unit == "mm"
            pixel_value = floor(pixel_value/gel_struct_in.x_res); %how many pixels to chop off, assuming y_res is mm's per pixel
    end
    gel_struct_out.profile = gel_struct_in.profile(:, 1:end-pixel_value);
    gel_struct_out.x_axis = gel_struct_in.x_axis(1:end-pixel_value);
elseif crop_direction == "left"
    if unit == "mm"
            pixel_value = floor(pixel_value/gel_struct_in.x_res); %how many pixels to chop off, assuming y_res is mm's per pixel
    end
    gel_struct_out.profile = gel_struct_in.profile(:, pixel_value:end);
    gel_struct_out.x_axis = gel_struct_in.x_axis(pixel_value:end);
else
    error("Crop direction not accurately specified. Should be top bottom left or right.")
end
gel_struct_out.x_axis = gel_struct_out.x_axis - gel_struct_out.x_axis(1);
gel_struct_out.y_axis = gel_struct_out.y_axis - gel_struct_out.y_axis(1);
end

