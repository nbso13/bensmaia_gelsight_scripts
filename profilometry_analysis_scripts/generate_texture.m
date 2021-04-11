function [texture_struct] = generate_texture(texture_type, dot_freq, dot_height, dot_diameter, res, window_size)
%generate_texture makes a dot pattern with specified characteristics in
%microns. Makes window_size squared matrix with resolution 'res' in microns per pixel. UNITS IN are
%microns UNITS OUT ARE mm. dot freq is dots per mm.
% IF GRATING: dot frequency is actually the period.
%profile

%put everything in pixels
real_dot_diameter = dot_diameter/res; %in pixels
dot_diameter = round(real_dot_diameter);
dot_radius = round(real_dot_diameter/2); 
period = ceil((1000* 1/dot_freq)/res); %in pixels
window_size = round(window_size/res); %now in pixels
radius_sq = dot_radius^2;


accepted_texture_list = "grating, dots, probe, rand_dots, gaussian_noise, flat, unit_step, unit_impulse";

if "dots" == texture_type
    %parabola equation given three points (two dot base and height), to give height given radius.
    x = [-dot_radius; 0; dot_radius]; 
    y = [0; dot_height; 0];
    f = fit(x,y,'poly2');
    %size of mat for 1 dot
    dot = zeros([dot_diameter, dot_diameter]);
    for i = 1:dot_diameter %for all points in the dot mat
        for j = 1:dot_diameter
            x_dist = (i-dot_radius)^2; %sq distance from center x and y
            y_dist = (j-dot_radius)^2;
            rad = x_dist+y_dist; %querry radius squared
            if rad>radius_sq %if outside dot base, height is 0
                dot(i,j) = 0;
            else
                %dot height at this point is y value of parabola defined by
                %f given querry radius (square root of rad).
                dot(i,j) = f(sqrt(rad)); 
            end 
        end
    end
    %make 1 dot

    %calculate space between dots
    buffer_dist = period - dot_diameter;
    
    if buffer_dist <0
        error("Buffer distance less than 0 not allowed. Enter dot frequency so that dist is greater than the sum of two dot radii")
    end
    buffer_mat_x = zeros([size(dot,1) buffer_dist]);
    buffer_mat_y = zeros([buffer_dist size(dot,2)+buffer_dist]);

    dot = [dot buffer_mat_x];
    dot = [dot; buffer_mat_y];

    num_dots_y = ceil(window_size/size(dot,1));
    num_dots_x = ceil(window_size/size(dot,2));
    texture = repmat(dot, [num_dots_y num_dots_x]);
    texture = texture(1:window_size, 1:window_size);
    num_dots = num_dots_x*num_dots_y;
    cyl_side_area = dot_height*2*pi*dot_radius; %approximating dot surface area with area of side face of cylinder with same radius as dot base
    total_side_area = num_dots * cyl_side_area;
    nondot_area = sum(sum(texture==0));
    total_area = nondot_area + total_side_area;
    
elseif "probe" == texture_type
    probe = zeros([dot_diameter, dot_diameter]);
    
    for i = 1:dot_diameter
        for j = 1:dot_diameter
            x_dist = (i-dot_radius)^2; %sq distance from center x and y
            y_dist = (j-dot_radius)^2;
            rad = x_dist+y_dist; %querry radius squared
            if rad>radius_sq %if outside dot base, height is 0
                probe(i,j) = 0;
            else
                %dot height at this point is determined by right triangle
                %(height = a, dist to center = b, dot radius = c
                probe(i,j) = sqrt(dot_radius^2-rad); 
            end 
        end
    end
    texture = probe;
    
elseif "rand_dots" == texture_type
    
    x = [-dot_radius; 0; dot_radius]; 
    y = [0; dot_height; 0];
    f = fit(x,y,'poly2');
    %size of mat for 1 dot
    dot = zeros([dot_diameter, dot_diameter]);
    for i = 1:dot_diameter %for all points in the dot mat
        for j = 1:dot_diameter
            x_dist = (i-dot_radius)^2; %sq distance from center x and y
            y_dist = (j-dot_radius)^2;
            rad = x_dist+y_dist; %querry radius squared
            if rad>radius_sq %if outside dot base, height is 0
                dot(i,j) = 0;
            else
                %dot height at this point is y value of parabola defined by
                %f given querry radius (square root of rad).
                dot(i,j) = f(sqrt(rad)); 
            end 
        end
    end
    %make 1 dot
    
    %along one axis, number of fittable dots.
    num_possible_dots = floor(window_size/(dot_diameter));
    disp(strcat("num possible dots: ", num2str(num_possible_dots^2)));
    %num dots a random number with mean n and std sqrt(n)
    num_dots = round(randn(1)*num_possible_dots + sqrt(num_possible_dots)*num_possible_dots);
    disp(strcat("num dots: ", num2str(num_dots)));
    %get random coordinates for dots
    x_locations = randi([1 num_possible_dots], 1, num_dots);
    y_locations = randi([1 num_possible_dots], 1, num_dots);
    
    texture = zeros(window_size, window_size);
    for i=1:num_dots
        x_first = (x_locations(i)-1)*dot_height*2+1;
        x_last = x_first+dot_height*2-1;
        y_first = (y_locations(i)-1) *dot_height*2+1;
        y_last = y_first+dot_height*2-1;
        %put the dot in the right place
        texture(y_first:y_last, x_first:x_last) = dot;
    end
    
elseif "grating" == texture_type
    period = round(dot_freq*1000/res); %in  pixels
    depth = dot_height;
    num_repeats = floor(window_size/period);
    one_mod = zeros(window_size, period);
    one_mod(:,1:floor(period/2)) = depth+ one_mod(:,1:floor(period/2));
    texture = repmat(one_mod, [1 num_repeats]);
    
elseif "gaussian_noise" == texture_type
    
elseif "flat"==texture_type
    texture=ones([window_size window_size])*dot_height;
elseif "unit_step" == texture_type
    texture_a=zeros([window_size/2 window_size]);
    texture_b=ones([window_size/2 window_size])*dot_height;
    texture=[texture_a;texture_b];
elseif "unit_impulse" == texture_type
    texture_a=zeros([((window_size/2)-1) window_size]);
    impulse = ones([1 window_size])*dot_height;
    texture_b = zeros([window_size/2 window_size]);
    texture = [texture_a; impulse; texture_b];
else
    error("texture type not recognized. Choose either" + accepted_texture_list);
end

texture_struct = struct;
texture_struct.type = texture_type;
texture_struct.profile = texture./1000; %height in mms
texture_struct.x_res = res/1000; %resolution in mms per pixel
texture_struct.y_res = res/1000;
texture_struct.y_axis = (1:size(texture,1))*res./1000; %axis in mms
texture_struct.x_axis = (1:size(texture,2))*res./1000; %axis in mms

if texture_type == "dots"
    texture_struct.sa = total_area*res/(10^(6)); %in mms sq
    texture_struct.radius = dot_radius*res/1000; %in mm
    texture_struct.dot_num = num_dots;

end



