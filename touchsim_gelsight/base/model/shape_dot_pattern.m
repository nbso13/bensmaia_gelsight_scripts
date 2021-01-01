function [shape,pin_offset] = shape_dot_pattern(shift, pins_per_mm, radius, dot_spacing, ...
    dot_height, noedge, random)
% [shape,pin_offset] = shape_dot_pattern(shift,pins_per_mm, radius, dot_spacing, dot_height, noedge, random)
% creates dot_pattern shape over a circular patch, with edge cutted borders
%
% inputs : (all are optional)
%          shift: shift towards x direction (in pins, default=0)
%          pins_per_mm : number of pins per mm (default=5)
%          radius : overall contact radius in mm (default=8)
%          dot_spacing : space between dot bases in x and y directions
%          (default=0.5 mm)
%          dot_height : dot height and diameter (default=1 mm)
%          noedge : flag to set edge_cutting or not (default=true)
%          random : flag to set location of dots to be random or aligned in
%          a grid (default = false)
% output : shape: 2D pin coordinates
%          pin_offset: indentation offset for each pin


if (nargin==0)|| isempty(shift)
    shift=0;
end
if (nargin<2) || isempty(pins_per_mm)
    pins_per_mm=5;
end
if (nargin<3) || isempty(radius)
    radius=8;
end
if (nargin<4) || isempty(dot_spacing)
    dot_spacing=0.5;
end
if (nargin<5) || isempty(dot_height)
    dot_height = 1;
end
if (nargin<6) || isempty(noedge)
    noedge=true;
end
if (nargin<7)
    random=0;
end

[x,y]=meshgrid(-radius:1/pins_per_mm:radius,-radius:1/pins_per_mm:radius);
r=hypot(x,y); x(r>radius)=nan; y(r>radius)=nan;

pin_offset=zeros(size(x));
z=x*pins_per_mm+shift+1000;
pin_offset(rem(z,period)<ridge_width)=1;

spanmm=1;
if(noedge)
    pin_offset=pin_offset.*(erf(-(r-radius+(spanmm))*2/(spanmm))+1)/2*depth;
end
pin_offset(pin_offset<0)=nan;

pin_offset=pin_offset(:);
x=x(:);
y=y(:);
shape=[x(~isnan(x)),y(~isnan(x))];
pin_offset=pin_offset(~isnan(x));
end

% 
% function [texture_struct] = generate_texture(texture_type, dot_freq, dot_height, window_size)
% %generate_texture makes a dot pattern with specified characteristics in
% %microns. Makes 4000x4000 matrix with resolution 1 micron. UNITS IN are
% %microns UNITS OUT ARE mm. dot freq is dots per mm.
% %profile
% 
% dot_num = dot_freq*window_size/1000;
% % disp("dot num:");
% % disp(dot_num);
% dot_spacing = ceil(1000* 1/dot_freq);
% 
% accepted_texture_list = "dots, rand_dots, gaussian_noise, flat, unit_step, unit_impulse";
% 
% if "dots" == texture_type
%     %size of mat for 1 dot
%     dot = zeros([dot_height*2, dot_height*2]);
%     height_sq = dot_height^2;
%     for i = 1:dot_height*2
%         for j = 1:dot_height*2
%             x_dist = (i-dot_height)^2;
%             y_dist = (j-dot_height)^2;
%             if x_dist+y_dist>height_sq
%                 dot(i,j) = 0;
%             else
%                 dot(i,j) = sqrt(height_sq - x_dist - y_dist);
%             end 
%         end
%     end
%     %make 1 dot
% 
%     %calculate space between dots
%     buffer_dist = dot_spacing - 2*dot_height;
%     
%     if buffer_dist <0
%         error("Buffer distance less than 0 not allowed. Enter dot frequency so that dist is greater than the sum of two dot radii")
%     end
%     buffer_mat_x = zeros([size(dot,1) buffer_dist]);
%     buffer_mat_y = zeros([buffer_dist size(dot,2)+buffer_dist]);
% 
%     dot = [dot buffer_mat_x];
%     dot = [dot; buffer_mat_y];
% 
%     num_dots_y = ceil(window_size/size(dot,1));
%     num_dots_x = ceil(window_size/size(dot,2));
%     texture = repmat(dot, [num_dots_y num_dots_x]);
%     texture = texture(1:window_size, 1:window_size);
%     
% elseif "rand_dots" == texture_type
%     
%     dot = zeros([dot_height*2, dot_height*2]);
%     height_sq = dot_height^2;
%     for i = 1:dot_height*2
%         for j = 1:dot_height*2
%             x_dist = (i-dot_height)^2;
%             y_dist = (j-dot_height)^2;
%             if x_dist+y_dist>height_sq
%                 dot(i,j) = 0;
%             else
%                 dot(i,j) = sqrt(height_sq - x_dist - y_dist);
%             end 
%         end
%     end
%     %make 1 dot
%     
%     %along one axis, number of fittable dots.
%     num_possible_dots = floor(window_size/(dot_height*2));
%     disp(strcat("num possible dots: ", num2str(num_possible_dots^2)));
%     %num dots a random number with mean n and std sqrt(n)
%     num_dots = round(randn(1)*num_possible_dots + sqrt(num_possible_dots)*num_possible_dots);
%     disp(strcat("num dots: ", num2str(num_dots)));
%     %get random coordinates for dots
%     x_locations = randi([1 num_possible_dots], 1, num_dots);
%     y_locations = randi([1 num_possible_dots], 1, num_dots);
%     
%     texture = zeros(window_size, window_size);
%     for i=1:num_dots
%         x_first = (x_locations(i)-1)*dot_height*2+1;
%         x_last = x_first+dot_height*2-1;
%         y_first = (y_locations(i)-1) *dot_height*2+1;
%         y_last = y_first+dot_height*2-1;
%         %put the dot in the right place
%         texture(y_first:y_last, x_first:x_last) = dot;
%     end
%     
