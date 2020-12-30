function [] = compare_lines(gel, no_gel, gel_val, no_gel_val, direction)
%compare_lines plots an overlapping plot of the two cross sections (used to
%compare gain values
gel_constant = 1;
figure;
subplot(2,2,1);
visualizeProfile(gel) %show gel with a line through it

gel_coord = floor(gel_val/gel.x_res)+1; %pixels = mm* pixels/mm
no_gel_coord = floor(no_gel_val/no_gel.x_res)+1; %pixels = mm* pixels/mm

if direction == 'v' %if vertical
    xline(gca, gel_val); %write line on gel
    no_gel_line = no_gel.profile(:, no_gel_coord); %get no_gel_vertical profile
    gel_line = gel.profile(:, gel_coord); % gel gel_vertical profile
    gel_ax =  gel.y_axis;
elseif direction == 'h'
    yline(gca, gel_val);
    no_gel_line = no_gel.profile(no_gel_coord, :); %get no_gel hor prof
    gel_line = gel.profile(gel_coord, :); % get gel hor prof
    gel_ax = gel.x_axis;
end

subplot(2,2,2);
visualizeProfile(no_gel) % now for no_gel
if direction == 'v'
    xline(gca, no_gel_val);
    no_gel_ax = no_gel.y_axis;
elseif direction == 'h'
    yline(gca, no_gel_val);
    no_gel_ax = no_gel.x_axis;
end

min_val_gel = min(gel_line);
min_val_no_gel = min(no_gel_line);
no_gel_line = no_gel_line - min_val_no_gel;
gel_line = gel_line - min_val_gel;

subplot(2,2,3)
hold on
plot(no_gel_ax, no_gel_line);
plot(gel_ax, gel_line*gel_constant);
end

