function [] = visTexture(shape,pin_offset, pin_per_mm)
%visTexture: Given shape and pin offset, visualizes texture using color map
% 
mino = min(shape, [], 'all');

shape = shape.*10 -mino*10+1;
texture = zeros(ceil(range(shape)));
y_axis = 1:size(texture, 1);
x_axis = 1:size(texture,2);
y_axis = y_axis./pin_per_mm;
x_axis = x_axis./pin_per_mm;

for i = 1:length(shape)
    texture(floor(shape(i,1)), floor(shape(i,2))) = pin_offset(i);
end


figure
title_str = "Visualized Texture";
imagesc(x_axis, y_axis, texture)
ax = gca;
ax.YDir = 'normal';
c = colorbar;
ylabel(c, 'mm');
caxis([0, max(texture, [], 'all')]);
title(title_str);
xlabel('mm'); ylabel('mm');
    
end

