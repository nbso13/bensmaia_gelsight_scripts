function shape = img2shape(img,loc,pins_per_mm,new_pins_per_mm)
% shape = img2shape(img,loc,pins_per_mm,new_pins_per_mm)

if nargin>3
    [a,b] = size(img);
    [x,y] = meshgrid(linspace(1,b,b),linspace(1,a,a));
    [X,Y] = meshgrid(linspace(1,b,b/pins_per_mm*new_pins_per_mm),linspace(1,a,a/pins_per_mm*new_pins_per_mm));
    img = interp2(x,y,img,X,Y);
end

[l1,l2] = find(~isnan(img));
shape = [l1 l2]-repmat(min([l1 l2]),length(l1),1);
shape = shape/pins_per_mm; 
for dim=1:2
    if isnan(loc(dim))
        shape(:,dim) = shape(:,dim)-(max(shape(:,dim))/2);
    else
        shape(:,dim) = shape(:,dim)+loc(dim);
    end
end
