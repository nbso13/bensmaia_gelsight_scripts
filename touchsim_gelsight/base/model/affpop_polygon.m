function a = affpop_polygon(polygon,density,class,varargin)
% a = affpop_polygon(polygon,density,class,varargin)
% Places receptors with a given density in cm^2 inside a polyon, whose
% coordinates are supplied in mm.

if nargin<3 || isempty(class)
    afftype = {'SA1','RA','PC'};
else
    afftype = {class};
end

if nargin<2 || isempty(density)
    density = 25;
end

if length(density)<length(afftype)
    density = repmat(density,1,length(afftype));
end

if nargin<1 || isempty(polygon)
    polygon = [0 0; 10 0; 10 10; 0 10];
end

a=AfferentPopulation();
for tt=1:length(afftype)
    locs = sample_random_shape(density(tt)/100,polygon);
    a = a.add_afferents(afftype{tt},locs,varargin{:});
end
