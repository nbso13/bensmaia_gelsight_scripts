% unityslope        --> Draws a line of unity slope in the current axes

function h = unityslope(ax)
if(nargin<1)
    ax=gca;
end
v = axis(ax); hold(ax,'on'); 
lb = min(v);
if lb < 0; lb = 0; end;
ub = v(2); d = 0.01*(ub-lb);
x = [lb:d:ub];
h = plot(ax,x,x,':k');
return