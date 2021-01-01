function x = denan(x,id)

if nargin<2
    x(isnan(x)) = [];
else
    x(x==id) = [];
end
