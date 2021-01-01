function out=GLM_nonlinfun(param,data)

x=data-param(1);

% rectified
x(x<0)=0;

% abs value
%x=abs(x);

out=param(2)*x.^param(3);
