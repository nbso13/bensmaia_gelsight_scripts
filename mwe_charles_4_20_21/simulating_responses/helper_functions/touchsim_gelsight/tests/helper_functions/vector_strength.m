function [r,ang] = vector_strength(angles,rev)
% [r,ang] = vector_strength(angles,rev)
% calculates vector strength for given phase angles

if nargin<2
    rev = 1;
end

if rev==1
    r = sqrt(sum(cos(angles))^2 + sum(sin(angles))^2)/length(angles);
    ang = circ_mean(angles);
else
    bimod_angles = mod(angles,pi)*2;
    r = sqrt(sum(cos(bimod_angles))^2 + sum(sin(bimod_angles))^2)/length(bimod_angles);
    ang = circ_mean(bimod_angles)/2;
    ang(2) = ang(1) + pi;
end
