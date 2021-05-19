function [dist_min,offset] = spkd_slide(st_cur,st_cmp,q,max_slide,step_size)
% dist_min = spkd_slide(st_cur,st_cmp,max_slide,step_size)

if nargin<5
    step_size = 0.001;
end

if nargin<4
    max_slide = 0.05;
end

offsets = -max_slide:step_size:max_slide;

dists = zeros(1,length(offsets));
for t=1:length(offsets)
    d_tmp = spkdl([st_cur; st_cmp+offsets(t)],[1 length(st_cur)+1],[length(st_cur) length(st_cur)+length(st_cmp)],q);
    dists(t) = d_tmp(2);
end

[dist_min,ind] = min(dists);
offset = offsets(ind);
