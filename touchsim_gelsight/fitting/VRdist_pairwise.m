function [dist,distvec,f,g] = VRdist_pairwise(sp1,sp2,tc)
% [dist,distvec,f,g] = VRdist_pairwise(sp1,sp2,tc)

kernel = exp(-(0:(tc*3))/tc);
kernel = kernel/sum(kernel);

f = conv(sp1,kernel,'same');
g = conv(sp2,kernel,'same');

distvec = (f-g).^2;

%dist = sum(distvec)/(sum(sp1)+sum(sp2))*10000;
dist = sum(distvec);
