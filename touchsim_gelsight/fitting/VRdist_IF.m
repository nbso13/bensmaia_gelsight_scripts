function [dist,distvec,f,g] = VRdist_IF(params,Stim,spvec,tc)
% [dist,distvec,f,g] = VRdist_MN(params,Stim,spvec,tc)

[~,predvec] = predict_IF(params,Stim);

kernel = exp(-(0:(tc*3))/tc);
kernel = kernel/sum(kernel);

f = conv(spvec,kernel,'same');
g = conv(predvec,kernel,'same');

distvec = (f-reshape(g,[],1)).^2;
dist = sum(distvec)/sum(spvec)*100;
