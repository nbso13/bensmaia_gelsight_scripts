function [pred,predvec] = predict_IF(params,stim,noisy)
% [pred,predvec] = predict_IF(params,stim,noisy)

if nargin<3
    aff.noisy = false;
else
    aff.noisy = noisy;
end

aff.parameters = params;
prop_struct.stat_comp = stim(:,1);
prop_struct.dyn_comp = stim(:,2);
prop_struct.sampling_frequency = 20000;
spikes = IF_neuron(aff,prop_struct);

edges = 0:1/5000:length(stim)/20000;
predvec = histc(spikes,edges)';
predvec(end) = [];
pred = find(predvec);
