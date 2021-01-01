% GLM neuron model
% inputs should be sampled at 5000 Hz

function spikes = GLM_neuron(aff,stat_comp,dyn_comp)

sf=5000; % sampling frequency is 5000 Hz

time_fac = sf/5000;
stimi = reshape(stat_comp,[],1);
dstimi = reshape(dyn_comp,[],1);
ddstimi = [diff(dstimi);0]*time_fac;

stim_expanded = rectify([stimi -stimi dstimi -dstimi ddstimi -ddstimi]);
if(aff.iPC)
    stim_expanded=stim_expanded(:,3:6);
else
    stim_expanded=stim_expanded(:,1:4);
end

plin=reshape(aff.parameters(1:end-3),200,4);
pnonlin=aff.parameters(end-2:end);

% linear pred
%lincoef=interp1(linspace(1,200,40),plin,1:200); % interpolate to sampling rate
linpred=sum(conv2(stim_expanded,plin,'same'),2);
% non lin part
nonlinpred=GLM_nonlinfun(pnonlin,linpred);
% poisson
spfit=nonlinpred>rand(size(nonlinpred));
spikes=find(spfit)/sf;
