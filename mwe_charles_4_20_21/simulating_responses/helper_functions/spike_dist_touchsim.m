function [dist_mat] = spike_dist_touchsim(spikes, spike_train_sim, ...
    time_length, good_neurons, neuron_identities, texture_num, speed)
%spike_dist_touchsim computes the spike distance metric between two input
%spike trains, one real and the other simulated.

qs = 1./[0.001, 0.005, 0.0010];
max_slide = 0.05;
step_size = 0.001;

good_logit = zeros(39,1);
for i = 1:length(good_logit)
    if ismember(i, good_neurons)
        good_logit(i) = 1;
    end
end
good_logit = logical(good_logit);
iPC = neuron_identities{1}; iRA = neuron_identities{2}; iSA = neuron_identities{3};
iPC = iPC(good_logit); iRA = iRA(good_logit); iSA = iSA(good_logit);

ind = speed/40;
spikes_mat = spikes{ind};
spikes_mat = squeeze(spikes_mat(texture_num,:, 1))'; %taking only the first trial
spikes_mat = spikes_mat(good_logit, :);
pcs = spikes_mat(iPC); ras = spikes_mat(iRA); sas = spikes_mat(iSA); 


real_spikes = {pcs, ras, sas};
sim_pcs = {}; sim_ras = {}; sim_sas = {};
pc_count = 1; ra_count = 1; sa_count = 1;

for i = 1:length(spike_train_sim)
    if spike_train_sim(i).afferent.class == "PC"
        sim_pcs{pc_count} = spike_train_sim(i).spikes;
        pc_count = pc_count + 1;
    elseif spike_train_sim(i).afferent.class == "RA"
        sim_ras{ra_count} = spike_train_sim(i).spikes;
        ra_count = ra_count + 1;
    elseif spike_train_sim(i).afferent.class == "SA1"
        sim_sas{sa_count} = spike_train_sim(i).spikes;
        sa_count = sa_count + 1;
    else
        error("aff type not recognized.")
    end
end

sim_spikes = {sim_pcs, sim_ras, sim_sas};
dist_mat = {};

for i = 1:3
    real_spike_trains = real_spikes{i};
    sim_spike_trains = sim_spikes{i};
    
    len_aff_real = size(real_spike_trains, 1);
    len_aff_sim = size(sim_spike_trains, 2);
    aff_dist_mat = NaN(len_aff_real, len_aff_sim);
    spike_counts = aff_dist_mat;
    
    for j = 1:len_aff_real
        for k = 1:len_aff_sim
            spikes_real = real_spike_trains{j}; spikes_real = spikes_real(spikes_real<time_length);
            spikes_sim = sim_spike_trains{k}; spikes_sim = spikes_sim(spikes_sim<time_length);
            
            [dist_min,~] = spkd_slide(spikes_real, spikes_sim', qs(i), max_slide, step_size);
            aff_dist_mat(j, k) = dist_min;
            spike_counts(j, k) = length(spikes_real) + length(spikes_sim);
        end
    end
    
    dist_mat{i} = aff_dist_mat./spike_counts;
end
          

end

