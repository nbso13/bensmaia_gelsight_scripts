function [] = tsPlotSpikes(spikes, time_length, good_neurons, names, neuron_identities, texture_num, speed)
%tsPlotSpikes wrapps plot_spikes in calculating color matrix

% pull spiketimes from neurons by afferent

good_logit = zeros(39,1);
for i = 1:length(good_logit)
    if ismember(i, good_neurons)
        good_logit(i) = 1;
    end
end
good_logit = logical(good_logit);

iPC = neuron_identities{1};
iRA = neuron_identities{2};
iSA = neuron_identities{3};

iPC = iPC(good_logit);
iRA = iRA(good_logit);
iSA = iSA(good_logit);
name_str = names{texture_num};

ind = speed/40;
spikes_mat = spikes{ind};
spikes_mat = squeeze(spikes_mat(texture_num,:, :));
spikes_mat = spikes_mat(good_logit, :);

pcs = spikes_mat(iPC, :); num_pcs = sum(iPC);
ras = spikes_mat(iRA, :); num_ras = sum(iRA);
sas = spikes_mat(iSA, :); num_sas = sum(iSA);

%taking only the first trial...
pcs = pcs(:,1);
ras = ras(:,1);
sas = sas(:,1);
spike_times = [sas;ras;pcs];

% taking out spikes after time length is up
for i = 1:length(spike_times)
    train = spike_times{i};
    train = train(train<time_length);
    spike_times{i} = train;
end


pc_col =  [255 127 0]/255;
ra_col =  [30 120 180]/255;
sa_col = [50 160 40]/255;

col_mat = vertcat(repmat(sa_col, num_sas, 1), repmat(ra_col, num_ras, 1), repmat(pc_col, num_pcs, 1));

plot_spikes(spike_times, 'color', col_mat)
end

