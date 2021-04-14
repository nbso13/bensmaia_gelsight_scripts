function [activities, averaged_spike_trains, space_vec] = pullRealActivities(rates, spikes, names, good_neurons, neuron_identities, texture_nums, speed, excludeNeurons)
%pullRealRates grabs real neural data, and calculated rates for each
%neuron, sorts into afferent class, calculates mean firing rate for each
%texture for each afferent, and returns a mean_rates array, where rows are
%textures and columns are afferent classes (PC, RA, SA) and entries are
%average firing rate. Next three columns are standard error on the mean of
%same afferent classes.

if texture_nums == 1
    error("for now, this doesn't work with only one texture.")
end
%determine "good neuron" indices
good_logit = zeros(39,1);
for i = 1:length(good_logit)
    if ismember(i, good_neurons)
        good_logit(i) = 1;
    end
end
good_logit = logical(good_logit);

%indices for each class
iPC = neuron_identities{1};
iRA = neuron_identities{2};
iSA = neuron_identities{3};

iPC = iPC(good_logit);
iRA = iRA(good_logit);
iSA = iSA(good_logit);

ind = speed/40;
spikes_mat = spikes{ind};
spikes_mat = squeeze(spikes_mat(texture_nums,:, :));
spikes_mat = spikes_mat(:, good_logit, :);


pcs = spikes_mat(:, iPC, :);
ras = spikes_mat(:, iRA, :);
sas = spikes_mat(:, iSA, :);
afferents = {pcs, ras, sas};

samp_period = 0.001; %seconds LATER MULT BY 80mm/s
dec_place = -1*log10(samp_period);
time_axis = 0:samp_period:1;
averaged_spike_trains = zeros(length(texture_nums), 3, length(time_axis)); % 3 being the number of aff classes (PC, RA, SA)

win_size = 11;
win = gausswin(win_size);
pad_size = (win_size-1)/2;

for i = 1:length(texture_nums)
    for j = 1:3
        aff = afferents{j};
        mean_train = zeros(1, length(time_axis));
        train_count = 0;
        for k = 1:size(aff, 2)
            for m = 1:4
                run = aff{i,k,m};
                if ~isempty(run) %if this afferent class for this texture for this neuron on this trial is not empty
                    run = round(run, dec_place);
                    run = run.*(10^dec_place) + 1; %turning spike times into indices
                    mean_train(run) = mean_train(run)+1;
                    train_count = train_count + 1;
                end
            end
        end
%         mean_train = (mean_train./train_count);
%         pad_before = mean_train(1:pad_size);
%         pad_after = mean_train(end-pad_size:end);
%         mean_train = [pad_before, mean_train, pad_after];
        mean_train = conv(mean_train, win, 'same');
        averaged_spike_trains(i, j, :) = mean_train;
    end
end

space_vec = time_axis.*speed; %now x vector is in mm.

rates_mat = rates{ind};
rates_mat = squeeze(rates_mat(texture_nums,:, :));
rates_mat = rates_mat(:, good_logit,:);

% rate_mat is three dims - texture num, neuron num, run num. Here we're
% collapsing the last two dimensions to give a cell array, each cell an
% afferent type, each cell texture_num by number of afferents.
full_rates = {reshape(rates_mat(:, iPC, :), [length(texture_nums), length(iPC(iPC))*4]),...
    reshape(rates_mat(:, iRA, :), [length(texture_nums), length(iRA(iRA))*4]),...
    reshape(rates_mat(:, iSA, :), [length(texture_nums), length(iSA(iSA))*4])};

for i = 1:length(full_rates) %for each afferent type
    aff_rates =  full_rates{i};
    if excludeNeurons
        aff_rates(aff_rates == 0) = nan;
    end
    mean_rates(:, i) = mean(aff_rates, 2, 'omitnan');
    mean_rates(:, i+3) = std(aff_rates, 0, 2, 'omitnan')/sqrt(length(~isnan(aff_rates)));
end
activities = struct;
activities.names = names;
activities.real = mean_rates;
end

