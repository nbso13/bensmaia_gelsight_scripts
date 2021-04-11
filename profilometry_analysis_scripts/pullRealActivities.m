function [activities] = pullRealActivities(rates, names, good_neurons, neuron_identities, texture_nums, speed, excludeNeurons)
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

%index into rates matrix according to speed and texture numbers
ind = speed/40;
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

