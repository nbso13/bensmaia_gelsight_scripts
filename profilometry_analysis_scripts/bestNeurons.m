function [mask] = bestNeurons(res, aff_pop, top_neurons, target_rate, index)
%topNeurons takes only the responses from top n neurons for the given
%afferent clsas
len_rates = length(res.rate);
indices = 1:len_rates;

if index == 1
    logit = aff_pop.iPC;
elseif index == 2
    logit = aff_pop.iRA;
elseif index == 3
    logit = aff_pop.iSA1;
else
    error("index does not refer to known afferent type")
end

%subract target rate from rates

logit = logit & ~(res.rate == 0)';

rate_array = res.rate(logit); % separate list of rates
rate_array = abs(rate_array - target_rate);
master_indices_array = indices(logit); % separate lists of indices of target afferent
[~, sorted_indices] = sort(rate_array, 'ascend'); %sorted rates of afferent class
top_indices = sorted_indices(1:top_neurons); %top indices
number_total_afferents = length(sorted_indices);
aff_indices_mask = zeros(1, number_total_afferents); % make mask to find the top master indices
aff_indices_mask(top_indices) = 1; %set mask to one in found top indices
aff_indices_mask = logical(aff_indices_mask);
top_master_indices = master_indices_array(aff_indices_mask); %take according to local mask

%make overall mask
overall_top_mask = zeros(1,len_rates);
overall_top_mask(top_master_indices) = 1;
mask = logical(overall_top_mask);
end

