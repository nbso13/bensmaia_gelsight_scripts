function [mask] = bestNeurons_match(res, aff_pop, target_rates, index)
%bestNeurons_match matches the best neurons by firing rate to real recorded
%affs
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
master_indices_array = indices(logit); % separate lists of indices of target afferents
number_total_afferents = length(rate_array);
aff_indices_mask = zeros(1, number_total_afferents); % make mask to find the top master indices
ind_list = [];
for i = 1:length(target_rates)
    rate_array_compare = abs(rate_array - target_rates(i));
    [~, sorted_indices] = sort(rate_array_compare, 'ascend'); %sorted rates of afferent class
    
    for j = 1:length(sorted_indices) %find best index not already on list
        ind = sorted_indices(j);
        if ~ismember(ind, ind_list)
            break
        end
    end
    
    ind_list = [ind_list, ind];
    aff_indices_mask(ind) = 1; %set mask to one in found top index for this neuron
end

aff_indices_mask = logical(aff_indices_mask);
top_master_indices = master_indices_array(aff_indices_mask); %take according to local mask

%make overall mask
overall_top_mask = zeros(1,len_rates);
overall_top_mask(top_master_indices) = 1;
mask = logical(overall_top_mask);

end

