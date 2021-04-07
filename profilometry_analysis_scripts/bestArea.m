function [mask] = bestArea(resp, loc, aff_pop, top_neurons, target_rate, index)
%bestArea takes n closest neurons to known firing rate for the given afferent class that are within the
%designated location (where the texture is in contact with the skin)

if index == 1
    logit = aff_pop.iPC;
elseif index == 2
    logit = aff_pop.iRA;
elseif index == 3
    logit = aff_pop.iSA1;
else
    error("index does not refer to known afferent type")
end

min_x = loc(1);
max_x = loc(2);
min_y = loc(3);
max_y = loc(4);

len_rates = length(resp.rate);
indices = 1:len_rates;
mask = zeros(1, len_rates);
affs = resp.affpop.afferents;

for i = indices %only take if 1) the right afferent 2) in x bounds and 3) in y bounds
    mask(i) = logit(i) && (min_x <= affs(i).location(1)) && ...
        (affs(i).location(1) < max_x) && ...  
        (min_y <= affs(i).location(2)) && ...
        (affs(i).location(2) < max_y);
end
    
mask = logical(mask);

% determined neurons in area
% now doing best match neurons. 

mask = mask & ~(resp.rate == 0)';
rate_array = resp.rate(mask); % separate list of rates
rate_array = abs(rate_array - target_rate);
master_indices_array = indices(mask); % separate lists of indices of target afferent
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

