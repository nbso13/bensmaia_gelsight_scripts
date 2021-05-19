function [mask] = randSampleAll(resp, aff_pop, index, samples)
%randSampleArea takes n neurons at random for the given afferent class that are within the
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

len_rates = length(resp.rate);
indices = 1:len_rates;
mask = zeros(1, len_rates);
affs = resp.affpop.afferents;

for i = indices %only take if 1) the right afferent 2) in x bounds and 3) in y bounds
    mask(i) = logit(i) && ~(0 == resp.rate(i));
end

total_chosen = sum(mask);
if samples > total_chosen
    warning("More samples requested than neurons in area.")
else
    rand_indices = randsample(total_chosen, samples);
    ind_count = 1;
    for i = 1:length(mask)
        if mask(i)
            if ~ismember(ind_count, rand_indices) 
                %if the id of this chosen neuron is not in the rand sample, take it out.
                mask(i) = 0;
            end
            ind_count = ind_count + 1;
        end
    end
end
mask = logical(mask);
end