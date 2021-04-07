function [new_res_coll, new_aff_pop] = chooseNeurons(resp, neuron_selection_modes, ...
        target_rates, top_neuron_number, aff_pop, loc)
%chooseNeurons selects neurons in the population based on their location or responses.



% make three masks and "or" them together at the end.
masks = zeros(length(resp.rate), 3);

for i = 1:3 %for each afferent
    if strcmp(neuron_selection_modes(i), "top")
        masks(:, i) = topNeurons(resp, top_neuron_number, aff_pop, i);
    elseif strcmp(neuron_selection_modes(i), "area")
        masks(:, i) = areaSelect(resp, loc, aff_pop, i);
    elseif strcmp(neuron_selection_modes(i), "best")
        masks(:, i) = bestNeurons(resp, aff_pop, top_neuron_number, target_rates(i), i);
    elseif strcmp(neuron_selection_modes(i), "best_area")
        masks(:, i) = bestArea(resp, loc, aff_pop, top_neuron_number, target_rates(i), i);
    else
        error("neuron selection mode not recognized")
    end

end

% 'or' together masks
final_mask  = (logical(masks(:,1)) | logical(masks(:,2))) | logical(masks(:,3));

% create new response collection and afferent population


%restricting responses
new_responses = resp.responses(final_mask);
%restricting affpop
new_aff_pop = resp.affpop;
new_aff_pop.afferents = new_aff_pop.afferents(final_mask);
new_res_coll = ResponseCollection(new_aff_pop, new_responses, resp.stimulus);
end

